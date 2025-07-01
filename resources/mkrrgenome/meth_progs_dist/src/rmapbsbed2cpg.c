/* rmapbsbed2cpg: to take rmapbs output and generate a list of
  individual CpG positions therefrom.

Process requires reading entire human genome, chromosome by chromosome
and complete reads (fasta & fastq format in PROG_VERSION 2.0).
So have linked list of reads, for each tile...

*/

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>

#include "bas_fns.h"
#include "wlu_fns.h"
#include "sqfl_fns.h"
#include "sqmat_fns.h"
#include "rmapbsbed2cpg.h"
#include "rbs_fns.h"
#include "mrg_fns.h"

/*#define PROG_VERSION 0.00 *//* Nov-2010 */
/* corrected out-by-one errors rmapbs: 1-Apr-2011 */
/* #define PROG_VERSION 1.0 */
/* modify to support HiSeq V2 & V3 headers */
/* #define PROG_VERSION 2.0 */
/* corrections for bsmap out-by-1 and truncated read info: Oct-2011 */
/* #define PROG_VERSION 2.1 */
/* incorporate code to check read mapping positions wrt RR genome positions: Oct-2011 */
/* #define PROG_VERSION 2.2 */
/* Allow for bsmap v1.2 output */
/* #define PROG_VERSION 2.3 */
/* don't run with no chromosomes: Oct-2011 */
#define PROG_VERSION 2.31

/* global debug & out style variables, for simplicity of access */
RBC_DBG debuglevel;

int err_msg(char *fmt,
            ...)
/* write user error message.  Return 0 for err return status */
{
va_list args;

va_start(args,fmt);
(void) vfprintf(stderr,fmt,args);
va_end(args);
return(0);
}

void err_msg_die(char *fmt,
                 ...)
/* write user error message then exit with error status */
{
va_list args;

va_start(args,fmt);
(void) vfprintf(stderr,fmt,args);
va_end(args);
exit(1);
}

void say_usage(FILE *fl,
               char *pnam)
{
fprintf(fl,"%s (v%.2f): generate list of CpG positions from rmapbs BED file\n",
          pnam,PROG_VERSION);
fputs("Options:\n",fl);
fputs("     -r <readfile> reads from Fasta fmt. <readfile>\n",fl);
fputs("     -R <readfile> reads from fastq <readfile>\n",fl);
fputs("     -v HiSeq V3 flowcell, V2 chemistry (def V2 GAII flowcell)\n",fl);
fputs("     -V HiSeq V3 flowcell, V3 chemistry (def V2 GAII flowcell)\n",fl);
fputs("     -b <bedfile> read match information from <bedfile> (rmapbs bed format)\n",fl);
fputs("     -B <bsmapfile> read match information from <bsmapfile> (bsmap V1.02 output format)\n",fl);
fputs("     -A <bsmapfile> read match information from <bsmapfile> (bsmap V1.2 output format)\n",fl);
fputs("     -p <posfile> read RR sections (mkrrgenome -m output) from <posfile>, count hits/misses\n",fl);
fputs("     -g <genomehead> dir and file string to locate genomic seq files by adding n.fa\n",fl);
}

RBC_READ_ELT *rbc_appndelt(RBC_READ_ELT **lstrt,
                           int xpx,
                           int ypx,
                           char *sq)
/* create and append a new element to *lstrt,
the seq is assumed to pre-exist and
is not created here. start, stop & sense values are
not set here, merely initialised. Return address of
new element */
{
RBC_READ_ELT *prev, *end_ptr;

if (lstrt != NULL)
  {
  prev = end_ptr = *lstrt;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nxtrdelt;
    }
  end_ptr = (RBC_READ_ELT *) getmemory(sizeof(RBC_READ_ELT),"Read element");
  end_ptr->nxtrdelt = NULL;
  end_ptr->xpix = xpx;
  end_ptr->ypix = ypx;
  end_ptr->pstart = end_ptr->pstop = 0;
  end_ptr->readfwd = 1;
  end_ptr->chrno = Chr_unk;
  end_ptr->readseq = sq;
  if (*lstrt == NULL)
    {
    *lstrt = end_ptr;
    end_ptr->prvrdelt = NULL;
    }
  else
    {
    prev->nxtrdelt = end_ptr;
    end_ptr->prvrdelt = prev;
    }
  return(end_ptr);
  }
else
  return(NULL);
}

void rbc_delrdelt(RBC_READ_ELT *ep,
                  RBC_READ_ELT **lstrt,
                  int clrseq)
/* delete ep from list *lstrt, if
clrseq then free that storage also */
{
RBC_READ_ELT *pt;

if (ep != NULL)
  {
  if ((pt = ep->prvrdelt) == NULL)
    *lstrt = ep->nxtrdelt;
  else
    pt->nxtrdelt = ep->nxtrdelt;
  if ((pt = ep->nxtrdelt) != NULL)
    pt->prvrdelt = ep->prvrdelt;
  if (clrseq && (ep->readseq != NULL))
    memfree(ep->readseq);
  memfree(ep);
  }
}

void rbc_clrallrdelts(RBC_READ_ELT **lstrt,
                      int clrseq)
  /* iteratively delete all of lstrt */
{
while (*lstrt != NULL)
  rbc_delrdelt(*lstrt,lstrt,clrseq);
}

int rbc_cntrdelts(RBC_READ_ELT *clst)
  /* recursively read list elements */
{
if (clst == NULL)
  return(0);
else
  return(rbc_cntrdelts(clst->nxtrdelt) + 1);
}

int rbc_cntmodrdelts(RBC_READ_ELT *clst)
  /* recursively count modified read list elements */
{
if (clst == NULL)
  return(0);
else
  if ((clst->chrno != Chr_unk) && (clst->pstart > 0) &&
       (clst->pstop > 0))
    return(rbc_cntmodrdelts(clst->nxtrdelt) + 1);
  else
    return(rbc_cntmodrdelts(clst->nxtrdelt));
}

RBC_READ_ELT *cpg_lastelt(RBC_READ_ELT *clst)
  /* iterate thru clst, returning
last element, NULL if none */
{
RBC_READ_ELT *ep;

if ((ep = clst) == NULL)
  return(NULL);
else
  {
  while (ep->nxtrdelt != NULL)
    ep = ep->nxtrdelt;
  return(ep);
  }
}

RBC_READ_ELT *rbc_elt4xypix(RBC_READ_ELT *rlst,
                            int xp,
                            int yp)
/* iterate thru rlst, returning the first element which
matches xp,yp */
{
RBC_READ_ELT *rp;

rp = rlst;
while (rp != NULL)
  if ((rp->xpix == xp) && (rp->ypix == yp))
    return(rp);
  else
    rp = rp->nxtrdelt;
return(NULL);
}

int rbc_readfqread(FILE *sfl,
                   char *hdrbuf,
                   int hbuflen,
                   char *sqbuf,
                   char *qulbuf,
                   int buflen)
/* assume sfl is at start of @ header line, read the header then
get the sequence, returning its length.  Skip quality lines, leaving
sfl on next @ header.  Write header to hdrbuf, seq to sqbuf &
qual to qulbuf if nonNULL */
{
int slen;

(void) rbc_readlntobuf(sfl,hdrbuf,hbuflen);
slen = rbc_readlntobuf(sfl,sqbuf,buflen);
(void) rbc_readlntobuf(sfl,NULL,0);
(void) rbc_readlntobuf(sfl,qulbuf,buflen);
return(slen);
}

int rbc_readnstorreads(RBC_RUNPARS *rpp,
                       SQFL_STRCT *sqf,
                       RBC_READ_ELT *rdlsts[])
/* rdlsts is the address of an array of MAXTILENO 
individual lists.  use the header line of sqf (fasta file,
already open)
to add the new read to the appropriate tile list. scan
the first 10 reads to check the read length, then
rewind the file and do the reads.
return the no of reads */
{
char *sqbuf;
int sno;
int tno;
int xpx;
int ypx;
int newlen;
RBC_READ_ELT *rdlstends[MAXTILENO+1];
char *hbuf;

rpp->srcreadlen = 0;
newlen = 1;
for (sno = 1; ((sno <= 10000) && (newlen > 0)); sno++)
  switch (rpp->readflform)
    {
    case RBC_readflfm_fastq:
      rpp->srcreadlen = imax((newlen = rbc_readfqread(sqf->sfl,NULL,0,NULL,NULL,0)),
                              rpp->srcreadlen);
      break;
    case RBC_readflfm_fasta:
    default:
      rpp->srcreadlen = imax((newlen = readsrcsq(sqf,NULL)),rpp->srcreadlen);
      break;
    }
if (rpp->srcreadlen <= 0)
  return(0);
else
  {
  for (tno = 0; tno < MAXTILENO; tno++)
    rdlstends[tno] = rdlsts[tno];
  sqbuf = (char *) getmemory(rpp->srcreadlen+1,"Sq buffer");
  sno = 0;
  switch (rpp->readflform)
    {
    case RBC_readflfm_fastq:
      rewind(sqf->sfl);
      hbuf = (char *) getmemory(128,"fastq hdr buf");
      while (rbc_readfqread(sqf->sfl,hbuf,128,sqbuf,NULL,rpp->srcreadlen) > 0)
        if ((tno = rbc_remaptileno(rpp->fcversion,rbc_hdrstr2tile(rpp,hbuf,&xpx,&ypx))) > 0)
          {
          sno++;
          rdlstends[tno-1] = rbc_appndelt(&rdlstends[tno-1],xpx,ypx,bas_strdup(sqbuf));
          if (rdlsts[tno-1] == NULL)
            rdlsts[tno-1] = rdlstends[tno-1];
          }
      break;
    case RBC_readflfm_fasta:
    default:
      sqfl_rewind(sqf);
      while ((newlen = readsrcsq(sqf,sqbuf)) > 0)
        {
        if ((tno = rbc_remaptileno(rpp->fcversion,rbc_hdrstr2tile(rpp,sqf->seqnam,&xpx,&ypx))) > 0)
          {
          sno++;
          rdlstends[tno-1] = rbc_appndelt(&rdlstends[tno-1],xpx,ypx,bas_strdup(sqbuf));
          if (rdlsts[tno-1] == NULL)
            rdlsts[tno-1] = rdlstends[tno-1];
          }
        }
      break;
    }
  return(sno);
  }
}

int rbc_readsrcfl(RBC_RUNPARS *rpp,
                  RBC_READ_ELT *rdlsts[],
                  FILE *sfl)
/* use fscanf to read successive lines from sfl.
for each, try to identify the header in rdlsts
and modify the values therein in accordance with
the sfl read.
Note: that bsmap may truncate reads, so need
to allow for variant positions/lengths for them.
This requires re-doing the read */
{
char chrstr[5];
int startpos;
int stoppos;
char hdrstr[65];
int something;
char sensestr[5];
int matcnt;
int tilno;
int xpx;
int ypx;
RBC_READ_ELT *rep;
int scnt;
char *readsq;
char rqual[5];
int qval;
char mmstr[129];
char *posstr;
int readsqlen;
char *qulstr;
char *othersq;
int ival;

matcnt = 0;
readsq = NULL;
switch (rpp->srcstyle)
  {
  case RBC_src_bsmap102:
    readsq = (char *) getmemory(rpp->srcreadlen+1,"bsmap sq read");
/* position string can be very long, esp for run of C's... */
    posstr = (char *) getmemory((rpp->srcreadlen*10),"pos string");
    while ((scnt = fscanf(sfl,"%s %s %s %s %d %s %d %s %s",&hdrstr[0],readsq,
                            &rqual[0],&chrstr[0],&startpos,&sensestr[0],
                            &qval,&mmstr[0],posstr)) != EOF)
      if ((scnt == 9) &&
           ((tilno = rbc_remaptileno(rpp->fcversion,rbc_hdrstr2tile(rpp,&hdrstr[0],&xpx,&ypx))) > 0)
           && ((rep = rbc_elt4xypix(rdlsts[tilno-1],xpx,ypx)) != NULL))
        {
        rep->readfwd = sensestr[0] == '+';
        rep->pstart = startpos;
        readsqlen = strlen(readsq);
        rep->pstop = startpos + readsqlen - 1;
        rep->chrno = rbc_str2chrno(&chrstr[0]);
        if (readsqlen != strlen(rep->readseq))
          {
          memfree(rep->readseq);
          rep->readseq = bas_strdup(readsq);
          }
        matcnt++;
        }
    memfree(readsq);
    memfree(posstr);
    break;
  case RBC_src_bsmap12:  /* has quality line included */
    readsq = (char *) getmemory(rpp->srcreadlen+1,"bsmap sq read");
/* position string can be very long, esp for run of C's... */
    qulstr = (char *) getmemory(rpp->srcreadlen+1,"bsmap 1.2 qual");
    othersq = (char *) getmemory(rpp->srcreadlen+1,"bsmap 1.2 ??");
    while ((scnt = fscanf(sfl,"%s %s %s %s %s %d %s %d %s %d %s",&hdrstr[0],readsq,
                            qulstr,&rqual[0],&chrstr[0],&startpos,&sensestr[0],
                            &qval,othersq,&ival,&mmstr[0])) != EOF)
      if ((scnt == 11) &&
           ((tilno = rbc_remaptileno(rpp->fcversion,rbc_hdrstr2tile(rpp,&hdrstr[0],&xpx,&ypx))) > 0)
           && ((rep = rbc_elt4xypix(rdlsts[tilno-1],xpx,ypx)) != NULL))
        {
        rep->readfwd = sensestr[0] == '+';
        rep->pstart = startpos;
        readsqlen = strlen(readsq);
        rep->pstop = startpos + readsqlen - 1;
        rep->chrno = rbc_str2chrno(&chrstr[0]);
        if (readsqlen != strlen(rep->readseq))
          {
          memfree(rep->readseq);
          rep->readseq = bas_strdup(readsq);
          }
        matcnt++;
        }
    memfree(readsq);
    memfree(qulstr);
    memfree(othersq);
    break;
  case RBC_src_rmapbs:
  default:
    bzero(&chrstr[0],5);
    while ((scnt = fscanf(sfl,"%s %d %d %s %d %s",&chrstr[0],
                            &startpos,&stoppos,&hdrstr[0],
                            &something,&sensestr[0])) != EOF)
      if ((scnt == 6) &&
           ((tilno = rbc_remaptileno(rpp->fcversion,rbc_hdrstr2tile(rpp,&hdrstr[0],&xpx,&ypx))) > 0)
           && ((rep = rbc_elt4xypix(rdlsts[tilno-1],xpx,ypx)) != NULL))
        {
        rep->pstart = startpos;
        rep->pstop = stoppos;
        rep->readfwd = sensestr[0] == '+';
        rep->chrno = rbc_str2chrno(&chrstr[0]);
        matcnt++;
        }
    break;
  }
return(matcnt);
}

int rbc_getgenomesqs(char *hdrstr,
                     char *seqsarr[],
                     int seqlens[],
                     int chrmax)
/* use hdrstr to create a series of file names, one
for each chromosome, open each as a Fasta sequence
file and create a buffer for each.  Return the
number of chromosomes processed */
{
char *sqfnam;
int chno;
int nblen;
SQFL_STRCT *chsqfl;
int rcnt;

rcnt = 0;
sqfnam = (char *) getmemory((nblen = strlen(hdrstr) + 16),"Sq file name buf");
for (chno = 1; chno <= chrmax; chno++)
  {
  snprintf(sqfnam,nblen-1,"%s%s.fa",hdrstr,rbc_chrno2str((RBC_CHRNO) chno,1));
  if ((chsqfl = sqfl_opnsqstrct(sqfnam,SFMT_fasta,"r")) != NULL)
    {
    seqlens[chno-1] = readsrcsq(chsqfl,NULL);
    if (debuglevel > RBC_dbg_on)
      {
      fprintf(stdout,"%s: %d res,",sqfnam,seqlens[chno-1]);
      fflush(stdout);
      }
    seqsarr[chno-1] = (char *) getmemory(seqlens[chno-1]+1,"Chr Buff");
    sqfl_rewind(chsqfl);
    (void) readsrcsq(chsqfl,seqsarr[chno-1]);
    sqfl_clssqstrct(chsqfl);
    if (debuglevel > RBC_dbg_on)
      {
      fprintf(stdout,"'%.10s...'\n",seqsarr[chno-1]);
      fflush(stdout);
      }
    rcnt++;
    }
  else
    {
    err_msg("Can't open chromosome file %s\n",sqfnam);
    seqlens[chno-1] = 0;
    }
  }
return(rcnt);
}

void rbc_rptchrout(FILE *ofl,
                   char c,
                   int ccnt)
/* put c out ccnt times to ofl */
{
while (ccnt-- > 0)
  fputc(c,ofl);
}

char *rbc_rrreadform2str(RBC_RRREAD_FORM rrdfm)
{
switch (rrdfm)
  {
  case RBC_rrread_onjoin:
    return("on join");
    break;
  case RBC_rrread_5p:
    return("5' end");
    break;
  case RBC_rrread_3p:
    return("3' end");
    break;
  case RBC_rrread_internal:
    return("within RR fragment");
    break;
  case RBC_rrread_nomap:
  default:
    return("unmapped");
    break;
  }
}

void rbc_prthdr4readelt(FILE *ofl,
                        int lno,
                        int tno,
                        RBC_READ_ELT *rlp)
/* generate a hdr type line at ofl for lane No lno, 
tileno tno & *rlp */
{
if (rlp != NULL)
  fprintf(ofl,"s%d_%d_%d_%d",lno,tno,rlp->xpix,rlp->ypix);
}

int rbc_methmismatchcnt(char *readsq,
                        char *genmsq,
                        int slen)
/* compare each residue in readsq with genmsq, returning
number of positions which mismatch, allowing for
C->T conversions */
{
char *rp;
char *gp;
char rbas;
int mmcnt;

mmcnt = 0;
rp = readsq;
gp = genmsq;
while (slen > 0)
  {
  if (toupper(rbas = *rp) == 'T')
    rbas = 'Y';
  if (!ssd_basmatch(rbas,*gp,BAS_iub))
    mmcnt++;
  rp++;
  gp++;
  slen--;
  }
return(mmcnt);
}

void rbc_nprtstr(FILE *fl,
                 char *str,
                 int nchrs)
/* put nchrs of str to fl,  any non-printable chars replaced by '*' */
{
char *sp;

sp = str;
while (nchrs > 0)
  {
  if (isprint(*sp))
    fputc(*sp,fl);
  else
    fputc('*',fl);
  sp++;
  nchrs--;
  }
}

void rbc_chkbedvsgenomseq(RBC_RUNPARS *rpp,
                          FILE *ofl,
                          int tileno,
                          RBC_READ_ELT *rdlst,
                          char *chrseqs[],
                          int chrlens[],
                          int maxcno)
/* taking each read that has valid, modified positions,
compare the read & genomic sequence to find methylated
and unmethylated CpGs. */
{
RBC_READ_ELT *rlp;
int rdlen;
char *lclchrbuf;
int lcblen;
char tchar;
char *lop;
char *hip;
int chrcnt;
int bcnt;
int sqpos;
int mmcnt;

rlp = rdlst;
lcblen = 0;
lclchrbuf = NULL;
while (rlp != NULL)
  {
  if ((rlp->chrno != Chr_unk) && (rlp->pstart > 0)
        && (rlp->pstop > 0))
    {
    rdlen = strlen(rlp->readseq);
    if (rdlen > lcblen)
      {
      if (lclchrbuf != NULL)
        memfree(lclchrbuf);
      lclchrbuf = (char *) getmemory(rdlen + 2,"lcl chr sq buf");
      lcblen = rdlen + 1;
      }
    if (debuglevel > RBC_dbg_on)
      {
      fprintf(ofl,"Chr%s: ",rbc_chrno2str(rlp->chrno,1));
      rbc_prthdr4readelt(ofl,1,tileno,rlp);
      fprintf(ofl,"\t%c\n",(rlp->readfwd?'+':'-'));
      rbc_rptchrout(ofl,' ',13);
      fprintf(ofl,"%.*s\n",rdlen,rlp->readseq);
      if (rpp->srcstyle == RBC_src_rmapbs)
        strncpy(lclchrbuf,(chrseqs[rlp->chrno-1] + rlp->pstart),rdlen);
      else
        strncpy(lclchrbuf,(chrseqs[rlp->chrno-1] + rlp->pstart - 1),rdlen);
      if (!rlp->readfwd)
        complmnt_seq(lclchrbuf,rdlen,BAS_exact);
      bcnt = 13;
      lop = rlp->readseq;
      hip = lclchrbuf;
      chrcnt = rdlen;
      while (chrcnt-- > 0)
        if (toupper(*lop++) == toupper(*hip++))
          bcnt++;
        else
          {
          rbc_rptchrout(ofl,' ',bcnt);
          bcnt = 0;
          fputc('*',ofl);
          }
      fputc('\n',ofl);
      if (rlp->readfwd)
        {
        sqpos = rlp->pstart;
        if (rpp->srcstyle == RBC_src_rmapbs)
          fprintf(ofl,"%10d:  %.*s :%d\n",rlp->pstart,rdlen,lclchrbuf,
                    (rlp->pstop-1));
        else
          fprintf(ofl,"%10d:  %.*s :%d\n",rlp->pstart,rdlen,lclchrbuf,
                    rlp->pstop);
        }
      else
        {
        sqpos = rlp->pstop - 1;
        if (rpp->srcstyle == RBC_src_rmapbs)
          fprintf(ofl,"%10d:  %.*s :%d\n",(rlp->pstop-1),rdlen,lclchrbuf,
                    rlp->pstart);
        else
          fprintf(ofl,"%10d:  %.*s :%d\n",rlp->pstop,rdlen,lclchrbuf,
                    rlp->pstart);
        }
      lop = rlp->readseq;
      hip = lclchrbuf;
      chrcnt = rdlen;
      while (chrcnt-- > 0)
        {
        if (toupper(*lop) != toupper(*hip))
          fprintf(ofl,"%c%c\t%d\n",toupper(*hip),toupper(*lop),sqpos);
        lop++;
        hip++;
        if (rlp->readfwd)
          sqpos++;
        else
          sqpos--;
        }
      rbc_rptchrout(ofl,'-',10);
      fputc('\n',ofl);
      }
    switch (rpp->ostyle)
      {
      case RBC_out_cpg:
      default:
        if (rpp->srcstyle == RBC_src_rmapbs) 
          if (rlp->readfwd)
            strncpy(lclchrbuf,(chrseqs[rlp->chrno-1] + rlp->pstart),rdlen+1);  /* need 1 more for CpG */
          else
            strncpy(lclchrbuf,(chrseqs[rlp->chrno-1] + rlp->pstart -1),rdlen+1);
        else
          strncpy(lclchrbuf,(chrseqs[rlp->chrno-1] + rlp->pstart - 1),rdlen+1);  /* need 1 more for CpG */
       *(lclchrbuf + rdlen + 1) = '\0';
       if (!rlp->readfwd)
          {
          lop = lclchrbuf;
          hip = lclchrbuf + rdlen-1;
          while (lop <= hip)
            {
            tchar = *lop;
            *lop = ssd_bascmplmnt(*hip,BAS_exact);
            *hip = ssd_bascmplmnt(tchar,BAS_exact);
            lop++;
            hip--;
            }
          }
        lop = rlp->readseq;
        hip = lclchrbuf;
        chrcnt = rdlen;
        if (rlp->readfwd)
          sqpos = rlp->pstart + 1;
        else
          sqpos = rlp->pstop - 1;
        if ((rpp->sqcheck >= 0) &&
             ((mmcnt = rbc_methmismatchcnt(rlp->readseq,lclchrbuf,rdlen)) > rpp->sqcheck))
          {
          fprintf(stderr,"SQ/Read %d mismatch ",mmcnt);
          rbc_prthdr4readelt(stderr,1,tileno,rlp);
          fprintf(stderr," %s mapped Chr %s at %d..%d\n",(rlp->readfwd?"+":"-"),
                    rbc_chrno2str(rlp->chrno,1),rlp->pstart,rlp->pstop);
          rbc_nprtstr(stderr,rlp->readseq,rdlen);
          fputc('\n',stderr);
          rbc_nprtstr(stderr,lclchrbuf,rdlen);
          fputc('\n',stderr);
          fputs("------\n",stderr);
          }
        while ((chrcnt-- > 0) && (*hip != '\0'))
          {
          if ((toupper(*hip) == 'C') && (toupper(*(hip+1)) == 'G'))  /* is CpG */
            {
            fprintf(ofl,"%s\t",rbc_chrno2str(rlp->chrno,1));
            if (debuglevel > RBC_dbg_none)
              {
              rbc_prthdr4readelt(ofl,1,tileno,rlp);
              fputc('\t',ofl);
              }
            if (rpp->srcstyle== RBC_src_rmapbs)
              fprintf(ofl,"%d\t%c\n",sqpos,((toupper(*lop)=='C')?'+':'-'));
            else
              if (rlp->readfwd)
                fprintf(ofl,"%d\t%c\n",sqpos-1,((toupper(*lop)=='C')?'+':'-'));
              else
                fprintf(ofl,"%d\t%c\n",sqpos+1,((toupper(*lop)=='C')?'+':'-'));
            }
          lop++;
          hip++;
          if (rlp->readfwd)
            sqpos++;
          else
            sqpos--;
          }
        break;
      }
    }
  rlp = rlp->nxtrdelt;
  }
}

int rbc_readmkrrgenposns(FILE *pfl,
                         MRG_REGN_ELT *rrposns[],
                         RBC_CHRNO chrlimit)
/* read a series of lines of form
Chrno\t<start>..<stop>\t(<len> bp) CpG: <cpgcnt>\n
from pfl, store relevant stuff in rrposns for appropriate
chromosome. return total No positions read */
{
MRG_REGN_ELT *rrposends[ChrY+1];
RBC_CHRNO cno;
int scnt;
char pbuf[33];
char *scp;
char dbuf1[17];
char dbuf2[17];
char dbuf3[5];
int cpgcnt;
char chrstr[5];

scnt = 0;
for (cno = Chr_unk; cno <= chrlimit; cno++)
  rrposends[cno] = NULL;
while (fscanf(pfl,"%s\t%s %s %s %s %d",&chrstr[0],&pbuf[0],&dbuf1[0],
                &dbuf2[0],&dbuf3[0],&cpgcnt) != EOF)
  {
  cno = rbc_str2chrno(&chrstr[0]);
  scp = &pbuf[0];
  while ((*scp != '.') && (*scp != '\0'))
    scp++;
  if ((scp - &pbuf[0]) < 30)
    {
    *scp = '\0';
    scp +=2;
    }
  rrposends[cno] = mrg_appndrgnelt(&rrposns[cno],(int) strtol(&pbuf[0],NULL,10),
                                     (int) strtol(scp,NULL,10),cpgcnt);
  if (rrposns[cno] == NULL)
    rrposns[cno] = rrposends[cno];
  scnt++;
  }
return(scnt);
}

RBC_RRREAD_FORM rbc_classfyreadvsrrbsfrag(RBC_READ_ELT *rdp,
                                          MRG_REGN_ELT *rep)
/* compare rdp with rep, return a read form value
for it */
{
if (rdp->pstart > rep->rstop)
  return(RBC_rrread_is3pto);
else
  if (rdp->pstop < rep->rstart)
    return(RBC_rrread_is5pto);
  else
    if (rdp->pstart == rep->rstart)
      return(RBC_rrread_5p);
    else
      if ((rdp->readfwd && (rdp->pstop == rep->rstop)) ||
        (!rdp->readfwd && (rdp->pstop == (rep->rstop+2))))  /* allow for C^CGG overhang */
        return(RBC_rrread_3p);
      else
        if (rbc_intinrng(rep->rstart+1,rdp->pstart,rep->rstop-1) &&
             rbc_intinrng(rep->rstart+1,rdp->pstop,rep->rstop-1))
          return(RBC_rrread_internal);
        else
          if (rbc_intinrng(rdp->pstart+1,rep->rstart,rdp->pstop-1) ||
               rbc_intinrng(rdp->pstart+1,rep->rstop,rdp->pstop-1))
            return(RBC_rrread_onjoin);
          else
            return(RBC_rrread_nomap);
}

void rbc_cmpmapposrrstrct(RBC_CHRNO maxchr,
                          MRG_REGN_ELT *rrposlsts[],
                          RBC_READ_ELT *rdlst,
                          int rrrcnts[])
/* scan rdlst and compare positions with rrposlsts
(these assumed to be in increasing order for each
chromosome) and count status information into
rrrcnts[] */
{
MRG_REGN_ELT *rep;
RBC_READ_ELT *rdlp;
RBC_RRREAD_FORM rfm;

rdlp = rdlst;
while (rdlp != NULL)
  {
  if ((rdlp->chrno >= Chr_unk) && (rdlp->chrno <= maxchr))
    {
    rep = rrposlsts[rdlp->chrno];
    while (rep != NULL)
      switch (rfm = rbc_classfyreadvsrrbsfrag(rdlp,rep))
        {
        case RBC_rrread_5p:
        case RBC_rrread_3p:
        case RBC_rrread_internal:
        case RBC_rrread_onjoin:
          rrrcnts[rfm]++;
          rep = NULL;
          break;
        case RBC_rrread_is3pto:
          rep = rep->nxtregn;
          break;
        case RBC_rrread_nomap:
        case RBC_rrread_is5pto: /* read is 5' to rrbs fragment, stop and note nomap */
          rrrcnts[RBC_rrread_nomap]++;
          rep = NULL;
          break;
      }
    }
  rdlp = rdlp->nxtrdelt;
  }    
}                          

int main(int argc,
         char *argv[])
{
int ap;
char op;
SQFL_STRCT *srcsq;
int ecnts;
int tilno;
RBC_READ_ELT *rdlsts[MAXTILENO+1];
FILE *srcfl;
int modcnt;
char *chromoseq[ChrY + 1];   /* sequences of each chromosome */
int chrlens[ChrY + 1];       /* their lengths */
char *genomsqhdrstr;         /* header string for genomic sequences */
int chrcnt;
RBC_RUNPARS rpars;
FILE *posfile;
MRG_REGN_ELT *rrposns[ChrY+1];
RBC_CHRNO cno;
int rrcnt;
MRG_REGN_ELT *mrp;
int prvend;
int offst;
int rrmapcnts[RBC_rrread_onjoin + 1];
RBC_RRREAD_FORM rrrdpt;

debuglevel = RBC_dbg_none;
rpars.ostyle = RBC_out_cpg;
rpars.srcstyle = RBC_src_rmapbs;
rpars.sqcheck = -1;
rpars.chkrrmaps = 0;
rrcnt = modcnt = ecnts = 0;
srcfl = posfile = NULL;
genomsqhdrstr = NULL;
rpars.fcversion = RBC_fcv_2;
rpars.readflform = RBC_readflfm_fasta;
for (ap = 1; ap < argc; ap++)
  if (*argv[ap] == '-')   /* an option */
    switch (op = *(argv[ap]+1))
      {
      case 'b':
        if (++ap > argc)
          err_msg_die("-%c needs file name\n",op);
        else
          if ((srcfl = fopen(argv[ap],"r")) == NULL)
            err_msg_die("Can't open .BED file %s\n",argv[ap]);
        rpars.srcstyle = RBC_src_rmapbs;
        break;
      case 'B':  /* bsmap input */
        if (++ap > argc)
          err_msg_die("-%c needs file name\n",op);
        else
          if ((srcfl = fopen(argv[ap],"r")) == NULL)
            err_msg_die("Can't open bsmap file %s\n",argv[ap]);
        rpars.srcstyle = RBC_src_bsmap102;
        break;
      case 'A':
        if (++ap > argc)
          err_msg_die("-%c needs file name\n",op);
        else
          if ((srcfl = fopen(argv[ap],"r")) == NULL)
            err_msg_die("Can't open bsmap file %s\n",argv[ap]);
        rpars.srcstyle = RBC_src_bsmap12;
        break;
      case 'R':
        rpars.readflform = RBC_readflfm_fastq;
        if (++ap > argc)
          err_msg_die("-%c needs name of fastq file\n",op);
        else
          if ((srcsq = sqfl_opnsqstrct(argv[ap],SFMT_raw,"r")) == NULL)
            err_msg_die("Can't open fastq format read file '%s'\n",argv[ap]);
        break;	  
      case 'r':
        if (++ap > argc)
          err_msg_die("-%c needs name of fasta file\n",op);
        else
          if ((srcsq = sqfl_opnsqstrct(argv[ap],SFMT_fasta,"r")) == NULL)
            err_msg_die("Can't open Fasta format read file '%s'\n",argv[ap]);
        break;	  
      case 'g':
        if (++ap > argc)
          err_msg_die("-%c needs header string\n",op);
        else
          genomsqhdrstr = bas_strdup(argv[ap]);
        break;
      case 'v':   /* V3 flowcell, V2 chemistry */
        rpars.fcversion = RBC_fcv_3_2;
        break;
      case 'V':   /* V3 flowcell, V3 chemistry */
        rpars.fcversion = RBC_fcv_3_3;
        break;
      case 'C':   /* enable read/chrom sq check: advise mismatches to stderr */
        if (++ap > argc)
          err_msg_die("-%c needs int value for mismatch check\n",op);
        else
          rpars.sqcheck = (int) atol(argv[ap]);
        break;
      case 'p': /* read a positionfile */
        if (++ap > argc)
          err_msg_die("-%c needs name of mkrrgenome position file\n",op);
        else
          if ((posfile = fopen(argv[ap],"r")) == NULL)
            err_msg_die("Can't open position file '%s'\n",argv[ap]);
          else
            rpars.chkrrmaps = 1;
        break;	  
      case 'd':   /* debug on */
        debuglevel = RBC_dbg_on;
        break;
      case 'D':   /* debug^2 on */
        debuglevel = RBC_dbg_serious;
        break;
      case 'h':
        say_usage(stdout,argv[0]);
        exit(0);
        break;
      default:
        err_msg("Unknown Option: '%s\n",argv[ap]);
        say_usage(stderr,argv[0]);
        exit(1);
        break;
      }
/* should be ready to go now */
if (genomsqhdrstr != NULL)
  {
  chrcnt = rbc_getgenomesqs(genomsqhdrstr,chromoseq,&chrlens[0],(int) ChrY);
  if (debuglevel > RBC_dbg_none)
    fprintf(stdout,"%d chromosomes read\n",chrcnt);
  }
if (chrcnt <= 0)
  {
  err_msg("No chromosomes found with %s<n>\n",genomsqhdrstr);
  exit(1);
  }
if (rpars.chkrrmaps && (posfile != NULL))
  {
  for (cno = Chr_unk; cno <= ChrY; cno++)
    rrposns[cno] = NULL;
  rrcnt = rbc_readmkrrgenposns(posfile,&rrposns[Chr_unk],ChrY);
  for (cno = Chr_unk; cno <= ChrY; cno++)
    {
    mrp = rrposns[cno];
    prvend = 0;
    while (mrp != NULL)
      {
      offst = mrp->rstart - prvend - 1;
      prvend = (mrp->rstop -= offst);
      mrp->rstart -= offst;
/*      if (debuglevel > RBC_dbg_none)
        fprintf(stdout,"%s: %d..%d\n",rbc_chrno2str(cno,0),mrp->rstart,mrp->rstop); */
      mrp = mrp->nxtregn;
      }
    }
  for (rrrdpt = RBC_rrread_nomap; rrrdpt <= RBC_rrread_onjoin; rrrdpt++)
    rrmapcnts[rrrdpt] = 0;
  }
if (srcsq != NULL)
  {
  for (tilno = 0; tilno < MAXTILENO; tilno++)
    rdlsts[tilno] = NULL;
  if (debuglevel > RBC_dbg_none)
    fprintf(stdout,"Loading reads...\n");
  if (((ecnts = rbc_readnstorreads(&rpars,srcsq,rdlsts)) > 0) &&
       (srcfl != NULL))
    {
    if (debuglevel > RBC_dbg_none)
      fprintf(stdout,"%d reads stored\n",ecnts);
    modcnt = rbc_readsrcfl(&rpars,rdlsts,srcfl);
    if (debuglevel > RBC_dbg_none)
      fprintf(stdout,"%d rmapbs matches read\n",modcnt);
    }
  }
if ((ecnts > 0) && (debuglevel > RBC_dbg_none))
  for (tilno = 0; tilno < MAXTILENO; tilno++)
    fprintf(stdout,"Tile %d: %d/%d reads modified\n",tilno+1,
             rbc_cntmodrdelts(rdlsts[tilno]),
             rbc_cntrdelts(rdlsts[tilno]));
if ((ecnts > 0) && (modcnt > 0))
  {
  for (tilno = 0; tilno < MAXTILENO; tilno++)
    {
    rbc_chkbedvsgenomseq(&rpars,stdout,tilno,rdlsts[tilno],chromoseq,
                           &chrlens[0],(int) ChrY);
    if (rrcnt > 0)   /* compare read map positions with rr genome fragment positions */
      rbc_cmpmapposrrstrct(ChrY,&rrposns[Chr_unk],rdlsts[tilno],&rrmapcnts[RBC_rrread_nomap]);
    }
  if (rrcnt > 0)
    {
    rbc_rptchrout(stdout,'#',10);
    fputs(" RRBS fragment hits:\n",stdout);
    for (rrrdpt = RBC_rrread_nomap; rrrdpt <= RBC_rrread_onjoin; rrrdpt++)
      fprintf(stdout,"#\t%d\t %s\n",rrmapcnts[rrrdpt],rbc_rrreadform2str(rrrdpt));
    rbc_rptchrout(stdout,'#',50);
    fputc('\n',stdout);
    }
  }
exit(0);
}
