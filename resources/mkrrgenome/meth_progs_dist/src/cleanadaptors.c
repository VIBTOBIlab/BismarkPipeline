/* cleanadaptors.c: attempt at writing something to scan for adaptor sequences
in Illumina reads.  PAS: May-2011 */

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>

#include "bas_fns.h"
#include "wlu_fns.h"
#include "sqfl_fns.h"
#include "sqmat_fns.h"
#include "rmapbsbed2cpg.h"
#include "bin_cnts.h"
#include "fsm_ops.h"

/* local defines */
/* #define PROG_VERSN 1.00 */
/* first release: June 2011 */
/* #define PROG_VERSN 1.01 */
/* don't write 0 length seqs: Jul-2011 */
/* #define PROG_VERSN 1.02 */
/* further refine length exclusion: Jul-2011 */
/* #define PROG_VERSN 1.03 */
/* allow N-filling of trims: Jul-2011 */
/* #define PROG_VERSN 1.04 */
/* allow quality lines starting with '>' */
#define PROG_VERSN 1.04
/* allow trimming back before adaptor: Jan-2012 */

#define MAX_ADAPTOR_LEN 256
#define MIN_LEAD_TRAIL 6
#define DEF_READBUFLEN 256
#define DEF_3PMARGIN 0
#define DEF_MATCHCRITERION 75.0

typedef struct CA_adhit    /* information about an adaptor hit */
  {
  int pos;                 /* where in seq line */
  FS_RESELT *resp;         /* result ptr */
  struct CA_adhit *nxthit;
  struct CA_adhit *prvhit;
  }
CA_ADHIT;   

typedef struct CA_adaptor   /* linked list element of original adaptor seqs */
  {
  char *adaptstr;          /* the sequence */
  int asqno;               /* its number */
  struct CA_adaptor *nxtadaptor;
  struct CA_adaptor *prvadaptor;
  }
CA_ADAPTOR;

typedef struct ca_adinfo  /* structure to contain an adaptor name and a pointer */
  {
  char *adname;
  CA_ADAPTOR *adptr;        /* pointer to the original full adaptor seq */
  }
CA_ADINFO;

typedef enum CA_fastqlntype
  {
  CA_fqlt_hdrline = 0,
  CA_fqlt_sqline,
  CA_fqlt_qulhdr,
  CA_fqlt_qualln,
  CA_fqlt_unk
  }
CA_FASTQLNTYPE;

typedef enum CA_outstyle
  {
  CA_out_fsmlist = 0,
  CA_out_listreads,
  CA_out_listhitreads,
  CA_out_trimhits,
  CA_out_nfill           /* put 'N's in trimmed region */
  }
CA_OUTSTYLE;

typedef enum CA_srctype
  {
  CA_src_fastq = 0,
  CA_src_fasta
  }
CA_SRCTYPE;

typedef struct CA_runpars
  {
  int trimcnt;
  CA_SRCTYPE srcstyle;
  CA_OUTSTYLE ostyle;
  int minleadtrail;
  int readflbuflen;
  int margn3p;
  int maxadlen;
  CA_ADAPTOR *adaptlst;
  int mismatches;
  int skipres;
  float matcriterion;
  int mintrimlen;     /* minimum length for output of trimmed seqs, 0=>no limit */
  int trimback;
  }
CA_RUNPARS;

/* local globals */

int debuglevel;
/* int trimcnt; */

/* code ... */

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
fprintf(fl,"%s v%.2f: scan Illumina reads for adaptor seqs: Trim FASTA/FASTQ files\n",
          pnam,PROG_VERSN);
fputs("Options:\n",fl);
fputs("     -i <adaptorfile> file of adaptor seqs 1/line)\n",fl);
fputs("         (def=\"AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATC\": 40bp of Universal adaptor)\n",fl);
fprintf(fl,"     -l max adaptor length (def=%d)\n",MAX_ADAPTOR_LEN);
fprintf(fl,"     -m min leading match for adaptor (def=%d)\n",MIN_LEAD_TRAIL);
fprintf(fl,"     -M margin at read end for trimming (def=%d)\n",DEF_3PMARGIN);
fprintf(fl,"     -R readfile buffer length (def=%d)\n",DEF_READBUFLEN);
fputs("     -s <skipres> skip <skipres> on each read before checking matches\n",fl);
fprintf(fl,"     -p <%%> %% match threshold with adaptor sequence for hit (def=%.1f%%)\n",
          DEF_MATCHCRITERION);
fputs("     -S enable single base mismatches (def=disallow)\n",fl);
fputs("     -L print fsm to stdout\n",fl);
fputs("     -F <readfile>: run scan on <readfile> FASTQ/FASTA fmt, trim matching ends\n",fl);
fputs("     -z <readfile>: as -F but don't check length of trimmed reads\n",fl);
fputs("     -x <lengthlimit>: only save trimmed reads exceeding length limit (def=1)\n",fl);
fputs("     -t <3'trimlength>: take further 3'trimlength bases before adaptor match (def=0)\n",fl);
fputs("     -N <readfile>: fill lines with 'N's rather than trimming\n",fl);
fputs("     -f <readfile>: run scan on <readfile> FASTQ/FASTA fmt, show all reads, indicate matches\n",fl);
fputs("     -H <readfile>: run scan on <readfile> FASTQ/FASTA fmt, indicate matches on hit reads only\n",fl);
fputs("         if <readfile> is '-', then use stdin\n",fl);
}

char tr_int2nares(int iv)
  /* return a nucleic acid residue for iv */
{
switch (iv)
  {
  case 0:
    return('a');
    break;
  case 1:
    return('c');
    break;
  case 2:
    return('g');
    break;
  case 3:
    return('t');
    break;
  default:
    return('?');
    break;
  }
}

int tr_nares2int(char res)
  /* return an int value 0..3 for res, -1 for unknown */
{
switch (toupper(res))
  {
  case 'A':
    return(0);
    break;
  case 'C':
    return(1);
    break;
  case 'G':
    return(2);
    break;
  case 'T':
  case 'U':
    return(3);
    break;
  default:
    return(-1);
    break;
  }
}

CA_ADAPTOR *ca_appndadaptor(CA_ADAPTOR **alst,
                            char *astr,
                            int asno)
/* append a new string element to *alst.
strdup astr in order that a conserved 
value is created.  Return address of new element
in case it is useful */
{
CA_ADAPTOR *prev, *end_ptr;

if (alst != NULL)
  {
  prev = end_ptr = *alst;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nxtadaptor;
    }
  end_ptr = (CA_ADAPTOR *) getmemory(sizeof(CA_ADAPTOR),"Adaptor elt");
  end_ptr->nxtadaptor = NULL;
  end_ptr->adaptstr = bas_strdup(astr);
  end_ptr->asqno = asno;
  if (*alst == NULL)
    {
    *alst = end_ptr;
    end_ptr->prvadaptor = NULL;
    }
  else
    {
    prev->nxtadaptor = end_ptr;
    end_ptr->prvadaptor = prev;
    }
  return(end_ptr);
  }
else
  return(NULL);
}

void ca_deladaptelt(CA_ADAPTOR *ep,
                    CA_ADAPTOR **lstrt)
/* delete ep from list *lstrt */
{
CA_ADAPTOR *pt;

if (ep != NULL)
  {
  if ((pt = ep->prvadaptor) == NULL)
    *lstrt = ep->nxtadaptor;
  else
    pt->nxtadaptor = ep->nxtadaptor;
  if ((pt = ep->nxtadaptor) != NULL)
    pt->prvadaptor = ep->prvadaptor;
  if (ep->adaptstr != NULL)
    memfree(ep->adaptstr);
  memfree(ep);
  }
}

void ca_clralladaptors(CA_ADAPTOR **lstrt)
  /* iteratively delete all of lstrt */
{
while (*lstrt != NULL)
  ca_deladaptelt(*lstrt,lstrt);
}

int ca_buildadaptorfsm(char *sqbuf,
                       int sqno,
                       CA_ADAPTOR *adpt,
                       CA_RUNPARS *rpars,
                       FS_FSMSTRCT *adfsm)
/* scan seq sqbuf as an adaptor sequence, making
FSM adfsm in accordance with parameters.
return No. of fsm entries */
{
int sqlen;
char *leadfragbuf;
int sqp;
char *lablstr;
SQ_RESTYPE sres;
int frgp;
char cacheres;
char subres;
int fsmcnt;
CA_ADINFO *adinfo;

fsmcnt = 0;
if ((sqlen = strlen(sqbuf)) > 0)
  {
  leadfragbuf = (char *) getmemory(sqlen+1,"lead seq buf");
  for (sqp = rpars->minleadtrail; sqp <= sqlen; sqp++)
    {
    (void) memcpy(leadfragbuf,sqbuf,sqp);
    *(leadfragbuf + sqp) = '\0'; 
    for (frgp = 0; frgp < sqp; frgp++)
      if (rpars->mismatches)
        {
        cacheres = *(leadfragbuf+frgp);
        for (sres = RES_a; sres <= RES_t; sres++)
          {
          *(leadfragbuf+frgp) = sqfl_restype2chr(sres);
          asprintf(&lablstr,"Ad_%d_l%d_%c%d%c",sqno,sqp,cacheres,(frgp+1),
                     sqfl_restype2chr(sres));
          adinfo = (CA_ADINFO *) getmemory(sizeof(CA_ADINFO),"Ad info");
          adinfo->adname = bas_strdup(lablstr);
          adinfo->adptr = adpt;
          fs_adddatprs(adfsm,leadfragbuf,adinfo);
          fsmcnt++;
          memfree(lablstr);
          }
        *(leadfragbuf+frgp) = cacheres;
        }
      else
        {
        asprintf(&lablstr,"Ad_%d_l%d",sqno,sqp);
        adinfo = (CA_ADINFO *) getmemory(sizeof(CA_ADINFO),"Ad info");
        adinfo->adname = bas_strdup(lablstr);
        adinfo->adptr = adpt;
        fs_adddatprs(adfsm,leadfragbuf,adinfo);
        fsmcnt++;
        memfree(lablstr);
        }
     }
  }
memfree(leadfragbuf);
return(fsmcnt);
}

int ca_readadaptorfile(FILE *sfl,
                       CA_RUNPARS *rpars,
                       FS_FSMSTRCT *adfsm)
/* scan sfl for adaptor sequences one seq/line.  Load
each original seq and all partitions starting from
position minleadtail up to minleadtail from end into
adfsm. return no of seqs read */
{
int sqno;
int sqlen;
char *sqbuf;
CA_ADAPTOR *alstptr;
char *sbp;

sqno = 0;
sqbuf = (char *) getmemory(rpars->maxadlen+1,"input seq buf");
while ((sbp = fgets(sqbuf,rpars->maxadlen,sfl)) != NULL)
  {
  while (*sbp != '\0')
    {
    if (*sbp == '\n')
      *sbp = '\0';
    sbp++;
    }
  if ((sqlen = strlen(sqbuf)) > 0)
    {
    sqno++;
    alstptr = ca_appndadaptor(&rpars->adaptlst,sqbuf,sqno);
    (void) ca_buildadaptorfsm(sqbuf,sqno,alstptr,rpars,adfsm);
    }
  }
memfree(sqbuf);
return(sqno);
}

void ca_nl(FILE *fl)
{
fputc('\n',fl);
}

void ca_prtctxt(FILE *fl,
                void *p)
{
(void) fs_lstcstr(fl,p);
}

int ca_inrankordr(CA_ADHIT *e1p,
                  CA_ADHIT *e2p)
/* if e1 preceeds e2 in rank order, then return 1. */
{
if (e1p == NULL)
  return(0);
else
  if (e2p == NULL)
    return(1);
  else
    return((e1p->resp == e2p->resp) && (e1p->pos < e2p->pos));
}

CA_ADHIT *ca_insertnewhit(CA_ADHIT **alst,
                          int lpos,
                          FS_RESELT *rep,
                          int (* preceeds_fn)(CA_ADHIT *e1p,
                                              CA_ADHIT *e2p))
/* scan through *alst, if exists, finding the first element
which is after lpos, insert a new element before it, returning
the address of the new element.  If we fall off the end of the list,
then append the new element there. */
{
CA_ADHIT *prvelt;
CA_ADHIT *ep;
CA_ADHIT *newp;

prvelt = NULL;
ep = *alst;
newp = (CA_ADHIT *) getmemory(sizeof(CA_ADHIT),"ADhitelt");
newp->pos = lpos;
newp->resp = rep;
while (ep != NULL)
  if ((*preceeds_fn)(newp,ep))
    {
    newp->nxthit = ep;
    if (ep->prvhit == NULL)
      {
      *alst = newp;
      newp->prvhit = NULL;
      }
    else
      {
      ep->prvhit->nxthit = newp;
      newp->prvhit = ep->prvhit;
      }
    ep->prvhit = newp;
    return(newp);
    }
  else
    {
    prvelt = ep;
    ep = ep->nxthit;
    }
/* haven't found it, stick on end */
newp->nxthit = NULL;
newp->prvhit = prvelt;
if (prvelt == NULL)
  *alst = newp;
else
  prvelt->nxthit = newp;
return(newp);
}

CA_ADHIT *ca_hit4pos5p(CA_ADHIT *hlst,
                       int hpos)
/* return the first hit element on hlst which
matches 5' hpos */
{
CA_ADHIT *hp;

hp = hlst;
while (hp != NULL)
  if (hp->pos == hpos)
    return(hp);
  else
    hp = hp->nxthit;
return(NULL);
}

int ca_len4fs_reselt(FS_RESELT *rep)
  /* return the data length for rep */
{
FS_DATPRELT *dpep;

if ((rep != NULL) && ((dpep = (FS_DATPRELT *) rep->action) != NULL))
  return(dpep->ldstr);
else
  return(0);
}

void ca_killhitelt(CA_ADHIT **alst,
                   CA_ADHIT *hep)
/* remove *hep from alst */
{
if (hep != NULL)
  {
  if (hep->prvhit == NULL)
    *alst = hep->nxthit;
  else
    hep->prvhit->nxthit = hep->nxthit;
  if (hep->nxthit != NULL)
    hep->nxthit->prvhit = hep->prvhit;
  memfree(hep);
  }
}

void ca_killhitlst(CA_ADHIT **alst)
  /* iteratively kill all of alst */
{
while (*alst != NULL)
  ca_killhitelt(alst,*alst);
}

CA_ADHIT *ca_hit4pos3p(CA_ADHIT *hlst,
                       int pos3p)
/* return the first hit element on hlst which
matches 3' pos3p */
{
CA_ADHIT *hp;
int p3p;

hp = hlst;
while (hp != NULL)
  {
  p3p = hp->pos + ca_len4fs_reselt(hp->resp) - 1;
  if (p3p == pos3p)
    return(hp);
  else
    hp = hp->nxthit;
  }
return(NULL);
}

CA_ADHIT *ca_hitat3pend(CA_ADHIT *hlst,
                        int sqlen,
                        int margin)
/* if hlst contains any hits which lie within
margin of end of seq line (sqlen long), then
return point to that hit element, else NULL */
{
CA_ADHIT *hp;

hp = hlst;
while (hp != NULL)
  if ((hp->pos + ca_len4fs_reselt(hp->resp)) >= (sqlen-margin))
    return(hp);
  else
    hp = hp->nxthit;
return(NULL);
}

float ca_hitpercntmatch(char *sqread,
                        CA_ADHIT *hitp,
                        int *olapcnt,
                        int *matcnt)
/* scan this hit vs sqread, counting
number of matched residues.  Return
% of matches by end of read or original
adaptor seq.  If olapcnt & matcnt are non-NULL
then return appropriate counts to them.
return 0.0 for issues */
{
int ocnt;
int mcnt;
char *ap;
char *sp;
FS_DATPRELT *dpep;
CA_ADINFO *adinfop;

if (hitp != NULL)
  {
  dpep = (FS_DATPRELT *) hitp->resp->action;
  adinfop = (CA_ADINFO *) dpep->stxt;
  ap = adinfop->adptr->adaptstr;
  sp = sqread + hitp->pos - 1;
  ocnt = mcnt = 0;
  while ((*ap != '\0') && (*sp != '\0'))
    {
    if (toupper(*ap++) == toupper(*sp++))
      mcnt++;
    ocnt++;
    }
  if (olapcnt != NULL)
    *olapcnt = ocnt;
  if (matcnt != NULL)
    *matcnt = mcnt;
  if (ocnt > 0)
    return((float) mcnt*100.0/ocnt);
  else
    return(0.0);
  }
else
  return(0.0);
}

int ca_hitmeetscriteria(char *sqread,
                         CA_ADHIT *hitp,
                         CA_RUNPARS *rpars)
/* return 1 if this match meets criteria for a match */
{
float mprcnt;
int p3gap;
FS_DATPRELT *dpep;
CA_ADINFO *adinfop;

mprcnt = ca_hitpercntmatch(sqread,hitp,NULL,NULL);
if ((hitp!= NULL) && (mprcnt >= rpars->matcriterion))
  if (rpars->margn3p <= 0)
    return(1);
  else
    {
    dpep = (FS_DATPRELT *) hitp->resp->action;
    adinfop = (CA_ADINFO *) dpep->stxt;
    p3gap = hitp->pos + strlen(adinfop->adptr->adaptstr) - strlen(sqread);
    return(p3gap <= rpars->margn3p);
    }
return(0);
}
int ca_ahitmeetscriteria(char *sqread,
                         CA_ADHIT *hlst,
                         CA_RUNPARS *rpars)
/* check all elements of hitlst and return if any
matches criteria for a match */
{
if (hlst == NULL)
  return(0);
else
  if (ca_hitmeetscriteria(sqread,hlst,rpars))
    return(1);
  else
    return(ca_ahitmeetscriteria(sqread,hlst->nxthit,rpars));
}

CA_ADHIT *ca_5pmostvalidhit(CA_ADHIT *hlst,
                            char *sqread,
                            CA_RUNPARS *rpars)
/* scan any hits in hlst, checking
for validity.  Return the 5'-most such
hit, NULL if none */
{
CA_ADHIT *h5p;
CA_ADHIT *hp;

h5p = NULL;
hp = hlst;
while (hp != NULL)
  {
  if (ca_hitmeetscriteria(sqread,hp,rpars))
    if ((h5p == NULL) || (h5p->pos > hp->pos))
      h5p = hp;
  hp = hp->nxthit;
  }
return(h5p);
}

void ca_rptchrout(FILE *dfl,
                  char outc,
                  int ccnt)
/* write ccnt outcX to dfl */
{
int cc;

cc = ccnt;
while (cc > 0)
  {
  fputc(outc,dfl);
  cc--;
  }
}

void ca_writedataout(FILE *dfile,
                     char *linbuf[],
                     CA_FASTQLNTYPE mxlt,
                     CA_ADHIT *hitlst,
                     CA_RUNPARS *rpars)
/* logic to write out the data from
linbuf[] according to requirements */
{
CA_ADHIT *hlp;
int bcnt;
FS_DATPRELT *dpep;
CA_FASTQLNTYPE fqlt;
int slen;
char *sp;
CA_ADINFO *adinfop;
char *ap;
int olapcnt;
int matcnt;
float pcntmat;
int origlen;

if ((rpars->ostyle != CA_out_listhitreads) || 
      ca_ahitmeetscriteria(linbuf[CA_fqlt_sqline],hitlst,rpars))
  {
  fqlt = CA_fqlt_hdrline;
  origlen = slen = strlen(linbuf[CA_fqlt_sqline]);
  if (((rpars->ostyle == CA_out_trimhits) || (rpars->ostyle == CA_out_nfill)) && 
        ((hlp = ca_5pmostvalidhit(hitlst,linbuf[CA_fqlt_sqline],rpars)) != NULL))
    {
    if ((slen = hlp->pos - 1) > 0)
      {
      rpars->trimcnt++;
      slen -= rpars->trimback;
      }
    }
  if ((slen >= 0) &&
        (((rpars->ostyle == CA_out_trimhits) || (rpars->ostyle == CA_out_nfill)) &&
        ((slen >= rpars->mintrimlen) || (rpars->mintrimlen <= 0)))
         || ((rpars->ostyle != CA_out_trimhits) && (rpars->ostyle != CA_out_nfill)))
    {
    while (fqlt <= mxlt)
      {
      switch (fqlt)
        {
        case CA_fqlt_hdrline:
        case CA_fqlt_qulhdr:
          fputs(linbuf[fqlt],dfile);
          fputc('\n',dfile);
          break;
        case CA_fqlt_sqline:
          bcnt = slen;
          sp = linbuf[fqlt];
          while (bcnt > 0)
            {
            fputc(*sp++,dfile);
            bcnt--;
            }
          if (rpars->ostyle == CA_out_nfill)
            ca_rptchrout(dfile,'N',(origlen - slen));
          fputc('\n',dfile);
/*        fprintf(dfile,"%*s\n",slen,linbuf[fqlt]); */
          if ((rpars->ostyle != CA_out_trimhits) && (rpars->ostyle != CA_out_nfill))
            {
            hlp = hitlst;
            while (hlp != NULL)
              {
              if (ca_hitmeetscriteria(linbuf[CA_fqlt_sqline],hlp,rpars))
                {
                bcnt = hlp->pos - 1;
                while (bcnt-- > 0)
                  fputc(' ',dfile);
                bcnt = ca_len4fs_reselt(hlp->resp);
                while (bcnt-- > 0)
                  fputc('^',dfile);
                fputc('\n',dfile);
                dpep = (FS_DATPRELT *) hlp->resp->action;
                adinfop = (CA_ADINFO *) dpep->stxt;
                bcnt = hlp->pos - 1;
                while (bcnt-- > 0)
                  fputc(' ',dfile);
                bcnt = ca_len4fs_reselt(hlp->resp);
                ap = adinfop->adptr->adaptstr;
                sp = linbuf[fqlt] + hlp->pos - 1;
                while ((*ap != '\0') && (*sp != '\0'))
                  {
                  if (toupper(*ap) == toupper(*sp++))
                    fputc(toupper(*ap),dfile);
                  else
                    fputc(tolower(*ap),dfile);
                  ap++;
                  }
                pcntmat = ca_hitpercntmatch(linbuf[fqlt],hlp,&olapcnt,&matcnt);
                fprintf(dfile," %d/%d match (%.1f%%) ",matcnt,olapcnt,pcntmat);
                ca_prtctxt(dfile,adinfop->adname);
                fputc('\n',dfile);
                }
              hlp = hlp->nxthit;
              }
            }
          break;
        case CA_fqlt_qualln:
          bcnt = slen;
          sp = linbuf[fqlt];
          while (bcnt > 0)
            {
            fputc(*sp++,dfile);
            bcnt--;
            }
          if (rpars->ostyle == CA_out_nfill)
            ca_rptchrout(dfile,';',(origlen - slen));
          fputc('\n',dfile);
/*        fprintf(dfile,"%*s\n",slen,linbuf[fqlt]); */
          break;
        }
      fqlt++;
      }
    }
  }
}

void ca_scanreads(FILE *rdfl,
                  FS_FSMSTRCT *adfsm,
                  CA_RUNPARS *rpars,
                  FILE *dfile)
/* scan thru rdfl as either fasta or fastq fmt, scanning
each read from skipres.  Note any matches to adfsm and
either tell dfile about it or trim reads to dfile */
{
char nc;
int lpos;
FS_RESELT *frp;
int hpos;
CA_ADHIT *hitlst;
CA_FASTQLNTYPE fqlt;
char *linbufs[CA_fqlt_unk];
char *dstp;
CA_FASTQLNTYPE mxlt;
CA_ADHIT *hlp;
int lcnt;

for (fqlt = CA_fqlt_hdrline; fqlt < CA_fqlt_unk; fqlt++)
  {
  linbufs[fqlt] = (char *) getmemory(rpars->readflbuflen+1,"ReadBuf");
  *linbufs[fqlt] - '\0';
  }
mxlt = CA_fqlt_hdrline;
fqlt = CA_fqlt_unk;
lpos = 1;
lcnt = 0;
hitlst = NULL;
while ((nc = fgetc(rdfl)) != EOF)
  {
  if (lpos == 1)
    {
    if (lcnt == 0)    /* classify file type at this point */
      switch (nc)
        {
        case '>':
          rpars->srcstyle = CA_src_fasta;
          break;
        case '@':
        default:
          rpars->srcstyle = CA_src_fastq;
          break;
        }
    switch (nc)
      {
      case '>':       /* fasta hdr line */
        if ((rpars->srcstyle == CA_src_fasta) && (lcnt%2 == 0))
          { /* on fasta header, put out prev data */
          if (lcnt > 0)
            {
            ca_writedataout(dfile,&linbufs[CA_fqlt_hdrline],mxlt,hitlst,rpars);
            ca_killhitlst(&hitlst);
            }
          fqlt = CA_fqlt_hdrline;
          fs_initrun(adfsm);
          }
        else  /* not header, so just process */
          {
          fqlt++;
          if (fqlt > mxlt)
            mxlt = fqlt;
          }
        break;
      case '@':
        if ((rpars->srcstyle == CA_src_fastq) && (lcnt%4==0))
          {        /* fastq hdr line, need to put out prev data */
          if (lcnt > 0)
            {
            ca_writedataout(dfile,&linbufs[CA_fqlt_hdrline],mxlt,hitlst,rpars);
            ca_killhitlst(&hitlst);
            }
          fqlt = CA_fqlt_hdrline;
          fs_initrun(adfsm);
          }
        else  /* not header, so just process */
          {
          fqlt++;
          if (fqlt > mxlt)
            mxlt = fqlt;
          }
        break;
      case '+':       /* fastq qual hdr line */
        fqlt = CA_fqlt_qulhdr;
        break;
      default:
        fqlt++;
        if (fqlt > mxlt)
          mxlt = fqlt;
        break;
      }
    dstp = linbufs[fqlt];
    *dstp = '\0';
    }
  lpos++;
  if (nc == '\n')   /* at end of line */
    {
    lcnt++;
    lpos = 1;
    }
  else
    switch (fqlt)
      {
      case CA_fqlt_sqline:
        if (lpos > rpars->skipres)
          if ((frp = fs_procchr(adfsm,nc,tr_nares2int)) != NULL)
            {
            hpos = lpos - ca_len4fs_reselt(frp);
            if (((hlp = ca_hit4pos5p(hitlst,hpos)) != NULL) &&
                (ca_len4fs_reselt(frp) > ca_len4fs_reselt(hlp->resp)))
              hlp->resp = frp;   /* only want longest match this pos */
            else
              (void) ca_insertnewhit(&hitlst,hpos,frp,ca_inrankordr);
            }
      case CA_fqlt_hdrline:
      case CA_fqlt_qulhdr:
      case CA_fqlt_qualln:
        *dstp = nc;
        dstp++;
        *dstp = '\0';
        break;
      default:
        break;
      }
  }
ca_writedataout(dfile,&linbufs[CA_fqlt_hdrline],mxlt,hitlst,rpars);
ca_killhitlst(&hitlst);
for (fqlt = CA_fqlt_hdrline; fqlt < CA_fqlt_unk; fqlt++)
  memfree(linbufs[fqlt]);
}

int main(int argc,
         char *argv[])
{
FILE *adsrcfl;
FS_FSMSTRCT *adfsm;
int adcnt;
int ap;
char op;
char *fendp;
FILE *readfl;
int mismatches;
CA_ADAPTOR *adaptlst;
CA_RUNPARS runpars;

debuglevel = 1;
adsrcfl = NULL;
adfsm = fs_initnewfsm(4,1,FS_inv_reset);
readfl = stdin;
runpars.readflbuflen = DEF_READBUFLEN;
runpars.minleadtrail = MIN_LEAD_TRAIL;
runpars.maxadlen = MAX_ADAPTOR_LEN;
runpars.margn3p = DEF_3PMARGIN;
runpars.mismatches = 0;
runpars.adaptlst = NULL;
runpars.trimcnt = 0;
runpars.ostyle = CA_out_listreads;
runpars.skipres = 0;
runpars.matcriterion = DEF_MATCHCRITERION;
runpars.mintrimlen = 1;
runpars.trimback = 0;
runpars.srcstyle = CA_src_fastq;
for (ap = 1; ap < argc; ap++)
  if (*argv[ap] == '-')   /* an option */
    switch (op = *(argv[ap]+1))
      {
      case 'L':     /* print/list fsm */
        runpars.ostyle = CA_out_fsmlist;
        break;
      case 'l':    /* adaptor max len */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.maxadlen = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (runpars.maxadlen <= 0)
              err_msg_die("Invalid max adaptor length '%s'\n",argv[ap]);
          }
        break;
      case 'm':    /* min lead/trail */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.minleadtrail = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (runpars.minleadtrail <= 0)
              err_msg_die("Invalid min lead/trail length '%s'\n",argv[ap]);
          }
        break;
      case 'M':    /* 3p margin for trim */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.margn3p = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (runpars.margn3p <= 0)
              err_msg_die("Invalid 3' trim margin '%s'\n",argv[ap]);
          }
        break;
      case 'R':    /* read buff length */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.readflbuflen = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (runpars.readflbuflen <= 0)
              err_msg_die("Invalid read file buffer length '%s'\n",argv[ap]);
          }
        break;
      case 'p':    /* match criterion percentage */
        if (++ap > argc)
          err_msg_die("-%c needs float\n",op);
        else
          {
          runpars.matcriterion = strtof(argv[ap],&fendp);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert float '%s' for -%c\n",argv[ap],op);
          else
            if (runpars.matcriterion < 0.0)
              err_msg_die("Invalid match criterion value '%s'\n",argv[ap]);
          }
        break;
      case 's':    /* skip residues each read */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.skipres = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (runpars.skipres <= 0)
              err_msg_die("Invalid read skip length '%s'\n",argv[ap]);
          }
        break;
      case 'S':    /* disable single base mismatches */
        runpars.mismatches = 1;
        break;
      case 'i':    /* user name for adaptor file */
        if (++ap > argc)
          err_msg_die("-%c needs adaptor seq file name\n",op);
        else
          if ((adsrcfl = fopen(argv[ap],"r")) == NULL)
            err_msg_die("Can't open adaptor seq file '%s'\n",argv[ap]);
        break;
      case 'H':    /* report hit matches only */
      case 'F':    /* Trim reads in file */
      case 'f':    /* Report all reads matches & other */
      case 'z':    /* include 0 length trimmed reads */
      case 'N':    /* N-fill lines */
        switch (op)
          {
          case 'H':
            runpars.ostyle = CA_out_listhitreads;
            break;
          case 'F':
            runpars.ostyle = CA_out_trimhits;
            break;
          case 'z':
            runpars.ostyle = CA_out_trimhits;
            runpars.mintrimlen = 0;
            break;
          case 'N':
            runpars.ostyle = CA_out_nfill;
            break;
          case 'f':
          default:
            runpars.ostyle = CA_out_listreads;
            break;
          }
        if (++ap > argc)
          err_msg_die("-%c needs read file name\n",op);
        else
          if (strcmp(argv[ap],"-") != 0)
            if ((readfl = fopen(argv[ap],"r")) == NULL)
              err_msg_die("Can't open read file '%s'\n",argv[ap]);
        break;
      case 'x':    /* set minimum trimmed seq length */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.mintrimlen = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          }
        break;
      case 't':    /* trimback */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.trimback = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          }
        break;
      case 'd':   /* debug on */
        debuglevel = 1;
        break;
      case 'D':   /* debug^2 on */
        debuglevel = 2;
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
if (adsrcfl == NULL)
  {
  (void) ca_appndadaptor(&runpars.adaptlst,"AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATC",1);
  adcnt = ca_buildadaptorfsm(runpars.adaptlst->adaptstr,runpars.adaptlst->asqno,
                               runpars.adaptlst,&runpars,adfsm);
  }
else
  adcnt = ca_readadaptorfile(adsrcfl,&runpars,adfsm);
(void) fs_bldfsm(adfsm,WLU_CASEIND,0,0,tr_int2nares,fs_chkinpstr,fs_shed2lurec);
if (runpars.ostyle == CA_out_fsmlist)
  {
  fs_lststttbl(stdout,adfsm);
  fs_lststrpr(stdout,adfsm->prlst,1,ca_prtctxt,ca_nl);
  }
else
  ca_scanreads(readfl,adfsm,&runpars,stdout);
exit(0);
}
