/* scan_cpg_depth: take appropriately-formatted input files
(chrno, position, strand/status where strand/status = '+' or '-')
and count the number of +/- hits on each CpG for all or selected
chromosomes */

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <math.h>

#include "bas_fns.h"
#include "wlu_fns.h"
#include "sqfl_fns.h"
#include "sqmat_fns.h"
#include "rmapbsbed2cpg.h"
#include "rbs_fns.h"
#include "fsm_ops.h"

/* local defines */
/* the length of sequence for each linked list cluster */
#define CHR_CLUSTER_DEF 1000

typedef struct SCD_cpgelt  /* element in linked list of CpGs */
  {
  int cpgpos;              /* position this CpG */
  int metcnt;              /* count methylated CpGs */
  int unmetcnt;            /* count unmethylated */
  struct SCD_cpgelt *nxtcpgelt;    /* forward link */
  struct SCD_cpgelt *prvcpgelt;    /* rev link */
  }
SCD_CPGELT;

typedef enum SCD_datum     /* different things we can return for SCD_CPGEELT list */
  {
  SCD_datm_count,
  SCD_datm_maxpos,
  SCD_datm_metcnt,
  SCD_datm_unmetcnt,
  SCD_datm_totcnt,
  SCD_datm_maxmetcnt,
  SCD_datm_maxunmetcnt,
  SCD_datm_maxtotcnt,
  SCD_datm_metunmetcpgpos,
  SCD_datm_metcpgpos,
  SCD_datm_unmetcpgpos
  }
SCD_DATUM;

/* global debug & out style variables, for simplicity of access */
RBC_DBG debuglevel;
int allowoutby1;   /* set in order to accept out-by-one positions */

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
fprintf(fl,"%s: scan for CpG read depth for all or selected chromosomes\n",pnam);
fputs("Options:\n",fl);
fputs("     -r <posfile> read <posfile> as set of chr posit strand/meth\n",fl);
fputs("     -g <genomehead> dir and file string to locate genomic seq files by adding n.fa\n",fl);
fprintf(fl,"     -C <m> use cluster size of <m> for chromosome positions, Def=%d\n",CHR_CLUSTER_DEF);
fputs("     -l list each CpG to stdout with counts\n",fl);
fputs("     -c <n> restrict to Chromosome <n> 1..22,X,Y. Def=all\n",fl);
fputs("     -p Permit out-by-one positions (e.g. Bismark complementary strand CpGs) def=don't\n",fl);
fputs("     -S generate statistics (range, mean, std deviation, etc. for counts\n",fl);
fputs("     -H generate histogram of counts\n",fl);
fputs("     -m list missed CpG lines\n",fl);
fputs("     -n list CpG hits (Nonmisses)\n",fl);
fputs("     -z omit zero count from histogram\n",fl);
}

SCD_CPGELT *scd_appndcpgelt(SCD_CPGELT **lstrt,
                            int cposn)
/* create and append a new element to *lstrt,
init counts to 0 & set position
Return address of new element */
{
SCD_CPGELT *prev, *end_ptr;

if (lstrt != NULL)
  {
  prev = end_ptr = *lstrt;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nxtcpgelt;
    }
  end_ptr = (SCD_CPGELT *) getmemory(sizeof(SCD_CPGELT),"CpG elt");
  end_ptr->nxtcpgelt = NULL;
  end_ptr->cpgpos = cposn;
  end_ptr->metcnt = end_ptr->unmetcnt = 0;
  if (*lstrt == NULL)
    {
    *lstrt = end_ptr;
    end_ptr->prvcpgelt = NULL;
    }
  else
    {
    prev->nxtcpgelt = end_ptr;
    end_ptr->prvcpgelt = prev;
    }
  return(end_ptr);
  }
else
  return(NULL);
}

void scd_delcpgelt(SCD_CPGELT *ep,
                   SCD_CPGELT **lstrt)
/* delete ep from list *lstrt */
{
SCD_CPGELT *pt;

if (ep != NULL)
  {
  if ((pt = ep->prvcpgelt) == NULL)
    *lstrt = ep->nxtcpgelt;
  else
    pt->nxtcpgelt = ep->nxtcpgelt;
  if ((pt = ep->nxtcpgelt) != NULL)
    pt->prvcpgelt = ep->prvcpgelt;
  memfree(ep);
  }
}

void scd_clrallcntbins(SCD_CPGELT **lstrt)
  /* iteratively delete all of lstrt */
{
while (*lstrt != NULL)
  scd_delcpgelt(*lstrt,lstrt);
}

int scd_iterlstscan(SCD_CPGELT *clst,
                    SCD_DATUM datm)
/* iteratively traverse list clst, returning the appropriate
count/total */
{
int retval;
SCD_CPGELT *lp;

retval = 0;
lp = clst;
while (lp != NULL)
  {
  switch (datm)
    {
    case SCD_datm_count:
      retval++;
      break;
    case SCD_datm_maxpos:
      retval = imax(retval,lp->cpgpos);
      break;
    case SCD_datm_metcnt:
      retval += lp->metcnt;
      break;
    case SCD_datm_unmetcnt:
      retval += lp->unmetcnt;
      break;
    case SCD_datm_totcnt:
      retval += lp->metcnt + lp->unmetcnt;
      break;
    case SCD_datm_maxmetcnt:
      retval = imax(retval,lp->metcnt);
      break;
    case SCD_datm_maxunmetcnt:
      retval = imax(retval,lp->unmetcnt);
      break;
    case SCD_datm_maxtotcnt:
      retval = imax(retval,(lp->metcnt+lp->unmetcnt));
      break;
    case SCD_datm_metunmetcpgpos:
      retval++;
      break;
    case SCD_datm_metcpgpos:
      if (lp->metcnt > 0)
        retval++;
      break;
    case SCD_datm_unmetcpgpos:
      if (lp->unmetcnt > 0)
        retval++;
      break;
    default:
      break;
    }
  lp = lp->nxtcpgelt;
  }
return(retval);
}

int scd_datm4allchr(SCD_CPGELT *chrcpglsts[],
                    RBC_CHRNO mxchr,
                    SCD_DATUM datm)
/* sum datum test to all chromosomes to mxchr */
{
RBC_CHRNO chno;
int retval;

retval = 0;
for (chno = Chr1; chno <= mxchr; chno++)
  retval += scd_iterlstscan(chrcpglsts[chno-1],datm);
return(retval);
}

int scd_cntcpgelts(SCD_CPGELT *clst)
  /* recursively read list elements */
{
return(scd_iterlstscan(clst,SCD_datm_count));
}

SCD_CPGELT *scd_lastcpgelt(SCD_CPGELT *clst)
  /* iterate thru clst, returning
last element, NULL if none */
{
SCD_CPGELT *ep;

if ((ep = clst) == NULL)
  return(NULL);
else
  {
  while (ep->nxtcpgelt != NULL)
    ep = ep->nxtcpgelt;
  return(ep);
  }
}

int int_in_rng(int b1,
               int x,
               int b2)
/* return 1 if b1 <= x <= b2 or
b2 <= x <= b1 */
{
if (b1 <= b2)
  return((b1 <= x) && (x <= b2));
else
  return(int_in_rng(b2,x,b1));
}

SCD_CPGELT *scd_cpgelt4pos(SCD_CPGELT *blst,
                           int posn)
/* return a pointer to the bin that contains posn.
Only looks on blst, doesn't attempt to look more widely.
NULL if none */
{
SCD_CPGELT *bp;

bp = blst;
while (bp != NULL)
  if ((bp->cpgpos == posn) || (allowoutby1 && (abs(bp->cpgpos-posn) == 1)))
    return(bp);
  else
    bp = bp->nxtcpgelt;
/* fell off end, return NULL */
return(NULL);
}

int scd_pos2luindex(int clstrsiz,
                    int spos,
                    int widenby1)
/* return a lookup list index for spos */
{
if (clstrsiz != 0)
  {
  if (widenby1)
    spos = imax(1,spos-1);
  if ((spos % clstrsiz) == 0)
    return(imax(1,((int)( spos/clstrsiz) - 1)));
  else
    return((int) spos/clstrsiz);
  }
else
  return(0);
}

SCD_CPGELT **scd_mklookuplst(SCD_CPGELT *cpgelst,
                            int clstrsiz)
/* examine max position in cpgeltlst and
return a lookup list based on clstrsize to
find entry positions fast.  It is assumed here that
cpgelst is in ascending order of position */
{
int maxpos;
int ecnt;
SCD_CPGELT **lulst;
SCD_CPGELT *clp;
int luarrp;

if (clstrsiz == 0)
  return(NULL);
else
  {
  maxpos = scd_iterlstscan(cpgelst,SCD_datm_maxpos);
  ecnt = 1 + scd_pos2luindex(clstrsiz,maxpos,0);
  lulst = (SCD_CPGELT **) getmemory(ecnt*sizeof(SCD_CPGELT*),"List Lookup Array");
  for (luarrp = 0; luarrp < ecnt; luarrp++)
    *(lulst + luarrp) = NULL;
  clp = cpgelst;
  while (clp != NULL)
    {
    luarrp = scd_pos2luindex(clstrsiz,clp->cpgpos,0);
    if (*(lulst + luarrp) == NULL)    /* haven't seen this one yet */
      *(lulst + luarrp) = clp;
    clp = clp->nxtcpgelt;
    }
  return(lulst);
  }
}

SCD_CPGELT *scd_pos2cpgelt(int chrlen,
                           SCD_CPGELT **clulist,
                           int clstrsiz,
                           int posn)
/* cpgeltlst is an array of cpgelt pointers, one
for each clstrsize positions in chrlen.  look
up the appropriate cpg position and return a pointer to
the element for that position if it exists. NULL otherwise */
{
SCD_CPGELT *strtp;
int clstrno;

if ((posn <= chrlen) && (clstrsiz > 0))
  {
  strtp = *(clulist + (clstrno = scd_pos2luindex(clstrsiz,posn,allowoutby1)));
  while (strtp != NULL)
    if ((strtp->cpgpos == posn) ||
          (allowoutby1 && (abs(strtp->cpgpos-posn) == 1)))
      return(strtp);
    else
      if (strtp->cpgpos > posn+1) /* missed it, stop */
        return(NULL);
      else
        strtp = strtp->nxtcpgelt;
  return(NULL);
  }
else
  return(NULL);
}

int scd_readsrcfl(SCD_CPGELT *celsts[],
                  int chrlens[],
                  SCD_CPGELT **clulists[],
                  RBC_CHRNO mxchr,
                  int clstrsiz,
                  FILE *sfl,
                  RBC_CHRNO uchrno,
                  int *miscnt,
                  FILE *misfl,
                  FILE *hitfl)
/* use fscanf to read successive lines from sfl.
If uchrno is nonzero, then only do that chromosome.

if miscnt is non-NULL then return to it the number of
missed counts */
{
char chrstr[5];
int sqpos;
char sensestr[5];
int matcnt;
int chrno;
int scnt;
SCD_CPGELT *bp;

matcnt = 0;
if (miscnt != NULL)
  *miscnt = 0;
while ((scnt = fscanf(sfl,"%s %d %s",&chrstr[0],&sqpos,&sensestr[0])) != EOF)
  if ((scnt == 3) && ((chrno = rbc_str2chrno(&chrstr[0])) > Chr_unk) &&
        (chrno <= mxchr) && ((uchrno == 0) || (chrno == uchrno)))
    {
    if ((bp = scd_pos2cpgelt(chrlens[chrno-1],clulists[chrno-1],clstrsiz,sqpos)) != NULL)
      {
      if (sensestr[0] == '+')
        bp->metcnt++;
      else
        bp->unmetcnt++;
      matcnt++;
      if (hitfl != NULL)  /* tell it about this hit */
        fprintf(hitfl,"Hit: C%s\t%d\t%c\n",rbc_chrno2str(chrno,1),sqpos,sensestr[0]);
      }
    else
      {
      if (miscnt != NULL)
        *miscnt += 1;
      if (misfl != NULL)
        fprintf(misfl,"Missed: C%s\t%d\t%c\n",rbc_chrno2str(chrno,1),sqpos,sensestr[0]);
      }
    }
return(matcnt);
}

int scd_chknreadsrcfl(SCD_CPGELT *celsts[],
                      int chrlens[],
                      SCD_CPGELT **clulists[],
                      RBC_CHRNO mxchr,
                      int clstrsiz,
                      FILE *sfl,
                      RBC_CHRNO uchrno,
                      int *miscnt,
                      FILE *misfl,
                      FILE *hitfl)
/* check sfl for openness, then call bc_readsrcfl, returning number
elements read */
{
if (sfl == NULL) /* can't do anything */
  return(0);
else
  return(scd_readsrcfl(celsts,chrlens,clulists,mxchr,clstrsiz,sfl,uchrno,miscnt,misfl,hitfl));
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

int scd_scangenomesqs(char *hdrstr,
                      int chrmax,
                      RBC_CHRNO uchrno,
                      FS_FSMSTRCT *cpgfsmp,
                      int clstrsiz,
                      SCD_CPGELT *chrcpglsts[],
                      int chrlens[])
/* use hdrstr to create a series of file names, one
for each chromosome, open each as a Fasta sequence
file and create a buffer for each.  Return the
number of chromosomes processed.  Generate cpglsts for each
chromosome scanned, adding observed CpG positions thereto */
{
char *sqfnam;
int chno;
int nblen;
SQFL_STRCT *chsqfl;
int rcnt;
char nxtres;
int cpos;
FS_RESELT *frp;
FS_DATPRELT *dpep;
int clstlen;
int clstp;
SCD_CPGELT *clstep;
int hitpos;
int clstrno;
int slen;
int hcnt;

rcnt = 0;
sqfnam = (char *) getmemory((nblen = strlen(hdrstr) + 16),"Sq file name buf");
for (chno = 1; chno <= chrmax; chno++)
  if ((uchrno == 0) || (uchrno == chno))
    {
    snprintf(sqfnam,nblen-1,"%s%s.fa",hdrstr,rbc_chrno2str((RBC_CHRNO) chno,1));
    if ((chsqfl = sqfl_opnsqstrct(sqfnam,SFMT_fasta,"r")) != NULL)
      {
      chrlens[chno-1] = slen = readsrcsq(chsqfl,NULL);
      clstlen = (int) chrlens[chno-1]/clstrsiz + 1,
/*      *chrcpglsts[chno-1] = (SCD_CPGELT *) getmemory(sizeof(SCD_CPGELT*)*clstlen,
                                              "Chr cluster list");
      *clstep = (SCD_CPGELT *) getmemory(sizeof(SCD_CPGELT*)*clstlen,
                                              "Chr cluster list ends");
      for (clstp = 0; clstp < clstlen; clstp++)
        {
        *(clstep + clstp) = NULL;
        *(chrcpglsts[chno-1]+clstp) = NULL;
        } */
      sqfl_rewind(chsqfl);
      (void) sqfl_skipsqflhdr(chsqfl);
      fs_initrun(cpgfsmp);
      cpos = 0;
      hcnt = 0;
      clstep = NULL;
      while ((nxtres = sqfl_getnxtres(chsqfl)) != '\0')
        {
        cpos++;
        if ((frp = fs_procchr(cpgfsmp,nxtres,tr_nares2int)) != NULL)
          {
          dpep = (FS_DATPRELT *) frp->action;
          hitpos = cpos - dpep->ldstr + 1;
/*          clstrno = (int) hitpos/clstrsiz;
          *(clstep+clstrno) = scd_appndcpgelt(&*(clstep+clstrno),hitpos);
          if (*(chrcpglsts[chno-1] + clstrno) == NULL)
            *(chrcpglsts[chno-1]+clstrno) = *(clstep+clstrno); */
          hcnt++;
          clstep = scd_appndcpgelt(&clstep,hitpos);
          if (chrcpglsts[chno-1] == NULL)
            chrcpglsts[chno-1] = clstep;
          }
        }
      rcnt++;
      }
    else
      err_msg("Can't open chromosome file %s\n",sqfnam);
    }
memfree(sqfnam);
return(rcnt);
}

double scd_medianxcntsarr(int *cnts,
                          int maxcnt)
/* cnts is an array of integers 0..maxcnt.  total
them and work out the median. */
{
int ttl;
int cp;
int medcnt;
int medcp1;
int tcnt;
int nxtnz;

cp = 0;
ttl = 0;
while (cp <= maxcnt)
  {
  ttl += *(cnts+cp);
  cp++;
  }
medcp1 = 0;
medcnt = (int) ttl/2.0;
if ((ttl % 2) == 0)
  medcp1 = medcnt + 1;
else
  medcnt++;
tcnt = 0;
cp = 0;
while (cp <= maxcnt)
  {
  if (tcnt >= medcnt)
    {
    cp--;
    if ((medcp1 == 0) /* odd */ ||
          (tcnt != medcnt)) /* even but don't need to look for next cnt bucket */
      return((double) cp);
    else             /* even, but need to look for next filled count bucket */
      {
      nxtnz = cp + 1;
      while ((nxtnz <= maxcnt) && (*(cnts+nxtnz) == 0))
        nxtnz++;
      return((cp + imin(nxtnz,maxcnt))/2.0);
      }
    }
  tcnt += *(cnts+cp);
  cp++;
  }
return((double) maxcnt);
}

int scd_modexcntsarr(int *cnts,
                     int maxcnt)
/* cnts is an array of integers 0..maxcnt. scan
for largest no and return position as the mode */
{
int biggest;
int cp;
int bigp;

biggest = 0;
cp = 0;
bigp = 0;
while (cp <= maxcnt)
  {
  if (*(cnts+cp) > biggest)
    {
    biggest = *(cnts+cp);
    bigp = cp;
    }
  cp++;
  }
return(bigp);
}

int main(int argc,
         char *argv[])
{
int ap;
char op;
int ecnts;
FILE *srcfl;
SCD_CPGELT *chrcpglsts[ChrY];  /* bin lists for each chromosome */
int chrlens[ChrY];       /* length each chromosome */
char *genomsqhdrstr;         /* header string for genomic sequences */
int chrcnt;
RBC_CHRNO chrno;
int uchrno;    /* user has specified a chromosome, 0=>all */
int cpgcnt;
FS_FSMSTRCT *cpgfsmp;    /* ptr to CpG searching fsm */
int clstrsiz;            /* size of clusters */
int listoutput;
int statsoutput;
SCD_CPGELT *cp;
int clstrno;
int clstp;
SCD_CPGELT **clulists[ChrY];
int maxcount;
int sumx;
int sumx2;
int *distn;
int dptr;
int curcnt;
int mincnt;
int histgrmoutput;
int biggstcnt;
int nast;
int ndigs;
int zcnt;
int hdigs;
char *srcflnam;
int missedincnts;
FILE *missfl;
int matcnts;
int omitzerocnt;
FILE *hitfl;

debuglevel = RBC_dbg_none;
allowoutby1 = 0;      /* don't by default */
ecnts = 0;
srcfl = NULL;
genomsqhdrstr = NULL;
uchrno = 0;
listoutput = statsoutput = histgrmoutput = 0;
clstrsiz = CHR_CLUSTER_DEF;
cpgfsmp = fs_initnewfsm(4,1,FS_inv_ignor);
srcflnam = "";
omitzerocnt = 0;
hitfl = missfl = NULL;
for (ap = 1; ap < argc; ap++)
  if (*argv[ap] == '-')   /* an option */
    switch (op = *(argv[ap]+1))
      {
      case 'r':
        if (++ap > argc)
          err_msg_die("-%c needs file name\n",op);
        else
          if ((srcfl = fopen(argv[ap],"r")) == NULL)
            err_msg_die("Can't open read position file %s\n",argv[ap]);
          else
            srcflnam = argv[ap];
        break;
      case 'g':
        if (++ap > argc)
          err_msg_die("-%c needs header string\n",op);
        else
          genomsqhdrstr = bas_strdup(argv[ap]);
        break;
      case 'd':   /* debug on */
        debuglevel = RBC_dbg_on;
        break;
      case 'D':   /* debug^2 on */
        debuglevel = RBC_dbg_serious;
        break;
      case 'c':   /* a chromosome */
        if (++ap > argc)
          err_msg_die("-%c needs a chromosome identifier (1..20,X,Y)\n",op);
        else
          if ((uchrno = rbc_str2chrno(argv[ap])) == Chr_unk)
            err_msg_die("Can't determine Chromosome '%s'\n",argv[ap]);
        break;
      case 'l':         /* give listing */
        listoutput = 1;
        break;
      case 'S':         /* statistics output */
        statsoutput = 1;
        break;
      case 'H':         /* histogram */
        histgrmoutput = 1;
        break;
      case 'm':         /* write misses to stdout */
        missfl = stdout;
        break;
      case 'n':        /* write nonmisses to stdout */
        hitfl = stdout;
        break;
      case 'p':        /* permit out-by-one position errors */
        allowoutby1 = 1;
        break;
      case 'z':        /* omit zero count bin from histogram */
        omitzerocnt = 1;
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
if ((genomsqhdrstr != NULL) && (listoutput || statsoutput || histgrmoutput)
      && (clstrsiz > 0))
  {
  fs_adddatprs(cpgfsmp,"CG","CG");
  (void) fs_bldfsm(cpgfsmp,WLU_CASEIND,0,0,tr_int2nares,fs_chkinpstr,fs_shed2lurec);
  for (chrno = Chr1; chrno <= ChrY; chrno++)
    {
    chrcpglsts[chrno-1] = NULL;
    chrlens[chrno-1] = 0;
    clulists[chrno-1] = NULL;
    }
  (void) scd_scangenomesqs(genomsqhdrstr,ChrY,uchrno,cpgfsmp,clstrsiz,&chrcpglsts[0],&chrlens[0]);
  if (debuglevel > RBC_dbg_none)
    for (chrno = Chr1; chrno <= ChrY; chrno++)
      if (chrcpglsts[chrno-1] != NULL)
        fprintf(stdout,"C%s %d elts\n",rbc_chrno2str(chrno,1),
                  scd_iterlstscan(chrcpglsts[chrno-1],SCD_datm_count));
  for (chrno = Chr1; chrno <= ChrY; chrno++)
    if (chrcpglsts[chrno-1] != NULL)
      clulists[chrno-1] = scd_mklookuplst(chrcpglsts[chrno-1],clstrsiz);
  if (!(matcnts = scd_chknreadsrcfl(&chrcpglsts[0],&chrlens[0],&clulists[0],ChrY,clstrsiz,srcfl,uchrno,
                           &missedincnts,missfl,hitfl)))
    err_msg_die("Run failed at sourcefile read\n");
  if (listoutput)
    {
    fprintf(stdout,"Source:%s\nChromosomes:%s<n>.fa\nhit counts =%d, missed counts=%d, OutByOne=%s\n",
              srcflnam,genomsqhdrstr,matcnts,missedincnts,(allowoutby1?"yes":"no"));
    for (chrno = Chr1; chrno <= ChrY; chrno++)
      {
      cp = chrcpglsts[chrno-1];
      while (cp != NULL)
        {
        fprintf(stdout,"C%s\t%d\t+=%d\t-=%d\tsum=%d\n",rbc_chrno2str(chrno,1),
                  cp->cpgpos,cp->metcnt,cp->unmetcnt,(cp->metcnt+cp->unmetcnt));
        cp = cp->nxtcpgelt;
        }
      }
    }
  if (statsoutput)
    {
    maxcount = 0;
    fprintf(stdout,"Source:%s\nChromosomes:%s<n>.fa\nhit counts =%d, missed counts=%d, OutByOne=%s\n",
              srcflnam,genomsqhdrstr,matcnts,missedincnts,(allowoutby1?"yes":"no"));
    fprintf(stdout,"TotMethCpG=%d at %dMCpGs, TotUnMethCpG=%d at %dUnMCpGs, (CpGs Meth+UnMeth=%d)\n",
              scd_datm4allchr(chrcpglsts,ChrY,SCD_datm_metcnt),
              scd_datm4allchr(chrcpglsts,ChrY,SCD_datm_metcpgpos),
              scd_datm4allchr(chrcpglsts,ChrY,SCD_datm_unmetcnt),
              scd_datm4allchr(chrcpglsts,ChrY,SCD_datm_unmetcpgpos),
              scd_datm4allchr(chrcpglsts,ChrY,SCD_datm_metunmetcpgpos));
    for (chrno = Chr1; chrno <= ChrY; chrno++)
      maxcount = imax(maxcount,scd_iterlstscan(chrcpglsts[chrno-1],SCD_datm_maxtotcnt));
    distn = (int *) getmemory((maxcount + 1) * sizeof(int),"Distribution");
    sumx = sumx2 = cpgcnt = 0;
    mincnt = maxcount;
    for (dptr = 0; dptr <= maxcount; dptr++)
      *(distn + dptr) = 0;
    for (chrno = Chr1; chrno <= ChrY; chrno++)
      {
      cp = chrcpglsts[chrno-1];
      while (cp != NULL)
        {
        cpgcnt++;
        curcnt = cp->metcnt + cp->unmetcnt;
        sumx += curcnt;
        sumx2 += curcnt * curcnt;
        (*(distn+curcnt))++;
        mincnt = imin(mincnt,curcnt);
        cp = cp->nxtcpgelt;
        }
      }
    fprintf(stdout,"Counts: min=%d max=%d mean=%.2f sdev=%.2f\n",
              mincnt,maxcount,(double)sumx/cpgcnt,sqrt((sumx2 - sumx*sumx/cpgcnt)/(cpgcnt-1)));
    fprintf(stdout,"Median = %.1f mode=%d\n",scd_medianxcntsarr(distn,maxcount),
              scd_modexcntsarr(distn,maxcount));
    }
  if (histgrmoutput)
    {
    fprintf(stdout,"Source:%s\nChromosomes:%s<n>.fa\nhit counts =%d, missed counts=%d, OutByOne=%s\n",
              srcflnam,genomsqhdrstr,matcnts,missedincnts,(allowoutby1?"yes":"no"));
    maxcount = 0;
    for (chrno = Chr1; chrno <= ChrY; chrno++)
      maxcount = imax(maxcount,scd_iterlstscan(chrcpglsts[chrno-1],SCD_datm_maxtotcnt));
    distn = (int *) getmemory((maxcount + 1) * sizeof(int),"Distribution");
    sumx = sumx2 = cpgcnt = 0;
    mincnt = maxcount;
    for (dptr = 0; dptr <= maxcount; dptr++)
      *(distn + dptr) = 0;
    for (chrno = Chr1; chrno <= ChrY; chrno++)
      {
      cp = chrcpglsts[chrno-1];
      while (cp != NULL)
        {
        cpgcnt++;
        curcnt = cp->metcnt + cp->unmetcnt;
        sumx += curcnt;
        sumx2 += curcnt * curcnt;
        (*(distn+curcnt))++;
        mincnt = imin(mincnt,curcnt);
        cp = cp->nxtcpgelt;
        }
      }
    biggstcnt = 0;
    if (omitzerocnt)
      dptr = 1;
    else
      dptr = 0;
    while (dptr <= maxcount)
      {
      biggstcnt = imax(*(distn+dptr),biggstcnt);
      dptr++;
      }
    dptr = 0;
    ndigs = bas_digitsin(biggstcnt);
    hdigs = bas_digitsin(maxcount);
    zcnt = 0;
    while (dptr <= maxcount)
      {
      if (*(distn+dptr) > 0)
        zcnt = 0;
      else
        zcnt++;
      if (zcnt <= 3)
        {
        fprintf(stdout,"%*d %*d |",hdigs,dptr,ndigs,*(distn+dptr));
        if (omitzerocnt && (dptr == 0))
          nast = 50 - bas_digitsin(*(distn)) - 1;
        else
          nast = (int)((*(distn+dptr)*50)/ biggstcnt);
        while (nast > 0)
          {
          fputc('*',stdout);
          nast--;
          }
        if (omitzerocnt & (dptr == 0))
          fprintf(stdout,">%d",*(distn));
        fputc('\n',stdout);
        }
      else
        if (zcnt == 4)
          fputs("[...]\n",stdout);
      dptr++;
      }
    }
  }
exit(0);
}
