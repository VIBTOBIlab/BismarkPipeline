/* mkrrgenome: to generate reduced representation genome files for all
or selected chromosomes */

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
#include "bin_cnts.h"
#include "fsm_ops.h"

/* Local defines */
/* version: first OK version: 1-Apr-2011 */
/* #define PROG_VERSION 1.0 */
/* generate CpG position list: Nov-2011 */
#define PROG_VERSION 1.1

#define MRG_DEFMIN 40
#define MRG_DEFMAX 220

typedef enum MRG_outstyle
  {
  MRG_out_rrgenomes = 0,  /* usual mode: generate RR genome files */
  MRG_out_rrlist,         /* tabbed list of valid RR fragment positions */
  MRG_out_cpginrng,       /* tabbed list of CpGs in valid RR fragments */
  MRG_out_allcpg          /* tabbed list of all CpGs */
  }
MRG_OUTSTYLE;

/* global debug & out style variables, for simplicity of access */
RBC_DBG debuglevel;
MRG_OUTSTYLE ostyle;
int glblreadlen;          /* global value for readlength */

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
fprintf(fl,"%s (v%.2f): generate reduced representation genome files for MspI digests\n",
          pnam,PROG_VERSION);
fputs("Options:\n",fl);
fputs("     -g <genomehead> dir and file string to locate genomic seq files by adding n.fa\n",fl);
fputs("     -c <C_no> restrict activity to chromosome C_no (1..22,X,Y), def = all\n",fl);
fprintf(fl,"     -M <j,k> scan for restricted rep fragment sizes between j & k residues. (def=%d..%d)\n",
          MRG_DEFMIN,MRG_DEFMAX);
fputs("     -m <j,k>   ditto, produce tab delimited list of positions to stdout\n",fl);
fputs("     -p <j,k> generate tabbed list of CpG positions for j-k fragments\n",fl);
fputs("     -P generate tabbed list of all CpG positions\n",fl);
fputs("     -G <desthead> dir & file string for output genome files (-M), completed with n.fa\n",fl);
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

int bc_fragsizeok(int fragmin,
                  int fragmax,
                  int fragsize)
/* return if fragsize fragment is acceptable,
noting that zero limits will always return true. */
{
if ((fragmin != 0) && (fragmax != 0))
  return(rbc_intinrng(fragmin,fragsize,fragmax));
else
  return(1);
}

int rbc_getnscangenomesqs(MRG_OUTSTYLE ostyle,
                          char *hdrstr,
                          FS_FSMSTRCT *cpgfsmp,
                          char *seqsarr[],
                          int seqlens[],
                          MRG_REGN_ELT *rrfrags[],
                          int chrmax,
                          int uchrno,
                          int minfrag,
                          int maxfrag)
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
char nsc;
int cpos;
FS_RESELT *frp;
FS_DATPRELT *dpep;
int prevmat;
MRG_REGN_ELT *lstend;
int thismat;
int cpgcnt;
MRG_REGN_ELT *cpglstp;
MRG_REGN_ELT *cpglstend;
MRG_REGN_ELT *cp;

rcnt = 0;
sqfnam = (char *) getmemory((nblen = strlen(hdrstr) + 16),"Sq file name buf");
for (chno = 1; chno <= chrmax; chno++)
  if ((uchrno == 0) || (uchrno == chno))
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
      if (ostyle == MRG_out_rrgenomes)
        seqsarr[chno-1] = (char *) getmemory(seqlens[chno-1]+1,"Chr Buff");
      sqfl_rewind(chsqfl);
      cpos = 0;
      prevmat = 1;
      lstend = NULL;
      fs_initrun(cpgfsmp);
      cpgcnt = 0;
      (void) sqfl_skipsqflhdr(chsqfl);
      cpglstp = cpglstend = NULL;
      while ((nsc = sqfl_getnxtres(chsqfl)) != '\0')
        {
        if ((frp = fs_procchr(cpgfsmp,nsc,tr_nares2int)) != NULL)   /* have a match */
          {
          dpep = (FS_DATPRELT *) frp->action;
          switch (dpep->ldstr)
            {
            case 2:        /* CpG */
              cpgcnt++;
              switch (ostyle)
                {
                case MRG_out_allcpg:
                  fprintf(stdout,"%s\t%d\n",rbc_chrno2str(chno,1),cpos);
                  break;
                case MRG_out_cpginrng:
                  cpglstend = mrg_appndrgnelt(&cpglstend,cpos,cpos,1);
                  if (cpglstp == NULL)
                    cpglstp = cpglstend;
                  break;
                default:
                  break;
                }
              break;
            case 4:       /* CCGG */
              thismat = cpos - dpep->ldstr + 2;
/*          if (debuglevel > RBC_dbg_none)
            fprintf(stdout,"CCGG@%d, Frag=%d\n",thismat,thismat-prevmat+1); */
              if (bc_fragsizeok(minfrag,maxfrag,thismat - prevmat + 1))
                {
                switch (ostyle)
                  {
                  case MRG_out_cpginrng:
                    cp = cpglstp;
                    while (cp != NULL)
                      {
                      fprintf(stdout,"%s\t%d\n",rbc_chrno2str(chno,1),cp->rstart);
                      cp = cp->nxtregn;
                      }
                    break;
                  case MRG_out_rrgenomes:
                  case MRG_out_rrlist:
                    lstend = mrg_appndrgnelt(&lstend,prevmat,thismat,cpgcnt);
                    if (rrfrags[chno-1] == NULL)
                      rrfrags[chno-1] = lstend;
                    break;
                  default:
                    break;
                  }
                }
              mrg_clrallregnelts(&cpglstp);
              cpglstp = cpglstend = NULL;
              prevmat = thismat + 1;
              cpgcnt = 1;       /* is a cpg in CCGG */
              break;
            }
          }
        if (ostyle == MRG_out_rrgenomes)
          {
          *(seqsarr[chno-1] + cpos) = nsc;
          *(seqsarr[chno-1] + cpos + 1) = '\0';
          }
        cpos++;
        }
      if (bc_fragsizeok(minfrag,maxfrag,cpos - prevmat + 1))
        (void) mrg_appndrgnelt(&lstend,prevmat,cpos,cpgcnt);
      sqfl_clssqstrct(chsqfl);
      if ((debuglevel > RBC_dbg_on) && (ostyle == MRG_out_rrgenomes))
        {
        fprintf(stdout,"'%.10s...'\n",seqsarr[chno-1]);
        fflush(stdout);
        }
      rcnt++;
      mrg_clrallregnelts(&cpglstp);
      }
    else
      {
      err_msg("Can't open chromosome file %s\n",sqfnam);
      seqlens[chno-1] = 0;
      }
    }
memfree(sqfnam);
return(rcnt);
}

void mrg_putreducedrepseq(char *hdrstr,
                          char *seqbuf,
                          int chrlen,
                          MRG_REGN_ELT *rrfragp,
                          int chno,
                          int minfrag,
                          int maxfrag)
/* use hdrstr to create a file name for thisp
chromosome, open as a Fasta output sequence
file and write regions in rrfragp list there to. */
{
char *sqfnam;
int nblen;
SQFL_STRCT *chsqfl;
int cpos;
MRG_REGN_ELT *rrfp;
int rc;
int lcnt;

sqfnam = (char *) getmemory((nblen = strlen(hdrstr) + 16),"Sq dest file name buf");
snprintf(sqfnam,nblen-1,"%s%s.fa",hdrstr,rbc_chrno2str((RBC_CHRNO) chno,1));
if ((chsqfl = sqfl_opnsqstrct(sqfnam,SFMT_fasta,"w")) != NULL)
  {
  rrfp = rrfragp;
  fprintf(chsqfl->sfl,">%srr reduced repr %d..%d for MspI digest Chr%s %d/%dbp CpG: %d\n",
            rbc_chrno2str((RBC_CHRNO) chno,1),minfrag,maxfrag,
            rbc_chrno2str((RBC_CHRNO) chno,1),mrg_sumregnelts(rrfragp),chrlen,
            mrg_sumcpgregnelts(rrfragp));
  rc = lcnt = 0;
  while (rrfp != NULL)
    {
    cpos = rrfp->rstart;
    while (cpos <= rrfp->rstop)
      {
      sqfl_putres(chsqfl,*(seqbuf + cpos - 1),rc,&lcnt);
      cpos++;
      }
    rrfp = rrfp->nxtregn;
    }
  sqfl_termsqfl(chsqfl,"",&lcnt);
  sqfl_clssqstrct(chsqfl);
  }
memfree(sqfnam);
}

void bc_commalst2ints(char *str,
                      int *v1,
                      int *v2)
/* expect str to contain two base10 integers separated by a comma.
attempt to read these and return as v1 & v2.  zero return on
failures */
{
char *ep;
char *fp;

*v1 = *v2 = 0;
*v1 = (int) strtol(str,&ep,10);
if (ep != str)
  *v2 = strtol(ep+1,&fp,10);
}

int main(int argc,
         char *argv[])
{
int ap;
char op;
SQFL_STRCT *srcsq;
int ecnts;
int modcnt;
char *chromoseq[ChrY + 1];   /* sequences of each chromosome */
int chrlens[ChrY + 1];       /* their lengths */
char *genomsqhdrstr;         /* header string for genomic sequences */
char *destsqhdrstr;          /* header string for destination seq files */
int chrcnt;
int minfrag;
int maxfrag;
FS_FSMSTRCT *cpgfsmp;    /* ptr to CpG searching fsm */
int uchrno;
RBC_CHRNO chrp;
MRG_REGN_ELT *rrregions[ChrY + 1];
MRG_REGN_ELT *rlstp;

debuglevel = RBC_dbg_none;
ostyle = MRG_out_rrgenomes;
modcnt = ecnts = 0;
destsqhdrstr = genomsqhdrstr = NULL;
minfrag = MRG_DEFMIN;
maxfrag = MRG_DEFMAX;
cpgfsmp = NULL;
uchrno = 0;
for (chrp = Chr1; chrp <= ChrY; chrp++)
  {
  chromoseq[chrp-1] = NULL;
  chrlens[chrp-1] = 0;
  rrregions[chrp-1] = NULL;  
  }
for (ap = 1; ap < argc; ap++)
  if (*argv[ap] == '-')   /* an option */
    switch (op = *(argv[ap]+1))
      {
      case 'c':   /* a chromosome */
        if (++ap > argc)
          err_msg_die("-%c needs a chromosome identifier (1..20,X,Y)\n",op);
        else
          if ((uchrno = rbc_str2chrno(argv[ap])) == Chr_unk)
            err_msg_die("Can't determine Chromosome '%s'\n",argv[ap]);
        break;
      case 'g':
        if (++ap > argc)
          err_msg_die("-%c needs header string\n",op);
        else
          genomsqhdrstr = bas_strdup(argv[ap]);
        break;
      case 'G':
        if (++ap > argc)
          err_msg_die("-%c needs header string\n",op);
        else
          destsqhdrstr = bas_strdup(argv[ap]);
        break;
      case 'M':
      case 'm':
      case 'p':
        switch (op)
          {
          case 'p':
            ostyle = MRG_out_cpginrng;
            break;
          case 'm':
            ostyle = MRG_out_rrlist;
            break;
          case 'M':
          default:
            ostyle = MRG_out_rrgenomes;
            break;
          }
        if (++ap < argc)
          if (*argv[ap] != '-')
            bc_commalst2ints(argv[ap],&minfrag,&maxfrag);
          else
            ap--;
        cpgfsmp = fs_initnewfsm(4,1,FS_inv_ignor);
        fs_adddatprs(cpgfsmp,"CCGG","CCGG");
        fs_adddatprs(cpgfsmp,"CG","CG");
        (void) fs_bldfsm(cpgfsmp,WLU_CASEIND,0,0,tr_int2nares,fs_chkinpstr,fs_shed2lurec);
        break;
      case 'P':   /* all CpG list */
        ostyle = MRG_out_allcpg;
        cpgfsmp = fs_initnewfsm(4,1,FS_inv_ignor);
        fs_adddatprs(cpgfsmp,"CG","CG");
        (void) fs_bldfsm(cpgfsmp,WLU_CASEIND,0,0,tr_int2nares,fs_chkinpstr,fs_shed2lurec);
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
if ((genomsqhdrstr != NULL) && (cpgfsmp != NULL))
  {
  chrcnt = rbc_getnscangenomesqs(ostyle,genomsqhdrstr,cpgfsmp,&chromoseq[0],&chrlens[0],
                                   &rrregions[0],(int) ChrY,uchrno,minfrag,maxfrag);
  if (debuglevel > RBC_dbg_none)
    fprintf(stdout,"%d chromosomes read\n",chrcnt);
  for (chrp = Chr1; chrp <= ChrY; chrp++)
    if ((uchrno == 0) || (uchrno == chrp))
      {
      rlstp = rrregions[chrp-1];
      switch (ostyle)
        {
        case MRG_out_rrlist:
          while (rlstp != NULL)
            {
            fprintf(stdout,"%s\t%d..%d (%d bp) CpG: %d\n",rbc_chrno2str(chrp,1),rlstp->rstart,
                      rlstp->rstop,(rlstp->rstop - rlstp->rstart + 1),rlstp->cpgcnt);
            rlstp = rlstp->nxtregn;
            }
          break;
        case MRG_out_rrgenomes:
          if (destsqhdrstr == NULL)
            (void) err_msg_die("Need -G <destinationheader>\n");
          else
            mrg_putreducedrepseq(destsqhdrstr,chromoseq[chrp-1],chrlens[chrp-1],
                                   rlstp,chrp,minfrag,maxfrag);
          break;
        default:
          break;
        }
      }
  }
else
  {
  if (genomsqhdrstr == NULL)
    (void) err_msg("No genome file header string (-g option) used\n");
  if (cpgfsmp == NULL)
    (void) err_msg("No fragment size range specified (-M or -m options)\n");
  exit(1);
  }
exit(0);
}
