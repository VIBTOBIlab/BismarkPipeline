/* rbs_fns.c: set of routines common to a number
of the bisulphite sequence programs: made into a 
separate file to simplify maintenance

Peter Stockwell 26-Sep-2010 */

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

char *rbc_chrno2str(RBC_CHRNO cno,
                    int terse)
{
if (terse)
  switch (cno)
    {
    case Chr1:
      return("1");
      break;
    case Chr2:
      return("2");
      break;
    case Chr3:
      return("3");
      break;
    case Chr4:
      return("4");
      break;
    case Chr5:
      return("5");
      break;
    case Chr6:
      return("6");
      break;
    case Chr7:
      return("7");
      break;
    case Chr8:
      return("8");
      break;
    case Chr9:
      return("9");
      break;
    case Chr10:
      return("10");
      break;
    case Chr11:
      return("11");
      break;
    case Chr12:
      return("12");
      break;
    case Chr13:
      return("13");
      break;
    case Chr14:
      return("14");
      break;
    case Chr15:
      return("15");
      break;
    case Chr16:
      return("16");
      break;
    case Chr17:
      return("17");
      break;
    case Chr18:
      return("18");
      break;
    case Chr19:
      return("19");
      break;
    case Chr20:
      return("20");
      break;
    case Chr21:
      return("21");
      break;
    case Chr22:
      return("22");
      break;
    case ChrX:
      return("X");
      break;
    case ChrY:
      return("Y");
      break;
    case Chr_unk:
    default:
      return("0");
      break;
    }
else
  switch (cno)
    {
    case Chr1:
      return("Chr1");
      break;
    case Chr2:
      return("Chr2");
      break;
    case Chr3:
      return("Chr3");
      break;
    case Chr4:
      return("Chr4");
      break;
    case Chr5:
      return("Chr5");
      break;
    case Chr6:
      return("Chr6");
      break;
    case Chr7:
      return("Chr7");
      break;
    case Chr8:
      return("Chr8");
      break;
    case Chr9:
      return("Chr9");
      break;
    case Chr10:
      return("Chr10");
      break;
    case Chr11:
      return("Chr11");
      break;
    case Chr12:
      return("Chr12");
      break;
    case Chr13:
      return("Chr13");
      break;
    case Chr14:
      return("Chr14");
      break;
    case Chr15:
      return("Chr15");
      break;
    case Chr16:
      return("Chr16");
      break;
    case Chr17:
      return("Chr17");
      break;
    case Chr18:
      return("Chr18");
      break;
    case Chr19:
      return("Chr19");
      break;
    case Chr20:
      return("Chr20");
      break;
    case Chr21:
      return("Chr21");
      break;
    case Chr22:
      return("Chr22");
      break;
    case ChrX:
      return("ChrX");
      break;
    case ChrY:
      return("ChrY");
      break;
    case Chr_unk:
    default:
      return("Chr??");
      break;
    }
}

RBC_CHRNO rbc_str2chrno(char *cstr)
  /* return the most appropriate chromosome No for cstr */
{
char *ep;
int cno;

cno = (int) strtol(cstr,&ep,10);
if ((ep != cstr) && (cno != 0))
  return((RBC_CHRNO) cno);
else  /* couldn't read any of string as integer */
  switch (toupper(*cstr))
    {
    case 'X':
      return(ChrX);
      break;
    case 'Y':
      return(ChrY);
      break;
    default:
      return(Chr_unk);
      break;
    }
}

int rbc_remaptileno(RBC_FLOWCELLVERSN fcv,
                    int fastqtileno)
/* remap flowcell V3 tiles to 1..24.  V2 just returns
the value */
{
int surface;  /* 1,2 */
int swath;  /* 1,2,3 */
int tdiv;

switch (fcv)
  {
  case RBC_fcv_3_2:
    surface = (int) fastqtileno/1000;
    fastqtileno %= 1000;
    swath = (int) fastqtileno/100;
    tdiv = fastqtileno%100;
    return(16*(surface-1)+8*(swath-1)+tdiv);
    break;
  case RBC_fcv_3_3:
    surface = (int) fastqtileno/1000;
    fastqtileno %= 1000;
    swath = (int) fastqtileno/100;
    tdiv = fastqtileno%100;
    return(24*(surface-1)+8*(swath-1)+tdiv);
    break;
  case RBC_fcv_2:
  default:
    return(fastqtileno);
    break;
  }
}

int rbc_intinrng(int b1,
                 int i,
                 int b2)
/* 1 if b1 <= i <= b2 or b2 <= i <= b1 */
{
if (b1 <= b2)
  return((b1 <= i) && (i <= b2));
else
  return(rbc_intinrng(b2,i,b1));
}

int rbc_modnwithn(int base,
                  int i)
/* return 1,2,...n, on a base of n */
{
int rtv;

if (base <= 0)
  return(0);
else
  if ((rtv = (i % base)) == 0)
    return(base);
  else
    return(rtv);
}

int rbc_invrtremaptileno(RBC_FLOWCELLVERSN fcv,
                         int tno)
/* tno is a tile number in range 1..120.  reverse the
remapping back to the original number.  Check that
tno is in appropriate range and return 0 if not */
{
int surface;
int swath;
int wosurf;

switch (fcv)
  {
  case RBC_fcv_3_2:
    if (rbc_intinrng(1,tno,32))
      {
      surface = (int) (tno - 1)/(8*2);
      wosurf = tno % 16;
      swath = (int)(wosurf + 7)/8;
      return(surface*1000 + swath*100 + rbc_modnwithn(8,tno));
      }
    else
      return(0);
    break;
  case RBC_fcv_3_3:
    if (rbc_intinrng(1,tno,48))
      {
      surface = (int) (tno - 1)/(8*3);
      wosurf = tno % 24;
      swath = (int)(wosurf + 7)/8;
      return(surface*1000 + swath*100 + rbc_modnwithn(8,tno));
      }
    else
      return(0);
    break;
  case RBC_fcv_2:
  default:
    if (rbc_intinrng(1,tno,MAXTILENO))
      return(tno);
    else
      return(0);
    break;
  }
}

int rbc_cntflds(char *buf,
                char fldsep)
/* scan buff for fldsep, return No of distinct fields
delimited by fldsep.  Note that leading or terminal 
occurrences are not considered */
{
char *bp;
int fcnt;
char prv;

bp = buf;
prv = '\0';
if (*bp == fldsep)
  fcnt = 0;
else
  fcnt = 1;
while ((prv = *bp) != '\0')
  {
  if (*bp == fldsep)
    fcnt++;
  bp++;
  }
if (prv == fldsep)
  fcnt--;
return(fcnt);
}

char *rbc_skiptofld(char *buf,
                    char fldsep,
                    int fldno)
/* jump to field no fldno in buf, based on
occurrences of fldsep.  return the point immediately
following the preceding delimiter, so if two successive
delimiters exist, the next will be found.
return NULL if the field doesn't exist */
{
int fcnt;
char *bp;

bp = buf;
if (*bp == fldsep)
  bp++;
fcnt = 1;
while (*bp != '\0')
  if (fcnt == fldno)
    return(bp);
  else
    {
    if (*bp == fldsep)
      fcnt++;
    bp++;
    }
return(NULL);
}

int rbc_hdrstr2tile(RBC_RUNPARS *rpp,
                    char *hstr,
                    int *xpx,
                    int *ypx)
/* decompose hstr to tile, xpixel,ypixel values and
return tile no if successful, else 0.  If xpx,ypx are non-NULL,
put the pixel positions to them. */
{
int xp;
int yp;
int tno;
char *sp;
char *ep;
int lno;
int fldcnt;

switch (rpp->readflform)
  {
  case RBC_readflfm_fasta:
    sp = hstr;
    if (*sp == '\0')    /* empty */
      return(0);
    else
      {
      if (*sp == '>')
        sp++;
      sp++;
      lno = (int) strtol(sp,&ep,10);
      if (ep == sp)
        return(0);
      else
        {
        sp = ep + 1;;
        tno = (int) strtol(sp,&ep,10);
        if (ep == sp)
          return(0);
        else
          {
          sp = ep + 1;
          xp = (int) strtol(sp,&ep,10);
          if (sp == ep)
            return(0);
          else
            {
            sp= ep + 1;
            yp = (int) strtol(sp,&ep,10);
            if (sp == ep)
              return(0);
            else
              {
              if (xpx != NULL)
                *xpx = xp;
              if (ypx != NULL)
                *ypx = yp;
              return(tno);
              }
            }
          }
        }
      }
    break;
  case RBC_readflfm_fastq:
    fldcnt = rbc_cntflds(hstr,':');
    if ((sp = rbc_skiptofld(hstr,':',(fldcnt-2))) != NULL)
      {
      tno = (int) strtol(sp,&ep,10);
      if (ep == sp)
        return(0);
      else
        {
        if ((sp = rbc_skiptofld(hstr,':',(fldcnt-1))) != NULL)
          {
          xp = (int) strtol(sp,&ep,10);
          if (ep == sp)
            return(0);
          else
            {
            if (xpx != NULL)
              *xpx = xp;
            if ((sp = rbc_skiptofld(hstr,':',fldcnt)) != NULL)
              {
              yp = (int) strtol(sp,&ep,10);
              if (ep == sp)
                return(0);
              else
                {
                if (ypx != NULL)
                  *ypx = yp;
                return(tno);
                }
              }
            else
              return(0);
            }
          }
        else
          return(0);
        }
      }
    else
      return(0);
    break;
  }         
}

int rbc_readlntobuf(FILE *sfl,
                    char *buf,
                    int buflen)
/* read chars from sfl up to end of line, if not '\0' or
EOF, put into buf if it is not null up to buflen,
return no chars read */
{
int nc;
char *dp;
int rlen;

rlen = 0;
if ((dp = buf) != NULL)
  *dp = '\0';
while (((nc = fgetc(sfl)) != (int) '\n') && (nc != EOF))
  {
  if ((buf != NULL) && (rlen < buflen))
    {
    *dp = (char) nc;
    dp++;
    *dp = '\0';
    }
  rlen++;
  }
return(rlen);
}

