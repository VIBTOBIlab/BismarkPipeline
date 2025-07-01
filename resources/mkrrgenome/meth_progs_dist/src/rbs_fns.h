/* rbs_fns.h: function headers for bisulphite programs
common routines */

char *rbc_chrno2str(RBC_CHRNO cno,
                    int terse);

RBC_CHRNO rbc_str2chrno(char *cstr);
  /* return the most appropriate chromosome No for cstr */

int rbc_remaptileno(RBC_FLOWCELLVERSN fcv,
                    int fastqtileno);
/* remap flowcell V3 tiles to 1..24.  V2 just returns
the value */

int rbc_intinrng(int b1,
                 int i,
                 int b2);
/* 1 if b1 <= i <= b2 or b2 <= i <= b1 */

int rbc_modnwithn(int base,
                  int i);
/* return 1,2,...n, on a base of n */

int rbc_invrtremaptileno(RBC_FLOWCELLVERSN fcv,
                         int tno);
/* tno is a tile number in range 1..120.  reverse the
remapping back to the original number.  Check that
tno is in appropriate range and return 0 if not */

int rbc_cntflds(char *buf,
                char fldsep);
/* scan buff for fldsep, return No of distinct fields
delimited by fldsep.  Note that leading or terminal 
occurrences are not considered */

char *rbc_skiptofld(char *buf,
                    char fldsep,
                    int fldno);
/* jump to field no fldno in buf, based on
occurrences of fldsep.  return the point immediately
following the preceding delimiter, so if two successive
delimiters exist, the next will be found.
return NULL if the field doesn't exist */

int rbc_hdrstr2tile(RBC_RUNPARS *rpp,
                    char *hstr,
                    int *xpx,
                    int *ypx);
/* decompose hstr to tile, xpixel,ypixel values and
return tile no if successful, else 0.  If xpx,ypx are non-NULL,
put the pixel positions to them. */

int rbc_readlntobuf(FILE *sfl,
                    char *buf,
                    int buflen);
/* read chars from sfl up to end of line, if not '\0' or
EOF, put into buf if it is not null up to buflen,
return no chars read */
