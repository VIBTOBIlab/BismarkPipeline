# illum2fasta.awk: read illumina text files and convert header & seq lines to fasta format.
#  shorten headers.
#
# to run: awk -f illum2fasta.awk s1_trim75.txt > s1_trim75.fa
#
BEGIN { if (hdrstr=="") hdrstr = "s1"; ecnt=0; prvishdr=0; }

/^@/ { split($0,colonsep,":");
  split(colonsep[5],chatchsep,"#");
  printf(">%s_%s_%s_%s\n",hdrstr,colonsep[3],colonsep[4],chatchsep[1]); \
  prvishdr = 1; \
  }
/^[^@]/ { if (prvishdr)
  {
  printf("%s\n",$0);
  prvishdr=0;
  ecnt++;
  }
}
