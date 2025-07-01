# validate_fastq.py

import sys

# Check the length of the reads
logfile = sys.argv[1]
line_count = 0
with open(logfile, "r") as file:
    for line in file:
        if line != "\n":
            line_count += 1

if line_count > 1:
    sys.exit('''
        ############### FASTQ VALIDATION ERROR ###############
        Reads with variable lengths detected in the input fastq file (run fastqc and see *lengths.txt for more details).
        Demultiplexing was possibly done with an adaptor specified in the SampleSheet.csv.
        Trimming with the <hard> option, or demultiplexing without an adaptor sequence is recommended.
        #######################################################
    ''')