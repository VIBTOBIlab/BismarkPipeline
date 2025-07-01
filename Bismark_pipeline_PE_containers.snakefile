# Parallelisation options
import multiprocessing
import sys
import itertools
import os
import collections
import json
import glob
cpuCount = (multiprocessing.cpu_count() - 2)

# Import config and make report
configfile: "config_bismark.yaml"
report: "report/workflow.rst"

if config["sequencer"] == "NextSeq":
    OpticalDupsPixelDistance = 100
elif config["sequencer"] == "HiSeq":
    OpticalDupsPixelDistance = 2500
elif config["sequencer"] == "NovaSeq":
    OpticalDupsPixelDistance = 12000
else:
    sys.exit("Specify sequencer (NextSeq, HiSeq, NovaSeq)")

if config["dualindex"] == False:
    ReadRegex = '[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+)_[0-9]+:[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+'
else:
    ReadRegex = '[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+)_[0-9]+:[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+\+[a-zA-Z0-9]+'

if config["zerocov"]== True:
    zerocov = "--zero_based"


zerocov = ""

fastqfolder = config["fastqfolder"].rstrip("/")
IDS, = glob_wildcards(fastqfolder + "/{sample}.fastq.gz")
R1IDS, = glob_wildcards(fastqfolder + "/{sample}_R1_001.fastq.gz")

rule all:
    input:
        "multiqc_report.html"

if config["trimming"] == "hard":
    rule adaptortrimming:
        input:
            fq1 = fastqfolder + "/{sample}_R1_001.fastq.gz",
            fq2 = fastqfolder + "/{sample}_R2_001.fastq.gz"
        output:
            temp("trimmed_reads/{sample}_R1_001_val_1.fq.gz"),
            temp("trimmed_reads/{sample}_R2_001_val_2.fq.gz")
        params:
            logfile = "{sample}_lengths.txt"
        container:
            "docker://quay.io/biocontainers/trim-galore:0.6.6--0"
        shell:
            """
            trim_galore --paired {input.fq1} {input.fq2} --fastqc --gzip --clip_R1 3 --clip_R2 4 --three_prime_clip_R1 3 --three_prime_clip_R2 3 -o ./trimmed_reads
            """

elif config["trimming"] == "soft":
    rule adaptortrimming:
        input:
            fq1 = fastqfolder + "/{sample}_R1_001.fastq.gz",
            fq2 = fastqfolder + "/{sample}_R2_001.fastq.gz"
        output:
            temp("trimmed_reads/{sample}_R1_001_val_1.fq.gz"),
            temp("trimmed_reads/{sample}_R2_001_val_2.fq.gz")
        params:
            logfile = "{sample}_lengths.txt"
        container:
            "docker://quay.io/biocontainers/trim-galore:0.6.6--0"
        shell:
            """
            zcat {input.fq1} | awk "{{if(NR%4==2) print length(\$1)}}" | sort -n | uniq -c > {wildcards.sample}_lengths.txt
            python {workflow.basedir}/validate_fastqs.py {wildcards.sample}_lengths.txt
            trim_galore --paired {input.fq1} {input.fq2} --fastqc --clip_R1 3 --clip_R2 4 --three_prime_clip_R1 1 --three_prime_clip_R2 1 --rrbs --gzip -o ./trimmed_reads
            """

elif config["trimming"] == "none":
    rule adaptortrimming:
        input:
            fq1 = fastqfolder + "/{sample}_R1_001.fastq.gz",
            fq2 = fastqfolder + "/{sample}_R2_001.fastq.gz"
        output:
            temp("trimmed_reads/{sample}_R1_001_val_1.fq.gz"),
            temp("trimmed_reads/{sample}_R2_001_val_2.fq.gz")
        params:
            logfile = "{sample}_lengths.txt"
        container:
            "docker://quay.io/biocontainers/trim-galore:0.6.6--0"
        shell:
            """
            trim_galore --paired {input.fq1} {input.fq2} --fastqc --gzip -o ./trimmed_reads
            """

elif not config["trimming"]:
    sys.exit("Specify trimming (hard or soft or none)")


rule mapping_bismark:
    input:
        genome = config["humangenome"],
        fq1 = "trimmed_reads/{sample}_R1_001_val_1.fq.gz",
        fq2 = "trimmed_reads/{sample}_R2_001_val_2.fq.gz"
    output:
        bam = temp("mapped_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam"),
        dir = temp(directory("bismarkTemp_{sample}")),
    container:
        "docker://quay.io/biocontainers/bismark:0.24.0--hdfd78af_0"
    params:
        tempdir = "bismarkTemp_{sample}",
        id = "{sample}",
        threads = 4
    shell:
        "bismark --samtools_path /usr/local/bin/ --nucleotide_coverage --temp_dir {params.tempdir} --multicore {params.threads} {input.genome} -1 {input.fq1} -2 {input.fq2} -o ./mapped_reads/"

rule query_sort_samtools:
    input:
        bam = "mapped_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam"
    output:
        bam = temp("query_sorted_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam"),
    container:
        "docker://quay.io/biocontainers/samtools:1.14--hb421002_0"
    params:
        threads = 4
    shell:
        "samtools sort -n -@ {params.threads} {input.bam} -o {output.bam}; "

rule rm_optical_dups:
    input: "query_sorted_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam"
    output: temp("rm_optical_dups/{sample}_R1_001_val_1_bismark_bt2_pe.bam")
    container:
        "docker://quay.io/biocontainers/picard:2.21.6--0"
    params:
        OpticalDupsPixelDistance = OpticalDupsPixelDistance,
        ReadRegex = ReadRegex
    shell:
        "picard MarkDuplicates "
        "I={input} "
        "O={output} "
        "M=rm_optical_dups/{wildcards.sample}_MarkDuplicates.txt "
        "OPTICAL_DUPLICATE_PIXEL_DISTANCE={params.OpticalDupsPixelDistance} "
        "REMOVE_SEQUENCING_DUPLICATES=true "
        "TAGGING_POLICY=OpticalOnly "
        "READ_NAME_REGEX='{params.ReadRegex}' "
        "ASSUME_SORT_ORDER=queryname COMPRESSION_LEVEL=0; "

rule sort_samtools:
    input:
        bam = "rm_optical_dups/{sample}_R1_001_val_1_bismark_bt2_pe.bam"
    output:
        bam = "sorted_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam",
        bai = "sorted_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam.bai"
    container:
        "docker://quay.io/biocontainers/samtools:1.14--hb421002_0"
    params:
        threads = 4
    shell:
        "samtools sort -@ {params.threads} {input.bam} -o {output.bam}; "
        "samtools index {output.bam} {output.bai} "

rule rrbsmetrics:
    input:
        bam = "sorted_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam",
        bai = "sorted_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam.bai"
    output:
        txt = "logs/{sample}.rrbs_summary_metrics"
    container:
        "docker://quay.io/biocontainers/picard:2.21.6--0"
        #os.path.join(config["containers"],'rbase_3.6.2.sif')
    params:
        genomepath = config["humangenomefa"],
        out = "logs/{sample}"
    shell:
        "picard CollectRrbsMetrics "
        "I={input.bam} "
        "M={params.out}  "
        "R={params.genomepath}"

if config["genomebuild"] == "GRCh37" or config["genomebuild"] == "GRCh38" or config["genomebuild"] == "GRCm38" or config["genomebuild"] == "GRCm39":
    rule hsmetrics:
        input:
            bam = "sorted_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam",
            bai = "sorted_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam.bai"
        output:
            txt = "logs/{sample}_hs_metrics.txt"
        container:
            "docker://quay.io/biocontainers/picard:2.21.6--0"
        params:
            genomepath = config["humangenomefa"],
            targetpath = config["target"],
        shell:
            "picard CollectHsMetrics "
            "I={input.bam} "
            "O={output.txt} "
            "R={params.genomepath} "
            "TARGET_INTERVALS={params.targetpath} "
            "BAIT_INTERVALS={params.targetpath} "
elif not config["genomebuild"]:
    sys.exit("Specify genomebuild (GRCh37, GRCh38, GRCm38 or GRCm39)")

rule mapping_to_lambda_bismark:
    input:
        genome = config["lambdagenome"],
        fq1 = "trimmed_reads/{sample}_R1_001_val_1.fq.gz",
        fq2 = "trimmed_reads/{sample}_R2_001_val_2.fq.gz"
    output:
        bam = "mapped_reads/lambda.{sample}_R1_001_val_1_bismark_bt2_pe.bam",
        tmpfolder = temp(directory("bismarkTempLambda_{sample}"))
    container:
        "docker://quay.io/biocontainers/bismark:0.23.1--hdfd78af_0"
    params:
        tempdir = "bismarkTempLambda_{sample}",
        smallbam = "mapped_reads/lambda.{sample}_subsample.bam",
        threads = 2
    shell:
        "bismark --samtools_path /usr/local/bin/ --nucleotide_coverage --prefix lambda --temp_dir {params.tempdir} --multicore {params.threads} {input.genome} -1 {input.fq1} -2 {input.fq2} -o ./mapped_reads/ ;"

rule extract_methylation_calls:
    input:
        "rm_optical_dups/{sample}_R1_001_val_1_bismark_bt2_pe.bam"
    output:
        "methylation_extractor_output/{sample}_R1_001_val_1_bismark_bt2_pe.bismark.cov.gz"
    container:
        "docker://quay.io/biocontainers/bismark:0.23.1--hdfd78af_0"
    params:
        threads = 6,
        zerocov = zerocov
    shell:
        "bismark_methylation_extractor --samtools_path /usr/local/bin/ -p {params.zerocov} --bedGraph --gzip --multicore {params.threads} {input} -o ./methylation_extractor_output/"

rule deduplicate:
    input:
        "rm_optical_dups/{sample}_R1_001_val_1_bismark_bt2_pe.bam"
    output:
        temp("dedupl_reads/{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated.bam")
    container:
        "docker://quay.io/biocontainers/bismark:0.23.1--hdfd78af_0"
    shell:
        "deduplicate_bismark --samtools_path /usr/local/bin/ -p --bam {input} --output_dir ./dedupl_reads"

rule extract_methylation_calls_deduplicated:
    input:
        "dedupl_reads/{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated.bam"
    output:
        "methylation_extractor_output/{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
    container:
        "docker://quay.io/biocontainers/bismark:0.23.1--hdfd78af_0"
    params:
        threads = 5
    shell:
        "bismark_methylation_extractor --samtools_path /usr/local/bin/ --bedGraph --gzip --multicore {params.threads} {input} -o ./methylation_extractor_output/"

rule rename:
    input:
        expand("methylation_extractor_output/lambda.{sample}_R1_001_val_1_bismark_bt2_pe.bismark.cov.gz", sample=R1IDS),
        expand("methylation_extractor_output/{sample}_R1_001_val_1_bismark_bt2_pe.{dupl}bismark.cov.gz", sample=R1IDS, dupl = ['deduplicated.', ''])
    output:
        expand("methylation_extractor_output/{sample}{dupl}.cov.gz", sample=R1IDS, dupl = ['dedupl', '']),
        expand("methylation_extractor_output/lambda.{sample}.cov.gz", sample=R1IDS)
    shell:
        "cd methylation_extractor_output/ ;"
        "rename _R1_001_val_1_bismark_bt2_pe.bismark.cov.gz .cov.gz *_R1_001_val_1_bismark_bt2_pe.bismark.cov.gz; "
        "rename _R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz dedupl.cov.gz *_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz;"
        "rm *txt.gz"

if config["genomebuild"] == "GRCh37" or config["genomebuild"] == "GRCh38" or config["genomebuild"] == "GRCm38" or config["genomebuild"] == "GRCm39":
    rule multiqc:
        input:
            expand("logs/{sample}_hs_metrics.txt", sample = R1IDS),
            expand("logs/{sample}.rrbs_summary_metrics" , sample = R1IDS),
            expand("methylation_extractor_output/{sample}{dupl}.cov.gz", sample=R1IDS, dupl = ['dedupl', '']),
            expand("methylation_extractor_output/lambda.{sample}.cov.gz", sample=R1IDS),
        output:
            report("multiqc_report.html")
        container:
            "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
        shell:
            "multiqc -f . ;"
elif not config["genomebuild"]:
    sys.exit("Specify genomebuild (GRCh37, GRCh38, GRCm38 or GRCm39)")
