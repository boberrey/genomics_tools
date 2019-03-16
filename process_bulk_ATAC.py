#!/usr/bin/env python

"""
Perform basic essential steps of ATAC-seq data processing

Note: Python 3

Ben Ober-Reynolds, adapted from Evan Boyle
Stanford University
"""

import pandas as pd
import os
import sys
import glob
#import json


include: "process_bulk_ATAC_config.py" 


# Make metadata file if it doesn't exist
if not os.path.exists(METADATA_FILE): make_meta(METADATA_FILE)

# Add execution directory to path if it isn't already there
if EXE_DIR not in sys.path: os.environ["PATH"] = EXE_DIR + os.pathsep + os.environ["PATH"]

# Load metadata
metadata = pd.read_table(METADATA_FILE, index_col = False)

sample_labels = metadata.Name.tolist()

# What genome are we working with
genome_name = "hg38"

rule all:
    input:
        # Per sample output files
        # These are listed in the order generated
        expand("output/fastqs/qc/{sample_label}_R1_untrimmed_fastqc.html", sample_label = sample_labels),
        expand("output/fastqs/qc/{sample_label}_R1_trimmed_fastqc.html", sample_label = sample_labels),
        expand("output/bams/deduped/{sample_label}.noMT.filtered.deduped.bam", sample_label = sample_labels),
        expand("output/bams/qc/counts/{sample_label}.counts.txt", sample_label = sample_labels),
        expand("output/beds/{sample_label}.insertions.bed.gz", sample_label = sample_labels),
        expand("output/peaks/{sample_label}_summits.bed", sample_label = sample_labels),
        expand("output/coverage_data/{sample_label}.insertions.bw", sample_label = sample_labels),
        #expand("output/beds/{sample_label}.insertions.bedGraph", sample_label = sample_labels),
        #expand("output/plots/TSS/{sample_label}_TSS_enrichment.pdf", sample_label = sample_labels),

        # per sample plots:
        expand("output/plots/qc/insert_size/{sample_label}_insert_size_histogram.pdf",sample_label = sample_labels),

        # Pooled stats, etc.
        "output/bams/qc/compiled_flagstats.txt",
        "output/bams/qc/compiled_idxstats.txt",
        "output/bams/qc/compiled_idxstats.mito_fraction.txt",
        "output/bams/qc/compiled_picard.dedup_metrics.txt",
        "output/plots/TSS/TSS_enrichments.txt",
        # Pooled plots:
        "output/plots/qc/compiled_flagstats.pdf",
        "output/plots/qc/compiled_idxstats.mito.pdf",

        # Whole experiment files
        "output/peaks/processed/filtered_fix_width_peaks.bed",
        "output/counts_matrix/counts_table.txt"

    output:
        "snakeATAC.txt"
    shell:
        "echo $(date) > {output};"
        "echo snake make stuff"


"""
Trim Nextera adapters using Skewer 
Version 0.2.2
"""

rule trim_adapters_skewer:
    input:
        left = lambda wildcards: glob.glob(metadata.loc[metadata.Name == wildcards.sample_label]["Read1"].values[0]),
        right = lambda wildcards: glob.glob(metadata.loc[metadata.Name == wildcards.sample_label]["Read2"].values[0])
    output:
        temp_left_cat = temp("output/fastqs/{sample_label}_skewer_R1.fastq.gz"),
        temp_right_cat = temp("output/fastqs/{sample_label}_skewer_R2.fastq.gz"),
        # Evan had these here but they are not necessary and they cause issues. 
        # Listing files in 'output' but renaming them during the shell command
        # causes snakemake to freak out
        # 
        #temp_left = "output/fastqs/trimmed/{sample_label}-trimmed-pair1.fastq.gz",
        #temp_right = "output/fastqs/trimmed/{sample_label}-trimmed-pair2.fastq.gz",
        #temp_log = "output/fastqs/trimmed/{sample_label}-trimmed.log"
        left = temp("output/fastqs/trimmed/{sample_label}_R1_trimmed.fastq.gz"),
        right = temp("output/fastqs/trimmed/{sample_label}_R2_trimmed.fastq.gz"),
        log = "output/fastqs/qc/{sample_label}.skewer.log"
    params:
        error_out_file = "error_files/{sample_label}_trim",
        run_time = "2:30:00",
        cores = "1",
        memory = "6000",
        job_name = "trimming"

    benchmark:
        "benchmarks/trimming/{sample_label}.txt"
    threads: 1
    shell:
        "cat {input.left} > {output.temp_left_cat};" # if there are multiple files to be combined
        "cat {input.right} > {output.temp_right_cat};"
        "skewer \
        -x CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
        -y CTGTCTCTTATACACATCTGACGCTGCCGACGANNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
        -m pe -t {threads} -f sanger\
        {output.temp_left_cat} {output.temp_right_cat} \
        -o output/fastqs/trimmed/{wildcards.sample_label} -z;"
        "mv output/fastqs/trimmed/{wildcards.sample_label}-trimmed-pair1.fastq.gz {output.left};"
        "mv output/fastqs/trimmed/{wildcards.sample_label}-trimmed-pair2.fastq.gz {output.right};"
        "mv output/fastqs/trimmed/{wildcards.sample_label}-trimmed.log {output.log};"

"""
Trim Nextera adapters using SeqPurge
Version 2018_11
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1069-7
"""

rule trim_adapters_seqpurge:
    input:
        left = lambda wildcards: glob.glob(metadata.loc[metadata.Name == wildcards.sample_label]["Read1"].values[0]),
        right = lambda wildcards: glob.glob(metadata.loc[metadata.Name == wildcards.sample_label]["Read2"].values[0])
    output:
        left = temp("output/fastqs/trimmed/{sample_label}_R1_trimmed.fastq.gz"),
        right = temp("output/fastqs/trimmed/{sample_label}_R2_trimmed.fastq.gz")
    params:
        error_out_file = "error_files/{sample_label}_trim",
        run_time="2:30:00",
        cores="1",
        memory="6000",
        job_name="trimming"
    benchmark: 
        "benchmarks/trimming/{sample_label}.txt"
    threads: 1
    shell:
        "SeqPurge \
        -a1 CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
        -a2 CTGTCTCTTATACACATCTGACGCTGCCGACGA \
        -threads {threads} \
        -out1 {output.left} -out2 {output.right} \
        -in1 {input.left} -in2 {input.right}"

ruleorder:  trim_adapters_skewer > trim_adapters_seqpurge
#ruleorder: trim_adapters_seqpurge > trim_adapters_skewer


"""
Run fastQC on the trimmed and untrimmed fastqs to get some information about
potential problems with fastqs
"""
rule fastqc_unmapped_trimmed:
    input:
        left = "output/fastqs/trimmed/{sample_label}_R1_trimmed.fastq.gz",
        right = "output/fastqs/trimmed/{sample_label}_R2_trimmed.fastq.gz"
    output:
        "output/fastqs/qc/{sample_label}_R1_trimmed_fastqc.html",
        "output/fastqs/qc/{sample_label}_R2_trimmed_fastqc.html",
        # stuff we don't really care about but want to eliminate when run is botched
        "output/fastqs/qc/{sample_label}_R1_trimmed_fastqc.zip",
        "output/fastqs/qc/{sample_label}_R2_trimmed_fastqc.zip"
    params:
        error_out_file = "error_files/{sample_label}_trim_fastqc",
        run_time="00:15:00",
        cores="1",
        memory="6000",
        job_name="fastqc"
    benchmark: 
        "benchmarks/fastqc/{sample_label}_trim.txt"
    shell:
        "fastqc {input.left} {input.right} --outdir=" + "output/fastqs/qc/"


rule fastqc_unmapped_untrimmed:
    input:
        left = lambda wildcards: glob.glob(metadata.loc[metadata.Name == wildcards.sample_label]["Read1"].values[0]),
        right = lambda wildcards: glob.glob(metadata.loc[metadata.Name == wildcards.sample_label]["Read2"].values[0]),
    output:
        temp_left = temp("output/fastqs/qc/{sample_label}_R1_untrimmed.fastq.gz"),
        temp_right = temp("output/fastqs/qc/{sample_label}_R2_untrimmed.fastq.gz"),
        lh = "output/fastqs/qc/{sample_label}_R1_untrimmed_fastqc.html",
        rh = "output/fastqs/qc/{sample_label}_R2_untrimmed_fastqc.html",
        # stuff we don't really care about but want to eliminate when run is botched
        lz = "output/fastqs/qc/{sample_label}_R1_untrimmed_fastqc.zip",
        rz = "output/fastqs/qc/{sample_label}_R2_untrimmed_fastqc.zip"
    params:
        error_out_file = "error_files/{sample_label}_untrim_fastqc",
        run_time="00:15:00",
        cores="1",
        memory="6000",
        job_name="fastqc"
    benchmark: 
        "benchmarks/fastqc/{sample_label}_untrim.txt"
    run:
        shell("cat {input.left} > {output.temp_left}"),
        shell("cat {input.right} > {output.temp_right}"),
        shell("fastqc {output.temp_left} {output.temp_right} --outdir=output/fastqs/qc/;")


"""
Map trimmed reads using Bowtie2
Version 2.2.6
Excludes mates separated by more than 2000 bp
Sorts and indexes the bam file afterwards using samtools
For info on sam (sequence alignment map) and bam (binary of sam):
https://training.h3abionet.org/postgraduate_workshop_2014/wp-content/uploads/2014/04/H3ABioNet_2014_NGS_8_SamFormat.pdf
"""
rule run_bowtie:
    input:
        # Adding the '.1.bt2' is necessary for snakemake to recognize the file
        idx = REFERENCE_FILE + ".1.bt2",
        left = "output/fastqs/trimmed/{sample_label}_R1_trimmed.fastq.gz",
        right ="output/fastqs/trimmed/{sample_label}_R2_trimmed.fastq.gz"
    output:
        bam = temp("output/bams/unprocessed/{sample_label}.bam"),
        idx = temp("output/bams/unprocessed/{sample_label}.bam.bai")
    params:
        error_out_file = "error_files/{sample_label}_bowtie",
        run_time = "4:59:00",
        cores = "8",
        memory = "8000",
        job_name = "bwt2"
    benchmark: "benchmarks/bowtie/{sample_label}.txt"
    threads: 8
    shell: 
        # -X 2000 # prevents mates separated by a lot 
        "bowtie2 \
        -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50\
        -X 2000 --threads {threads} \
        --rg-id {wildcards.sample_label} \
        --rg 'SM:{wildcards.sample_label}' \
        -x " + REFERENCE_FILE + " -1 {input.left} -2 {input.right} \
        | samtools view -b -S - \
        | samtools sort -o output/bams/unprocessed/{wildcards.sample_label}.bam -; "
        "samtools index output/bams/unprocessed/{wildcards.sample_label}.bam; "


"""
Library complexity is estimated using preseq:
http://smithlabresearch.org/software/preseq/
"""
rule estimate_library_complexity:
    input:
        bam = rules.run_bowtie.output.bam
    output:
        lc = "output/bams/qc/complexity/{sample_label}.extrapolated_yield.txt",
        c = "output/bams/qc/complexity/{sample_label}.downsampled_yield.txt"
    params:
        error_out_file = "error_files/{sample_label}_estimate_lc",
        run_time = "1:00:00",
        cores = "1",
        memory = "8000",
        job_name = "lc_extrap"
    benchmark: "benchmarks/preseq/{sample_label}.txt"
    threads: 1
    shell:
        "preseq lc_extrap -P -o {output.lc} -B {input.bam}; " +
        "preseq c_curve -P -s 100000 -o {output.c} -B {input.bam}"


"""
flagstats calculated with SAMTools
"""
rule calc_flagstats:
    input:
        bam = "output/bams/unprocessed/{sample_label}.bam",
        idx = "output/bams/unprocessed/{sample_label}.bam.bai"
    output:
        "output/bams/qc/flagstats/{sample_label}.flagstat.txt" 
    params:
        error_out_file="error_files/flagstats",
        run_time="00:05:00",
        cores="1",
        memory="3000",
        job_name="flagstat"
    shell:
        "samtools flagstat {input.bam} | awk '{{print \"{wildcards.sample_label}\\t\" $0}}' > {output};"


"""
idxstats calculated with SAMTools
"""
rule calc_idxstats:
    input:
        bam = "output/bams/unprocessed/{sample_label}.bam",
        idx = "output/bams/unprocessed/{sample_label}.bam.bai"
    output:
        "output/bams/qc/idxstats/{sample_label}.idxstats.txt" 
    params:
        error_out_file="error_files/idxstats",
        run_time="00:05:00",
        cores="1",
        memory="1000",
        job_name="idxstats"
    shell:
        "samtools idxstats {input.bam} | awk '{{print \"{wildcards.sample_label}\\t\" $0}}' > {output};"


"""
compile SAMTools flagstats of all samples into one table.
"""
rule plot_flagstats:
    input:
        expand("output/bams/qc/flagstats/{sample_label}.flagstat.txt", sample_label=sample_labels)
    output:
        table = "output/bams/qc/compiled_flagstats.txt",
        pdf = "output/plots/qc/compiled_flagstats.pdf"
    params:
        error_out_file="error_files/flagstat_plot",
        run_time="00:10:00",
        cores="1",
        memory="1000",
        job_name="plot_flagstat"
    shell:
        (
            "awk 'BEGIN {{OFS = \"\\t\"; print \"sample_label\",\"total\",\"secondary\","
            "\"supplementary\",\"duplicates\",\"mapped\",\"paired\",\"read1\",\"read2\","
            "\"proper_pair\",\"both_mapped\",\"singletons\",\"separate_chr\",\"separate_chr_mapq_above5\"}} "
            "FNR == 1 && NR != 1 {{print \"\"}} FNR == 1 {{printf $1}} {{printf \"\\t\" $2 }} "
            "END {{print \"\"}} ' {input} > {output.table};"
            "Rscript --vanilla {ATAC_TOOLS}/qc_boxplot.R {output.table} read_count {output.pdf}"
        )


"""
compile SAMTools idxstats of all samples into one table.
"""
rule plot_idxstats:
    input:
        expand("output/bams/qc/idxstats/{sample_label}.idxstats.txt", sample_label=sample_labels)
    output:
        qc_table = "output/bams/qc/compiled_idxstats.txt",
        mito_table = "output/bams/qc/compiled_idxstats.mito_fraction.txt",
        qc_pdf = "output/plots/qc/counts/compiled_idxstats.counts.pdf",
        mito_pdf = "output/plots/qc/compiled_idxstats.mito.pdf",
    params:
        error_out_file="error_files/idxstats_plot",
        run_time="00:10:00",
        cores="1",
        memory="1000",
        job_name="plot_idxstats"
    shell:
        (
            "awk 'BEGIN {{OFS = \"\\t\"; print \"sample_label\",\"chr\",\"ref_length\",\"mapped\",\"unmapped\"}} "
            "{{totals[$1] += $4}} $2 == \"chrM\" {{mito[$1] += $4}} {{print}} END "
            "{{print \"sample_label\",\"total_reads\",\"mito_reads\",\"mito_percent\" "
            "> \"{output.mito_table}\"; for(s in totals) {{print s,totals[s],mito[s],mito[s]/totals[s] * 100 "
            "> \"{output.mito_table}\"}} }}' {input} > {output.qc_table};"
            "Rscript --vanilla {ATAC_TOOLS}/qc_bargraph.R {output.mito_table} sample_label mito_percent {output.mito_pdf};"
            "Rscript --vanilla {ATAC_TOOLS}/qc_bargraph.R {output.mito_table} sample_label total_reads {output.qc_pdf};"
        )


"""
Remove mitochondrial and chrY reads using samtools
"""
rule rm_mito:
    input:
        bam = rules.run_bowtie.output.bam,
        idx = rules.run_bowtie.output.idx
    output:
        bam = temp("output/bams/noMT/{sample_label}.noMT.bam"),
        idx = temp("output/bams/noMT/{sample_label}.noMT.bam.bai")
    params:
        error_out_file = "error_files/{sample_label}_remove_mitochondrial_reads",
        run_time = "00:30:00",
        cores = "1",
        memory = "4000",
        job_name = "rm_mt_reads"
    threads: 1
    shell:
        (
            "samtools idxstats {input.bam} | cut -f 1 | grep -v chrM | grep -v chrY | "
            "xargs samtools view -b {input.bam} > {output.bam}; " # something like this
            "samtools index {output.bam}"
        )


"""
Filter bams of low quality or unmapped reads
"""
rule filter_bams:
    input: 
        bam = rules.rm_mito.output.bam,
        idx = rules.rm_mito.output.idx
    output: 
        bam = temp("output/bams/filtered/{sample_label}.noMT.filtered.bam"),
        #idx = temp("output/bams/filtered/{sample_label}.noMT.filtered.bam.bai")
    params:
        error_out_file="error_files/{sample_label}_filtered_bams",
        mapq_threshold="20",
        run_time="00:30:00",
        cores="1",
        memory="8000",
        job_name="filter_bams"
    threads: 1
    run:
        # -F 1804: exclude flag, exludes unmapped, next segment unmapped, secondary alignments, not passing platform q, PCR or optical duplicates
        # See for explaination: https://broadinstitute.github.io/picard/explain-flags.html
        # -f 2: flags to require, properly aligned
        # -q 30: exlude low MAPQ, set as parameter to adjust
        if BLACKLIST is None:
            shell(
                "samtools view -F 1804 -f 2 -q {params.mapq_threshold} -b {input.bam} > {output.bam}"
                )
        else:
            shell(
                "samtools view -F 1804 -f 2 -q {params.mapq_threshold} -b {input.bam} \
                | bedtools intersect -v -abam - -b " + BLACKLIST + " -wa > {output}"
                )



# PICARD is not good at handling memory, if crashes rerun with more memory...
# include below if not standard illumnia naming convetion. Tries to extract read location information.
# All bam files listed as imput here to prevent deleting until final bams generated.
rule rm_duplicates_picard: 
    input: # low mapping quality reads can also be removed
        bam = rules.filter_bams.output.bam,
        #idx = rules.filter_bams.output.idx,
        rm_mito_bam = rules.rm_mito.output.bam,
        rm_mito_idx = rules.rm_mito.output.idx,
        filter_bam = rules.filter_bams.output.bam,
        #filter_idx = rules.filter_bams.output.idx,
        bowtie_bam = rules.run_bowtie.output.bam,
        bowtie_idx = rules.run_bowtie.output.idx
    output:
        bam = "output/bams/deduped/{sample_label}.noMT.filtered.deduped.bam",
        idx = "output/bams/deduped/{sample_label}.noMT.filtered.deduped.bam.bai",
        raw_metrics = "output/picard/duplicates/raw/picard_dedup_metrics_{sample_label}.txt",
        parsed_metrics = "output/picard/duplicates/parsed/picard_dedup_metrics_{sample_label}.parsed.txt", 
    params:
        error_out_file =  "error_files/{sample_label}_picard_rmdup",
        run_time="01:00:00",
        cores="1",
        memory="40000",
        job_name="picard_rm_duplicate_reads"
    benchmark: "benchmarks/picard_MarkDuplicates/{sample_label}.txt"
    threads: 4
    shell: # -Xms4g # this seems to get the process killed... # WE CAN INCLUDE READ_NAME INFO if we have illumina reads...
        # There's some issue with the metrics file right now. It's getting formatted strangely and
        # not recognizing the library names for some reason. 
        # Original:
        # "grep -A 1 ESTIMATED_LIBRARY_SIZE {output.raw_metrics} | tail -1 \
        #| cut -f 9-10 | xargs echo {wildcards.sample_label} | tr ' ' $'\t' > {output.parsed_metrics};"
        # New:
        # "grep -A 2 ESTIMATED_LIBRARY_SIZE {output.raw_metrics} | tail -1 \
        # | cut -f 8-9 | xargs echo {wildcards.sample_label} | tr ' ' $'\t' > {output.parsed_metrics};"
        (
            "java -XX:ParallelGCThreads=3 -jar {PICARD_JAR} MarkDuplicates "
            "INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.raw_metrics} "
            "REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT READ_NAME_REGEX=null; "
            "samtools index {output.bam}; " # and index
            "grep -A 2 ESTIMATED_LIBRARY_SIZE {output.raw_metrics} | tail -1 "
            "| cut -f 8-9 | xargs echo {wildcards.sample_label} | tr ' ' $'\t' > {output.parsed_metrics};"
        )


"""
Count the reads in each category (unprocessed, mito removed, filtered, deduped)
Put the counts into a single file for each sample.
"""
rule count_bam_reads:
    input:
        "output/bams/unprocessed/{sample_label}.bam",
        "output/bams/noMT/{sample_label}.noMT.bam",
        "output/bams/filtered/{sample_label}.noMT.filtered.bam",
        "output/bams/deduped/{sample_label}.noMT.filtered.deduped.bam",
    output:
        "output/bams/qc/counts/{sample_label}.counts.txt"
    params:
        error_out_file="error_files/bam_counts",
        run_time="00:10:00",
        cores="1",
        memory="3000",
        job_name="count_bam"
    shell:
        "u=$(samtools flagstat output/bams/unprocessed/{wildcards.sample_label}.bam | head -1 | cut -f 1 -d ' ');"
        "M=$(samtools flagstat output/bams/noMT/{wildcards.sample_label}.noMT.bam | head -1 | cut -f 1 -d ' ');"
        "f=$(samtools flagstat output/bams/filtered/{wildcards.sample_label}.noMT.filtered.bam | head -1 | cut -f 1 -d ' ');"
        "d=$(samtools flagstat output/bams/deduped/{wildcards.sample_label}.noMT.filtered.deduped.bam | head -1 | cut -f 1 -d ' ');"
        "echo {wildcards.sample_label} $'\t' $u $'\t' $M $'\t' $f $'\t' $d > output/bams/qc/counts/{wildcards.sample_label}.counts.txt"


"""
Compile some metrics on deduping and filtering
"""
rule plot_duplicate_stats:
    input:
        expand("output/picard/duplicates/parsed/picard_dedup_metrics_{sample_label}.parsed.txt", sample_label=sample_labels)
    output:
        duplicate_table = "output/bams/qc/compiled_picard.dedup_metrics.txt",
        est_libsize_pdf = "output/plots/qc/compiled_picard_rmdup.est_libsize.pdf",
        duplicate_percent_pdf = "output/plots/qc/compiled_picard_rm_dup.duplicate_percent.pdf"
    params:
        error_out_file="error_files/bam_count_plot",
        run_time="00:10:00",
        cores="1",
        memory="1000",
        job_name="plot_bam"
    shell:
        "awk 'BEGIN{{OFS=\"\\t\";print \"sample_label\",\"percent_duplication\",\"estimated_library_size\"}} {{print}}' {input} > {output.duplicate_table};"
        "Rscript --vanilla {ATAC_TOOLS}/qc_bargraph.R {output.duplicate_table} sample_label percent_duplication {output.duplicate_percent_pdf};"
        "Rscript --vanilla {ATAC_TOOLS}/qc_bargraph.R {output.duplicate_table} sample_label estimated_library_size {output.est_libsize_pdf};"


"""
Plot a bunch of bargraphs to visualize various aspects of filtering/deduping
"""
rule plot_bam_reads:
    input:
        expand("output/bams/qc/counts/{sample_label}.counts.txt", sample_label=sample_labels)
    output:
        fraction_table = "output/bams/qc/compiled_counts.fraction.txt",
        count_table = "output/bams/qc/compiled_counts.txt",
        disjoint_table = "output/bams/qc/compiled_counts.disjoint.txt",
        qc_fraction_pdf = "output/plots/qc/counts/compiled_counts.qc_fraction.pdf",
        total_count_pdf = "output/plots/qc/counts/compiled_counts.total.pdf",
        post_filter_count_pdf = "output/plots/qc/counts/compiled_counts.post_filter.pdf",
        post_dedup_count_pdf = "output/plots/qc/counts/compiled_counts.post_dedup.pdf",
        post_mito_count_pdf = "output/plots/qc/counts/compiled_counts.post_mitochondria.pdf",
        disjoint_count_pdf = "output/plots/qc/counts/compiled_counts.disjoint.pdf",
        mito_fraction_pdf = "output/plots/qc/counts/compiled_counts.mito_fraction.pdf",
        filter_fraction_pdf = "output/plots/qc/counts/compiled_counts.filter_fraction.pdf",
        duplicate_fraction_pdf = "output/plots/qc/counts/compiled_counts.duplicate_fraction.pdf"
    params:
        error_out_file="error_files/bam_count_plot",
        run_time="00:10:00",
        cores="1",
        memory="1000",
        job_name="plot_bam"
    shell:
        (
            "awk 'BEGIN{{OFS=\"\\t\";print \"sample_label\",\"total\",\"post_mitochondria\",\"post_filter\",\"post_dedup\"}} "
            "{{print $1,$2,$3,$4,$5}}' {input} > {output.count_table};"
            "awk 'BEGIN{{OFS=\"\\t\";print \"sample_label\",\"mitochondria\",\"filtered\",\"duplicate\",\"informative\"}} "
            "{{print $1,$2 - $3,$3 - $4,$4 - $5,$5}}' {input} > {output.disjoint_table};"
            "awk 'BEGIN{{OFS=\"\\t\";print \"sample_label\",\"mitochondria\",\"filtered\",\"duplicate\"}} "
            "{{print $1,($2-$3)/$2,($3-$4)/$3,($4-$5)/$4}}' {input} > {output.fraction_table};"
            "Rscript --vanilla {ATAC_TOOLS}/qc_boxplot.R {output.fraction_table} fraction_reads_removed {output.qc_fraction_pdf};"
            "Rscript --vanilla {ATAC_TOOLS}/qc_bargraph.R {output.count_table} sample_label total {output.total_count_pdf};"
            "Rscript --vanilla {ATAC_TOOLS}/qc_bargraph.R {output.count_table} sample_label post_filter {output.post_filter_count_pdf};"
            "Rscript --vanilla {ATAC_TOOLS}/qc_bargraph.R {output.count_table} sample_label post_dedup {output.post_dedup_count_pdf};"
            "Rscript --vanilla {ATAC_TOOLS}/qc_bargraph.R {output.count_table} sample_label post_mitochondria {output.post_mito_count_pdf};"
            "Rscript --vanilla {ATAC_TOOLS}/qc_bargraph.R {output.fraction_table} sample_label filtered {output.filter_fraction_pdf};"
            "Rscript --vanilla {ATAC_TOOLS}/qc_bargraph.R {output.fraction_table} sample_label duplicate {output.duplicate_fraction_pdf};"
            "Rscript --vanilla {ATAC_TOOLS}/qc_bargraph.R {output.fraction_table} sample_label mitochondria {output.mito_fraction_pdf};"
            "Rscript --vanilla {ATAC_TOOLS}/qc_stackbargraph.R {output.disjoint_table} sample_label {output.disjoint_count_pdf};"
        )


"""
Plot the insert size distribution of each sample
"""
rule plot_insert_size:
    input:
        bam = rules.rm_duplicates_picard.output.bam,
        idx = rules.rm_duplicates_picard.output.idx
    output:
        histogram_plot = "output/plots/qc/insert_size/{sample_label}_insert_size_histogram.pdf",
    params:
        error_out_file="error_files/{sample_label}_insert_size_hist",
        run_time="20:00:00",
        cores="1",
        memory="6000",
        job_name="plot_insert_size",
        insert_window="600"
    threads: 1
    shell:
        (
            "Rscript --vanilla {ATAC_TOOLS}/insert_size_dist_plot.R {input.bam} {output.histogram_plot} {params.insert_window}"
        )


"""
Create an insertion bed file for each sample. 
Also adjusts the insertion site based on how the transposase sits on the DNA
"""
rule make_insertion_bed:
    input:
        bam = rules.rm_duplicates_picard.output.bam,
        idx = rules.rm_duplicates_picard.output.idx
    output:
        bed = "output/beds/{sample_label}.insertions.bed.gz"
    params:
        error_out_file="error_files/{sample_label}_bam2bed",
        run_time="20:00:00",
        cores="1",
        memory="6000",
        job_name="bam2bed"
    threads: 1
    benchmark: "benchmarks/make_bed/{sample_label}_bam2bed.txt"
    shell:
        # Many people adjust the insertion start and end to reflect the
        # actual position of the transposase during insertion of the adpaters
        # That's what's going on in the awk commands seen here.
        (
            "bedtools bamtobed -i {input.bam} | awk -F $\'\\t\' 'BEGIN {{OFS = FS}} "
            "$6 == \"+\" {{$2 = $2 + 4; $3 = $2 + 1; print}} $6 == \"-\" "
            "{{$3 = $3 - 5; $2 = $3 - 1; print $0}}' "
            "| sort -k1,1 -k2,2n | gzip -c > {output.bed};"
        )  
  



"""
Call peaks on Tn5-corrected insertion positions.
This pipeline calls peaks on each smaple individually, and then merges peaks
downstream. 
"""
rule run_MACS2:
    input:
        bed = rules.make_insertion_bed.output.bed
    output:
        narrowPeak = temp("output/peaks/{sample_label}_peaks.narrowPeak"),
        peak_xls = temp("output/peaks/{sample_label}_peaks.xls"),
        peak_bed = "output/peaks/{sample_label}_summits.bed"
    params:
        error_out_file = "error_files/{sample_label}_MACS2_bam",
        run_time = "00:59:59",
        cores = "1",
        memory = "8000",
        job_name = "macs2"
    benchmark: "benchmarks/macs2/{sample_label}.bed.txt"
    shell:
        (
            # Activate the python2 environemnt to run macs2
            "PS1=''; source " + P2_ACTIVATE + ";"+
            "macs2 callpeak -g {EFFECTIVE_GENOME_SIZE} --name {wildcards.sample_label} "
            "--treatment {input.bed} --outdir output/peaks --format BED --shift -75 "
            "--extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01;"
            "source " + P3_ACTIVATE
        )



"""
Generate fixed width peaks as described in Jeff's TCGA ATAC paper
"""
rule process_peaks:
    input:
        expand("output/peaks/{sample_label}_summits.bed", sample_label=sample_labels)
    output:
        merged_peaks = temp("output/peaks/processed/merged_fix_width_peaks.bed"),
        filtered_peaks = "output/peaks/processed/filtered_fix_width_peaks.bed"
    params:
        error_out_file="error_files/process_peaks",
        run_time="00:30:00",
        cores="1",
        memory="20000",
        job_name="fix_widths"
    benchmark: "benchmarks/process_peaks/process_peaks.txt"
    threads: 10
    shell:
        (
            "Rscript --vanilla {ATAC_TOOLS}/get_fixed_width_peaks.R output/peaks output/peaks/processed mm10 5 250 10;"
            "Rscript --vanilla {ATAC_TOOLS}/filter_peakset.R {output.merged_peaks} {output.filtered_peaks}"
        )


"""
Generate a counts matrix using the filtered, merged peak set
"""
# get_counts_matrix.R <peak_bed_file> <bam_file_dir> <output_filename> <n_cores> <genome> <normalization>
rule make_counts_matrix:
    input:
        peak_bed = "output/peaks/processed/filtered_fix_width_peaks.bed",
    output:
        counts_matrix = "output/counts_matrix/counts_table.txt"
    params:
        error_out_file="error_files/raw_counts_table",
        run_time="00:10:00",
        cores="10",
        memory="1000",
        job_name="plot_bam",
        genome=genome_name,
        normalization="housekeeping",
        bam_dir = "output/bams/deduped/"
    threads: 10
    shell:
        (
            "Rscript --vanilla {ATAC_TOOLS}/get_counts_matrix.R "
            "{input.peak_bed} {params.bam_dir} {output.counts_matrix} 10 {params.genome} {params.normalization}"
        )


"""
Make an insertion bigWig file for each sample
bigWigs will be made from final bam files scaled to 3e7 reads in peaks
"""
# generate_ATAC_signal_tracks.R <bam_file> <peaks_file> <output_filename> <genome> <bin_size> <normalization_style> <scale_factor_file>
rule make_insertion_bw:
    input:
        #bed = "output/beds/{sample_label}.insertions.bed.gz"
        bam = rules.rm_duplicates_picard.output.bam,
        idx = rules.rm_duplicates_picard.output.idx,
        peaks = rules.process_peaks.output.filtered_peaks,
        counts = rules.make_counts_matrix.output.counts_matrix
    output:
        bw = "output/coverage_data/{sample_label}.insertions.bw"
    params:
        error_out_file="error_files/{sample_label}_insertion_bw",
        run_time="00:30:00",
        cores="1",
        memory="2000",
        job_name="insertion_bw",
        genome=genome_name,
        bin_size='100',
        normalization="housekeeping",
        scale_factors="output/counts_matrix/scale_factors.txt"
    benchmark: "benchmarks/insertion_bw/{sample_label}.txt"
    threads: 1
    shell:
        (
            "Rscript --vanilla {ATAC_TOOLS}/generate_ATAC_signal_tracks.R {input.bam} {input.peaks} {output.bw} "
            "{params.genome} {params.bin_size} {params.normalization} {params.scale_factors}"
        )


"""
Make TSS insertion enrichment plots
"""
rule plot_insertion_profile:
    input:
        expand("output/bams/deduped/{sample_label}.noMT.filtered.deduped.bam", sample_label=sample_labels)
    output:
        #"output/plots/TSS/{sample_label}_TSS_enrichment.pdf",
        "output/plots/TSS/TSS_enrichments.txt"
    params:
        error_out_file="error_files/insertion_profile",
        run_time="00:30:00",
        cores="10",
        memory="6000",
        job_name="insertion_profile",
        genome=genome_name,
        window="2000",
        ncores="10"
    benchmark: "benchmarks/insertion_profiles/insertion_profile.txt"
    threads: 10
    shell:
        (
            "Rscript --vanilla {ATAC_TOOLS}/TSS_enrichment_plot.R output/bams/deduped output/plots/TSS "
            "{params.genome} {params.window} {params.ncores}"
        )
        





