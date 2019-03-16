# Configuration file 

import os
import sys
import glob
import re

# Input data

FASTQ_DIR = 'input/fastqs'
METADATA_FILE = 'input/metadata.txt'
#PCA_COLOR = "Agent"
#PCA_SHAPE = "None"
#CONTRAST_FILE = "input/contrasts.txt"

# Resources
REFERENCE_FILE = '/raid/shr/genomes/bowtie2/mm10/mm10'
P2_ACTIVATE = '/home/ben/venvs/py2k/bin/activate'
P3_ACTIVATE = '/home/ben/venvs/py2k/bin/activate'
CHROM_SIZES = '/raid/shr/gSizes/mm10.all.genomsize'
EFFECTIVE_GENOME_SIZE = 1.87*10**9 #hg19 = 2.7e9, hg38 = 3.05e9, mm9 = 1.87e9, sacCer3 = 1.2*10**7
#FASTQ_SCREEN_CONF = "/app/fastqscreen/fastq_screen_v0.4.4/fastq_screen.conf"
BLACKLIST = '/raid/shr/Downloaded_data/mm10_data/blacklist_kundaje/mm10-blacklist-liftedover-from-mm9.bed'

# this directory is added to the system path for the purposes of running snakeATAC
EXE_DIR = '/raid/app/ATAC_bin'

PICARD_JAR = '/raid/app/ATAC_bin/picard.jar'
#SNAKE_DIR = '/mnt/raid1/lab/greenleaf/snakeATAC'
#ATAC_TOOLS = SNAKE_DIR + '/atac_tools'
ATAC_TOOLS = "/home/ben/git_clones/ATAC_tools"
CHIP_INPUT_HEADER = "mESC_WCE_SRR5077675"


# metadata file

def make_meta(filename):
    """
    Generate a metadata file with the names of all the samples to be processed.
    Sample names are inferred from fastq files.
    Assumes single read files with SRA naming conventions (i.e. ends with _1.fastq*)
    """
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    # R1 files may have either *_R1*.fastq* or just *_1*.fastq* (This is not very specific unfortunately)
    #r1_files = list(map(os.path.abspath, re.findall(r"*_R?1\.fastq*")))
    r1_files = list(map(os.path.abspath, glob.glob(os.path.join(FASTQ_DIR,"*_1.fastq*"))))
    print(r1_files)
    if (len(r1_files) < 1):
        sys.exit("No fastqs with _R1 or _1.fastq found.")
    sample_labels = [os.path.basename(r1_file).split("_1")[0] for r1_file in r1_files]
    with open(filename, 'w') as outfile:
        outfile.write("\t".join(["Name","Read1"]) + "\n")
        for sample_label, r1_file in zip(sample_labels, r1_files):
            if len(sample_label) > 30:
                sample_label = sample_label[:20] + "..." + sample_label[-10:]
            outfile.write("\t".join([sample_label, r1_file]) + "\n")
