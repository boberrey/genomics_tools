# Configuration file 

import os
import sys
import glob

# Input data

FASTQ_DIR = 'input/fastqs'
METADATA_FILE = 'input/metadata.txt'
BEDS = {"TSS": '/raid/shr/Downloaded_data/mm10_data/TSS/refSeqmm10.TSS.bed6.merged.bed'}
#BEDS = {"TSS" : '/lab/greenleaf/snakeATAC/resources/sacCer3/sgd_tss.bed', \
#"CACGTG" : '/lab/greenleaf/snakeATAC/resources/sacCer3/cacgtg_sites.bed', \
#"NUC" : '/lab/greenleaf/snakeATAC/resources/sacCer3/nuc_positions.eab.bed',\
#"Reb1" : '/lab/greenleaf/snakeATAC/resources/sacCer3/reb1_organic.10m.80mM.bed',\
#"Abf1" : '/lab/greenleaf/snakeATAC/resources/sacCer3/abf1_organic.2p5m.80mM.bed' }
NARROWPEAKS = {}
BROADPEAKS = {}
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
ATAC_TOOLS = "."



# metadata file

def make_meta(filename):
    """
    Generate a metadata file with the names of all the samples to be processed.
    Sample names are inferred from fastq files.
    """
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    r1_files = list(map(os.path.abspath, glob.glob(os.path.join(FASTQ_DIR,"*_R1*.fastq*"))))
    if (len(r1_files) < 1):
        sys.exit("No fastqs with _R1 found.")
    r2_files = [os.path.join(os.path.dirname(r1_file), 
        os.path.basename(r1_file).replace('R1', 'R2')) for r1_file in r1_files]
    if  all([os.path.isfile(r2_file) for r2_file in r2_files]) is False:
        sys.exit("Not all matching _R2 files found.")
    sample_labels = [os.path.basename(r1_file).split("_R1")[0] for r1_file in r1_files]
    with open(filename, 'w') as outfile:
        outfile.write("\t".join(["Name","Read1","Read2"]) + "\n")
        for sample_label, r1_file, r2_file in zip(sample_labels, r1_files, r2_files):
            if len(sample_label) > 30:
                sample_label = sample_label[:20] + "..." + sample_label[-10:]
            outfile.write("\t".join([sample_label, r1_file, r2_file]) + "\n")
