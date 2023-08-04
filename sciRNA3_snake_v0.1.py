#######################################
# This is an sciRNA3 analysis snakemake script.
# Written my Matthew J. Brooks on August 3rd, 2023
# Modified from Diego snakefile of Junue sci3_RNA pipeline
# This runs on an HPC running SLURM
#######################################



#######################################
# Import config file and modules needed
#######################################

# Import modules
import glob
import os
import json
import pandas as pd

# Snakemake Base location
try:
	BASE_DIR=os.environ['BASE_DIR']
except KeyError:
	print("I can not locate your BASE_DIR directory.")
	pass

# Import configs
configfile: BASE_DIR + "/rna_config.json"

fq_dir = config["fq_dir"]
meta_file = config["meta_file"]
meta = pd.read_csv(config["meta_file"])
fq_R1_end = config["fastq_R1_end"]
fq_R2_end = config["fastq_R2_end"]

SAMPLES = meta.tolist()



#############################################################
# List of directories needed and end point files for analysis
#############################################################

# DIRS = ['trim/','BAMS/','star/','kallisto/','fastqc/','logs/']

RDS = expand(OUTPUT_PATH + "seurat_objects/{sample}/raw.Rds", sample=SAMPLES)
# FQC = expand("fastqc/{sample}{read}_fastqc.html", sample=SAMPLES, read=READS)
# KAL = expand("kallisto/{sample}/abundance.h5", sample=SAMPLES)
# IDX = expand("BAMS/{sample}.bam.bai", sample=SAMPLES)
# QCGB = ["QC/GeneBody/rseqc.Coverage.heatMap.pdf"]
# QCF = ["QC/Multi_FQ/Fastqc_multiqc_report.html"]
# QCS = ["QC/Multi_STAR/STAR_multiqc_report.html"]
# QCR = expand("count_rRNA/{sample}/rRNA.txt", sample=SAMPLES)

##############################
# Snakemake rules for analysis
##############################

localrules: all

rule all:
        input:  RDS #QCR + KAL + FQC + IDX + QCGB + QCF + QCS
        params:
                batch = config["job_all"]


rule UMI_attach:
    """
    # This rule attaches the UMI from read1 to the name in read2
    # The script takes in an input folder, a sample ID, an output folder,
    # an oligo-dT barcode file, a corresponding N5 barcode file, and
    # it pass the factors to the python script
    """
    input:
            r1 = fq_dir + "/{sample}" + fq_R1_end,
            r2 = fq_dir + "/{sample}" + fq_R2_end,
    output:
            OUTPUT_PATH + "UMI_attach/{sample}.R2.fastq.gz"
    log:    "logs/umi_attach.{sample}.log"
    params:
            rulename = "umi_attach",
            batch = config["job_umi"],
            script = config["script"],
            lig = config["lig_barcode"],
            oligo = config["oligo_barcode"],
            core = config["core"]

    script: "{params["script"]}"




# rule process_sample_step1:
#     input: r1=lambda w: SAMPLE_TO_PATH[w.sample][0],
#         r2=lambda w: SAMPLE_TO_PATH[w.sample][1]
#     output:
#         temp(OUTPUT_PATH + "UMI_attach/{sample}.R2.fastq.gz"),
#         temp(OUTPUT_PATH + "trimmed_fastq/{sample}.R2_trimmed.fq.gz"),
#         OUTPUT_PATH + "trimmed_fastq/{sample}.R2.fastq.gz_trimming_report.txt"
#     priority: 10
#     params:
#         error_out_file=BASE_PATH + "/slurm_files/RNA_S1_{sample}",
#         run_time="10:00:00", cores="1", memory="3", job_name="RNA_S1"
#     shell:
#         # This provides 1) fastq path, 2) sample id, 3) output folder path
#         # 4) GTF file, 5) STAR index file
#         "{BASE_PATH}scripts/human_mapping_step1.sh {input.r1} {input.r2} "
#         "{wildcards.sample} {OUTPUT_PATH} {HU_GTF} {HU_STAR_IDX} ; "












