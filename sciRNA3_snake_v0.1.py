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
#try:
#	BASE_DIR=os.environ['BASE_DIR']
#except KeyError:
#	print("I can not locate your BASE_DIR directory.")
#	pass

# Import configs
configfile: "sciRNA3_config.json"

out_dir = config["out_dir"]
fq_dir = config["fq_dir"]
meta_file = config["meta_file"]
meta = pd.read_csv(meta_file)
fq_R1_end = config["fastq_R1_end"]
fq_R2_end = config["fastq_R2_end"]

SAMPLES = meta["sample"].tolist()



#############################################################
# List of directories needed and end point files for analysis
#############################################################

# DIRS = ['trim/','BAMS/','star/','kallisto/','fastqc/','logs/']

#RDS = expand(out_dir + "seurat_objects/{sample}/raw.Rds", sample=SAMPLES)
UMI = expand(out_dir + "/UMI_attach/{sample}.R2.fastq.gz", sample=SAMPLES)
TRM = expand(out_dir + "/trimmed_fastq/{sample}.R2_trimmed.fq.gz", sample=SAMPLES)
SMF = expand(out_dir + "/filtered_sam/{sample}.sam", sample=SAMPLES)
RD2 = expand(out_dir + "/rmdup_sam_2/{sample}.sam", sample=SAMPLES)
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
        input:  UMI + SMF + RD2 #QCR + KAL + FQC + IDX + QCGB + QCF + QCS
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
            r2 = fq_dir + "/{sample}" + fq_R2_end
    output:
            out_dir + "/UMI_attach/{sample}.R2.fastq.gz"
    log:    "logs/umi_attach.{sample}.log"
    params:
            rulename = "umi_attach",
            batch = config["job_umi"],
            out_dir = config["out_dir"],
            script = config["script_umi"],
            lig = config["lig_barcode"],
            oligo = config["oligo_barcode"],
            core = config["core"]

    shell: """
    module unload python
    source /data/brooksma/conda/etc/profile.d/conda.sh
    conda activate sciRNA3
    mkdir -p {params.out_dir}/UMI_attach
    python {params.script} {input.r1} {input.r2} {wildcards.sample} {params.out_dir}/UMI_attach {params.lig} {params.oligo} {params.core}
    """

rule trim_fq:
   """
   # This rule adapter trims the UMI corrected fastq file
   """
   input:
            out_dir + "/UMI_attach/{sample}.R2.fastq.gz"
   output:
            out_dir + "/trimmed_fastq/{sample}.R2_trimmed.fq.gz"
   log:    "logs/trim.{sample}.log"
   version: 
           config["trimgalore"]
   params:
           rulename = "trim",
           batch = config["job_trim"],
           out_dir = config["out_dir"],
   shell: """
   module load trimgalore/{version} || exit 1
   
   ##########
   ## Trim the UMI corrected fastq file
   ##########
   trim_galore {input} \
   -a AAAAAAAA \
   --three_prime_clip_R1 1 \
   -o {params.out_dir}/trimmed_fastq
   """


rule star:
   """
   # This rule aligns to the genome using STAR
   """
   input:
            out_dir + "/trimmed_fastq/{sample}.R2_trimmed.fq.gz"
   output:
            sam = out_dir + "/star/{sample}/{sample}Aligned.out.sam"
   log:    "logs/star.{sample}.log"
   version:
           config["star"]
   params:
           rulename = "star",
           batch = config["job_star"],
           out_dir = config["out_dir"],
           ref = config["staridx"]
   shell: """
   module load STAR/{version} || exit 1
   
   ##########
   ## Align with STAR
   ##########
   mkdir -p /lscratch/${{SLURM_JOB_ID}}/STARtmp
   mkdir -p {params.out_dir}/star/{wildcards.sample}
   STAR --runThreadN ${{SLURM_CPUS_ON_NODE}} \
   --genomeDir {params.ref} \
   --readFilesIn {input} \
   --outFileNamePrefix {params.out_dir}/star/{wildcards.sample}/{wildcards.sample} \
   --outTmpDir=/lscratch/${{SLURM_JOB_ID}}/STARtmp/{wildcards.sample} \
   --readFilesCommand zcat \
   --outSAMstrandField intronMotif
   """


rule sam_filt:
   """
   # Filters alignment sam file with samtools
   """
   input:
            out_dir + "/star/{sample}/{sample}Aligned.out.sam"
   output:
           out_dir + "/filtered_sam/{sample}.sam"
   log:    "logs/sam_filt.{sample}.log"
   version:
           config["samtools"]
   params:
           rulename = "sam_filt",
           batch = config["job_samfilt"],
           out_dir = config["out_dir"]
   shell: """
   module load samtools/{version} || exit 1
   mkdir -p {params.out_dir}/filtered_sam
   samtools view \
   -bh -q 30 -F 4 \
   {input} | \
   samtools sort -@ ${{SLURM_CPUS_ON_NODE}} - | \
   samtools view -h - > \
   {output}
   """


rule rmdup1:
    """
    # This rule accepts as input a sorted sam file, an output sam file, and a mismatch rate.
    # It will then remove duplicates based on the barcode + UMI (edit distance <= 1),
    # and chromatin and start site, at the same time, it will output the duplication
    # number for each read.
    """
    input:
            out_dir + "/filtered_sam/{sample}.sam"
    output:
            sam = out_dir + "/rmdup_sam/{sample}.sam",
            csv = out_dir + "/rmdup_sam/{sample}.sam.csv",
            csv2 = out_dir + "/report/duplicate_read/{sample}.sam.csv"
    log:    "logs/rmdup1.{sample}.log"
    params:
            rulename = "rmdup1",
            batch = config["job_rmdup"],
            out_dir = config["out_dir"],
            script = config["script_rmdup"]

    shell: """
    ##########
    ## Run rmdup script
    mkdir -p {params.out_dir}/rmdup_sam
    python {params.script} {input} {output.sam} 0 \
    ##########
    ## Move csv output to report directory
    mkdir -p {params.out_dir}/report/duplicate_read
    cp {output.csv} {params.out_dir}/report/duplicate_read \
    """


rule rmdup2:
    """
    # This rule accepts as input a sam file previously without perfect match duplicates, an output sam file, and a mismatch rate.
    # It will then remove duplicates based on the barcode + UMI (edit distance <= 1),
    # and chromatin and start site, at the same time, it will output the duplication
    # number for each read.
    """
    input:
            out_dirsam = out_dir + "/rmdup_sam/{sample}.sam"
    output:
            sam = out_dir + "/rmdup_sam_2/{sample}.sam",
            csv = out_dir + "/rmdup_sam_2/{sample}.sam.csv"
    log:    "logs/rmdup2.{sample}.log"
    params:
            rulename = "rmdup2",
            batch = config["job_rmdup"],
            out_dir = config["out_dir"],
            script = config["script_rmdup"]

    shell: """
    ##########
    ## Run rmdup script
    mkdir -p {params.out_dir}/rmdup_sam_2
    python {params.script} {input} {output.sam} 1 \
    """


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












