#!/bin/sh

# Run with: sbatch --time=24:00:00 rna_submit_snakemake.sh

#####################################
# This script is the RNAseq submit script for a snakemake pipeline.
# This pipeline was created by Matthew J Brooks in November 2020
# The pipeline was based off the NNRL lab pipeline.
# This pipeline adapted to run on HPCs running SLURM
# This requires the snakefile RNAseq_v3.0.py and rna_config.json
#####################################

# Load module
module load python/3.7

# Export variables
#NOW=$(date +"%Y%m%d")
NOW='20210830'
export BASE_DIR="/data/brooksma/RNA-seq_outside/Liang_jun"
export WORK_DIR=${BASE_DIR}/${NOW}
SNAKEFILE=${BASE_DIR}/RNAseq_v3.0.py

# Make result directories and change into result directory
mkdir -p ${WORK_DIR}/logs
cd $WORK_DIR

# Snakemake command
echo "Get ready for snakemake..." >> logs/snakemake.%j.o
snakemake\
	--directory $WORK_DIR \
	--snakefile $SNAKEFILE \
	--jobname '{rulename}.{jobid}' \
	--rerun-incomplete \
	--nolock \
	--verbose \
	-k -p \
	-j 3000 \
	--stats rna_pipeline_${NOW}.stats \
	--cluster "sbatch --mail-type=FAIL -o logs/{params.rulename}.%j.o {params.batch}" \
	>& rna_pipeline_${NOW}.log

# Summary
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --summary

## DRY Run with Print out the shell commands that will be executed
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --dryrun -p -r
# snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --dryrun -p -r

#DAG
 # snakemake --directory $WORK_DIR --snakefile $SNAKEFILE  --dag | dot -Tpng > dag.png

#Rulegraph
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE  -n --forceall --rulegraph | dot -Tpng > rulegraph.png

# Mail Rulegraph and DAG to self
#  echo DAG |mutt -s "DAG" -a dag.png -a rulegraph.png -- brooksma@mail.nih.gov
