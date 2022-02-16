#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -o cluster_logs/slurm-%x-%j-%N.out

umask 002

BATCH=$1
BIOSAMPLES=$2
CONFIG=$3

if [ -z $3 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
 echo -e "\nUsage: $(basename $0) <batch_name> <biosample_csv> <config_yaml>\n"
 exit 0
fi

mkdir -p "batches/${BATCH}/"
LOCKFILE="batches/${BATCH}/process_batch.lock"

# add lockfile to directory to prevent multiple simultaneous jobs
lockfile -r 0 "${LOCKFILE}" || exit 1
trap "rm -f ${LOCKFILE}; exit" SIGINT SIGTERM ERR EXIT

# execute snakemake
snakemake --reason \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds \
    --configfile "${CONFIG}" \
    --config batch="${BATCH}" \
             biosamples="${BIOSAMPLES}" \
             scripts=workflow/scripts \
    --nolock \
    --local-cores 4 \
    --jobs 500 \
    --max-jobs-per-second 1 \
    --use-conda --conda-frontend mamba \
    --use-singularity --singularity-args '--nv ' \
    --latency-wait 120 \
    --cluster-config workflow/cluster.yaml \
    --cluster "sbatch --partition={cluster.partition} \
                      --cpus-per-task={cluster.cpus} \
                      --output={cluster.out} {cluster.extra} " \
    --snakefile workflow/Snakefile_capture
