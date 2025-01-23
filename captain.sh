#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=10-00:00:00
#SBATCH --parsable
#SBATCH -J "[[PIPENICKNAME]]"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output "[[WORKDIR]]/Pipe_runtime/[[SNAPSHOT]]/logs/Pipe_runtime.%j.out"
#SBATCH --error  "[[WORKDIR]]/Pipe_runtime/[[SNAPSHOT]]/logs/Pipe_runtime.%j.err"
set -euo pipefail

current_hostname=$(hostname)

# Check hostname to avoid running snakemake in the biowulf login node
# if [ "$current_hostname" == "biowulf.nih.gov" ]; then
#     echo "Current hostname is biowulf.nih.gov. Stopping the script."
#     exit 1
# else
#     module load snakemake
# fi
if [ "$current_hostname" == "lnode001" ]; then
    echo "Current hostname is login node. Stopping the script."
    exit 1
else
    # module load snakemake
    source activate base
    conda activate snake_env
fi


uid=$(uuidgen|cut -d '-' -f1)
echo 'Snakemake captain uid:' $uid
echo $uid > [[WORKDIR]]/Pipe_runtime/current_uid 
echo [[SNAPSHOT]] > [[WORKDIR]]/Pipe_runtime/current_snapshot

# Main process of pipeline
snakemake --latency-wait 120 --snakefile Snakefile -d "[[WORKDIR]]" \
  --cores 8 \
  --local-cores 8 \
  --jobs 999 \
  --latency-wait 120 all \
  --use-conda \
  --configfile="[[WORKDIR]]/Pipe_runtime/[[SNAPSHOT]]/config.yaml" \
  --printshellcmds \
  --jobname "{name}.{jobid}.uid_"$uid".sh" \
  --keep-going \
  --executor slurm \
  --restart-times 1 \
  --rerun-incomplete \
  --rerun-triggers mtime && \
snakemake -d "[[WORKDIR]]/Pipe_runtime/[[SNAPSHOT]]" --report "[[WORKDIR]]/Pipe_runtime/[[SNAPSHOT]]/Snakemake_Report.html" && \
echo Snakemake run $uid finished!
