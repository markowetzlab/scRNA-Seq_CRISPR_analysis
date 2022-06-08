#!/bin/bash
#SBATCH --partition general
#SBATCH --no-requeue
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --job-name jupyter-notebook
#SBATCH --output jupyter-notebook-%J.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexander.baker@cruk.cam.ac.uk

# get tunneling info
port=$(shuf -i8000-9999 -n1)
node=$(hostname -s)
user=$(whoami)
cluster=clust1-headnode-1.cri.camres.org

# print tunneling instructions jupyter-log
echo -e "
MacOS or linux terminal command to create your ssh tunnel
ssh -N -L ${port}:${node}:${port} baker02@clust1-headnode-1.cri.camres.org

Use a Browser on your local machine to go to:
localhost:${port}  (prefix w/ https:// if using password)
" > jupyter_notebook_info.txt

# load modules or conda environments here
source activate scRNA_env
jupyter lab --no-browser --port=${port} --ip=${node} --notebook-dir=../
