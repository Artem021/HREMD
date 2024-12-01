#!/usr/bin/sh
#SBATCH -D /s/ls4/users/artem_k/aREMD_xrd/HREMD
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=aREMD_job
#SBATCH --partition=hpc4-el7-3d


export PLUMED_NUM_THREADS=8
echo $PATH
export PATH=/s/ls4/users/artem_k/miniconda3/bin:/s/ls4/users/artem_k/miniconda3/condabin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin
echo $PATH
/usr/bin/apptainer exec /home/users/artem_k/lammps_container/lammps_plumed_py.sif /usr/local/bin/miniconda3/bin/python main.py
# apptainer exec /home/users/artem_k/lammps_container/lammps_plumed_py.sif /usr/local/bin/miniconda3/bin/python main.py


# echo $(which python)
# echo $(which conda)


# export PATH=/home/users/artem_k/miniconda3/bin:/home/users/artem_k/miniconda3/condabin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin
# export PATH=/s/ls4/users/artem_k/miniconda3/bin:/s/ls4/users/artem_k/miniconda3/condabin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin
# export LD_LIBRARY_PATH=/s/ls4/users/artem_k/miniconda3/lib
# export LIBRARY_PATH=/s/ls4/users/artem_k/miniconda3/lib


# source /s/ls4/users/artem_k/.bashrc
# source /s/ls4/users/artem_k/.bash_profile

# echo $(which python)
# echo $(which conda)


# echo $PLUMED_NUM_THREADS
# export PLUMED_NUM_THREADS=8
# echo $PLUMED_NUM_THREADS

# python main.py
# /s/ls4/users/artem_k/miniconda3/bin/python main.py