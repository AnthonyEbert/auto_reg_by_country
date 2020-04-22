#!/bin/bash
# 
# script to generate and submit independent jobs to HPC cluster.
#
# Author: David Warne (david.warne@qut.edu.au)
#         School of Mathematical Sciences
#         Science and Engineering Faculty
#         Queensland University of Technology
#
# NOTE: must be modified for site specific HPC configurations
#
# NOTE: submits N independent single code PBS jobs, if smc_abc_rw_par is used
#       then set ncpus to a larger number.
#
# NOTE: Due to significant variability in runtimes between countries, the 
#       PBS array job option is not used.

# for every country
for (( j=1;j<=252;j++ ))
do
    # write PBS submit script
    cat > sub << EOF
#!/bin/bash -l
#PBS -N COV19R$j
#PBS -l ncpus=1
#PBS -l walltime=24:00:00
#PBS -l mem=4gb

module load matlab
cd \$PBS_O_WORKDIR
mkdir results_iso3166_R$j
cp ./*.m results_iso3166_R$j
cd results_iso3166_R$j
matlab -r "country_id=$j ; run_smc_intensity_reg_cluster; exit();"
cd ../
EOF
    # submit
    qsub sub
done
