# #!/bin/bash -l
# #SBATCH -p shared,brc  # More resource than shared
# #SBATCH --output=/scratch/users/%u/000_sbatch_job_output/%j.out
# #SBATCH --time=5-0  # Maximum of 5 days (instead of default 1)
# #SBATCH --mem=10gb  # Mem of 10 GB (instead of default 1GB)
# #SBATCH --ntasks=4  # For 4 cores
# #SBATCH --nodes=1  # For use with --ntasks to keep cores on same node

REFLIST=cohorts_with_sample_sizes_vExample.txt
RESULTSVER=4

for i in $(seq 3 20); do
 COHORT=$(sed -n ${i}p ${REFLIST} | cut -f1)
 FILE=$(sed -n ${i}p ${REFLIST} | cut -f2)
 CASES=$(sed -n ${i}p ${REFLIST} | cut -f3)
 CTRLS=$(sed -n ${i}p ${REFLIST} | cut -f4)
 EFF=$(echo ${CASES} ${CTRLS} | awk '{print 4/((1/$1) + (1/$2))}')

 echo $i ${COHORT} ${EFF}

 if [[ "${FILE: -3}" == ".gz" ]]; then
  zcat ${FILE} |  awk -va=${EFF} '{if(NR==1){print $0,"N"}else{print $0,a}}' \
   > input_results/${COHORT}_results${RESULTSVER}
 else
  awk -va=${EFF} '{if(NR==1){print $0,"N"}else{print $0,a}}' ${FILE} \
   > input_results/${COHORT}_results${RESULTSVER}
 fi

done
