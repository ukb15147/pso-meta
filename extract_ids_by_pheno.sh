# #!/bin/bash -l
# #SBATCH --output=/scratch/users/%u/%j.out

## Submit using sbatch -p shared script.sh


V2_LIST=$1
CTV3_LIST=$2
OUT=$3

dos2unix ${V2_LIST}
dos2unix ${CTV3_LIST}

awk -F"\t" 'F==1{a2[$1]=$1}F==2{a3[$1]=$1}F==3{if($4 in a2 || $5 in a3){print $1}}' \
 F=1 ${V2_LIST} F=2 ${CTV3_LIST} \
 F=3 [path_on_hpc]/gp_clinical.txt \
 > ${OUT}

sort -u ${OUT} > ${OUT}.uniq

echo "Remember to exclude withdrawn participants. See directory derm_ukb/withdrawals for latest list"
