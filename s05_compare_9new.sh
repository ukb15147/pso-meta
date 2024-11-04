#!/bin/env bash
#$ -cwd
#$ -q HighMemLongterm.q,LowMemLongterm.q,LowMemShortterm.q

## load any modules you need 

module load bioinformatics/plink2/1.90b3.38
module load bioinformatics/king/2.0


## Script to compare relatedness between datasets
## Note that the null marker list does not contain any strand ambiguous SNPs so no need to handle
## Three datasets so far do not contain any strand mismatches

## The "master list" file has two columns: study acronym ("UKB" for UK Biobank) and filestem for the
##   associated PLINK files containing shared variant subsets (e.g. for UK Biobank it lists
##   [path_on_hpc]/ukb_step4_interim_common and at this location there are two sets of PLINK files:
##   ukb_step4_interim_common1.bed/bim/fam and ukb_step4_interim_common2.bed/bim/fam - the two versions
##   correspond to variants shared with different meta-analysis input studies). This script uses the
##   "1" version and there is an analogous script for the "2" version


#mkdir comparisons

## Loop through each pair of datasets
for i in $(seq 1 $(($(cat master_list_9new_kielErlangenAdj.txt | wc -l) - 1)));
 do for j in $(seq $(($i + 1)) $(cat master_list_9new_kielErlangenAdj.txt | wc -l)); do

  echo PAIR $i $j

  COHORT_I=$(sed -n "${i}p" master_list_9new_kielErlangenAdj.txt | cut -f1)
  COHORT_J=$(sed -n "${j}p" master_list_9new_kielErlangenAdj.txt | cut -f1)

  echo $i is ${COHORT_I}
  echo $j is ${COHORT_J}

  if test -f comparisons/DONE_${COHORT_I}_${COHORT_J}; then
   echo COMPARISON ${COHORT_I}_${COHORT_J} ALREADY DONE
   continue
  fi

  DATA_I=$(sed -n "${i}p" master_list_9new_kielErlangenAdj.txt | cut -f2)
  DATA_J=$(sed -n "${j}p" master_list_9new_kielErlangenAdj.txt | cut -f2)

  ## Copy data locally but prefix fam files with cohort
  if test -f comparisons/${COHORT_I}.fam; then
   echo ${COHORT_I} DATA ALREADY IN PLACE
  else
   awk -vc=${COHORT_I} '{print c"_"$1,c"_"$2,$3,$4,$5,$6}' ${DATA_I}1.fam > comparisons/${COHORT_I}.fam
   awk '{print $1,$1"_"$4,$3,$4,$5,$6}' ${DATA_I}1.bim > comparisons/${COHORT_I}.bim
   ln -s ${DATA_I}1.bed comparisons/${COHORT_I}.bed
  fi

  if test -f comparisons/${COHORT_J}.fam; then
   echo ${COHORT_J} DATA ALREADY IN PLACE
  else
   awk -vc=${COHORT_J} '{print c"_"$1,c"_"$2,$3,$4,$5,$6}' ${DATA_J}1.fam > comparisons/${COHORT_J}.fam
   awk '{print $1,$1"_"$4,$3,$4,$5,$6}' ${DATA_J}1.bim > comparisons/${COHORT_J}.bim
   ln -s ${DATA_J}1.bed comparisons/${COHORT_J}.bed
  fi

  ## Merge data
  plink \
  --allow-no-sex \
  --bfile comparisons/${COHORT_I} \
  --bmerge comparisons/${COHORT_J} \
  --make-bed \
  --out comparisons/${COHORT_I}_${COHORT_J}

  ## If there are strand issues between pairs of studies we'll pick up here and think again...
  ## See https://www.cog-genomics.org/plink/1.9/data#merge3
  if test -f comparisons/${COHORT_I}_${COHORT_J}-merge.missnp; then
   echo NEED TO REVIEW STRAND FOR ${COHORT_I}_${COHORT_J}-merge.missnp
   touch comparisons/REVIEW_STRANDS_FOR_${COHORT_I}_${COHORT_J}
   break
  fi

  ## Find no LRLD and no PS regions
  plink \
  --allow-no-sex \
  --bfile comparisons/${COHORT_I}_${COHORT_J} \
  --exclude range /users/k1191584/psoriasis/repository/LRLD_and_PS_loci_20161215.regions \
  --maf 0.1 \
  --make-bed \
  --out comparisons/${COHORT_I}_${COHORT_J}_noLRLD_noHits

  ## ... and identify pairwise genotype LD correlations...
  plink \
  --allow-no-sex \
  --bfile comparisons/${COHORT_I}_${COHORT_J}_noLRLD_noHits \
  --indep-pairwise 1500 150 0.2 \
  --out comparisons/${COHORT_I}_${COHORT_J}_noLRLD_noHits

  ## ... and remove SNPs that are in LD...
  plink \
  --allow-no-sex \
  --bfile comparisons/${COHORT_I}_${COHORT_J}_noLRLD_noHits \
  --exclude comparisons/${COHORT_I}_${COHORT_J}_noLRLD_noHits.prune.out \
  --extract comparisons/${COHORT_I}_${COHORT_J}_noLRLD_noHits.prune.in \
  --make-bed \
  --out comparisons/${COHORT_I}_${COHORT_J}_noLRLD_noHits_pruned

  ## Run relatedness check
  king -b comparisons/${COHORT_I}_${COHORT_J}_noLRLD_noHits_pruned.bed --kinship --related --degree 2 --prefix comparisons/${COHORT_I}_${COHORT_J}_kingmerge

  ## Extract lines which are cross-study (i.e. not both samples in same study)
  grep ${COHORT_I} comparisons/${COHORT_I}_${COHORT_J}_kingmerge.kin0 | grep ${COHORT_J} \
   > comparisons/${COHORT_I}_${COHORT_J}_kingmerge.kin1

  touch comparisons/DONE_${COHORT_I}_${COHORT_J}
 done
done
