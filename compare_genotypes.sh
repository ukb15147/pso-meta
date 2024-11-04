## Pseudocode describing the process of identifying and removing participants that are genetically related
## across different psoriasis GWAS studies
##
## The process is based on sharing a minimal number of genetic markers from contributing studies. Sharing of
## genotype data was not permitted for the three biobank studies included (Estonian Biobank, HUNT, UKB), so
## these were not directly compared against one another. They were each compared to the 15 other included
## studies by local analysts without sharing genotype data.




## Step 1 - Directly typed variants outside of known psoriasis regions and shared between contributing
##          studies were identified by inspection of variant lists (e.g. bim files).
##          Due to limited total overlap, this was done using smaller shared lists of variants.
##          For UKB: commonmarkers_new9cohorts.noknownloci (6863 variants)
##                   ukb_binary_v2.chrpos.commonmarkers.allmichigan.cohorts.noknownloci (4729 variants)

PLINK=ukb_step4_interim  ## Post-QC bed/bim/fam filestem
COMMON1=commonmarkers_new9cohorts.noknownloci
COMMON2=ukb_binary_v2.chrpos.commonmarkers.allmichigan.cohorts.noknownloci

## Identify marker names for SNPs in list 1

PLINK_BASE=$(basename ${PLINK})

awk 'F==1{b=$1":"$4;a[b]=b}F==2{if(!($1 in a))print $0}' \
 F=1 ${PLINK}.bim F=2 ${COMMON1} > ${PLINK_BASE}_common1.missing

awk 'F==1{a[$1]=$1}F==2{if($1":"$4 in a)print $2}' \
 F=1 ${COMMON1} F=2 ${PLINK}.bim > tmp.list1

plink \
--allow-no-sex \
--bfile ${PLINK} \
--extract tmp.list1 \
--make-bed \
--out ${PLINK_BASE}_common1

## Identify marker names for SNPs in list 2

awk 'F==1{a[$1":"$4]=$1":"$4}F==2{if(!($1 in a))print $0}' \
 F=1 ${PLINK}.bim F=2 ${COMMON2} > ${PLINK_BASE}_common2.missing

awk 'F==1{a[$1]=$1}F==2{if($1":"$4 in a)print $2}' \
 F=1 ${COMMON2} F=2 ${PLINK}.bim > tmp.list2

plink \
--allow-no-sex \
--bfile ${PLINK} \
--extract tmp.list2 \
--make-bed \
--out ${PLINK_BASE}_common2

rm tmp.list1 tmp.list2




## Step 2 - For each pair of studies, use PLINK to merge variant subsets, then KING to identify relationships
##          and extract cross-dataset relationships via grep
##
##  MAIN STEPS: See example script s05_compare_9new.sh

## Collect together all inter-dataset relationships:
cat comparisons/*kin1 > comparisons/all_pairs_9new.kin1




## Step 3 - Analyse all pairs of inter-dataset relationships and generate exclusion lists according to a
##          hierarchy of preferences: exclude to remove duplicated participants before comparing first and
##          second degree relatives; exclude from a biobank-type dataset ahead of an ascertained case-control
##          dataset due to direct phenotyping; exclude participants with multiple relationships ahead of those
##          with single relationships; otherwise prioritise datasets using higher-coverage arrays
##
##  MAIN STEPS: See R script s08_extract_exclusions_final.R




## Step 4 - Identify participants to exclude from each individual GWAS (for example, UK Biobank GWAS
##          indicated by acronym "UKB" as variable $i in the loop

for i in $(cat FINAL_final_exclusion.cohorts); do
  echo $i
  awk -v i=${i} '$2==i' FINAL_final_exclusion_list_20191126.txt |\
    sed "s/^${i}_//" > FINAL_final_exclusion_list_${i}.list
done
