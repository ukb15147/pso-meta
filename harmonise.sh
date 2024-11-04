## Steps to harmonise UK Biobank GWAS summary statistics against those of other meta-analysis studies
##
## - Run on HPC
## - Input files:
##    - PLINK GWAS summary stats
##       - PLINK2 output in logistic.hybrid format (column headings #CHROM, POS, ID, REF, ALT, A1, AX
##         A1_FREQ, A1_CASE_FREQ, A1_CTRL_FREQ, MACH_R2, FIRTH?, TESTOBS_CT, OR, LOG(OR)_SE, L95, U95
##         Z_STAT, P, ERRCODE)
##    - marker_panel_ref_refaltcohorts_exclUKB.txt
##       - list of reference alleles derived from other meta-analysis input datasets that were supplied
##         with ref allele (format 'chr:pos allele')
##    - marker_panel_ids_current.list
##       - list of variant IDs from previously processed meta-analysis input files. To be matched to so
##         that variants are coded with consistent IDs across studies
 

## 1. Get the ref and alt alleles together
## ---------------------------------------

(head -1 ../PS_ukb_r11_initial_case_control_MAF0_005_INFO7_chr1_v2.PHENO1.glm.logistic.hybrid;
 cat ../PS_ukb_r11_initial_case_control_MAF0_005_INFO7_chr*_v2.PHENO1.glm.logistic.hybrid |\
 grep -v ^#) | cut -f1,2,4,5 > summarise_alleles_results.txt

## Separate SNPs and non-SNPs

awk '{if($3 ~ /^[ACGT]$/ && $4 ~ /^[ACGT]$/){print $1,$2,$3,$4,$1":"$2"_"$3"_"$4}}' \
 summarise_alleles_results.txt > summarise_alleles_results_onlySNPs.txt 
awk '{if($3 !~ /^[ACGT]$/ || $4 !~ /^[ACGT]$/){print $1,$2,$3,$4,$1":"$2"_"$3"_"$4}}' \
 summarise_alleles_results.txt | sed 1d > summarise_alleles_results_nonSNPs.txt

rm summarise_alleles_results.txt


## Check for duplicated IDs
cut -f5 -d' ' summarise_alleles_results_onlySNPs.txt | sort | uniq -d > dup_ids_onlySNPs.txt
cut -f5 -d' ' summarise_alleles_results_nonSNPs.txt | sort | uniq -d > dup_ids_nonSNPs.txt

## There are none


## 2. Compare the UK Biobank alleles against a consensus list of those from other meta-analysis studies
## ----------------------------------------------------------------------------------------------------

## Check ref alleles same between SNP data from UKB and from combined panel

awk 'F==1{a[$1]=$2}F==2{if($1":"$2 in a && a[$1":"$2] != $3){print $1":"$2,a[$1":"$2],$3}}' \
 F=1 [path_on_hpc]/marker_panel_ref_refaltcohorts_exclUKB.txt \
 F=2 summarise_alleles_results_onlySNPs.txt \
 > ref_mismatches_snps_vsCombined.txt


## There are 88 differences to combined panel - all on chr 21 and 22. Checking a couple on dbSNP suggests
##  that UKB has the wrong reference allele. Drop ambiguous SNPs

awk 'BEGIN{a["A_T"]="A_T"; a["C_G"]="C_G"; a["G_C"]="G_C"; a["T_A"]="T_A"}{
 if($2"_"$3 in a){print $1}}' ref_mismatches_snps_vsCombined.txt \
 > ref_mismatches_exclude_v2.list

awk 'BEGIN{a["A_T"]="A_T"; a["C_G"]="C_G"; a["G_C"]="G_C"; a["T_A"]="T_A"}{
 if(!($2"_"$3 in a)){print $1}}' ref_mismatches_snps_vsCombined.txt \
 > ref_mismatches_swap_v2.list


## 3. Make merged UKB sumstats with chr:pos_allele_allele ID to match other meta-analysis input files
##---------------------------------------------------------------------------------------------------

echo "#CHR POS SNP A1 AX OR LOGOR_SE P" \
 > PS_ukb_r11_initial_case_control_MAF0_005_INFO7.logistic.hybrid_merged

cat ../PS_ukb_r11_initial_case_control_MAF0_005_INFO7_chr*_v2.PHENO1.glm.logistic.hybrid |\
 grep -v ^# |\
 awk -vc1="ref_mismatches_exclude_v2.list" -vc2="ref_mismatches_swap_v2.list" 'BEGIN{
  while(getline<c1){a[$1]=$1}
  while(getline<c2){b[$1]=$1}
 }{
  if($1>20) {
   if(!($1":"$2 in a)) {
    if(($1":"$2 in b) && ($4 ~ /[ACGT]/)) {
     ## The ref allele disagrees with BSTOP so ID needs to be fixed (swap ref/alt in ID)
     print $1,$2,$1":"$2"_"$5"_"$4,$6,$7,$15,$16,$20
    } else {
     print $1,$2,$1":"$2"_"$4"_"$5,$6,$7,$15,$16,$20
    }
   }
  } else { print $1,$2,$1":"$2"_"$4"_"$5,$6,$7,$15,$16,$20 }
 }' >> PS_ukb_r11_initial_case_control_MAF0_005_INFO7.logistic.hybrid_merged


## Make a comparison file
 
cat ../PS_ukb_r11_initial_case_control_MAF0_005_INFO7_chr*_v2.PHENO1.glm.logistic.hybrid |\
 grep -v ^# |\
 awk '{print $1,$2,$1":"$2"_"$4"_"$5,$6,$7,$15,$16,$20}' | sort > temp_comparison_file

sort PS_ukb_r11_initial_case_control_MAF0_005_INFO7.logistic.hybrid_merged > temp_main_file_sorted

diff temp_main_file_sorted temp_comparison_file > temp_diff_check_changes

## Fixes look fine. After checking can delete temp files
rm temp*


## 4. Reformat columns into final meta-analysis input format (incl final ID fixes)
##--------------------------------------------------------------------------------

MARKERS=[path_on_hpc]/marker_panel_ids_current.list

head -1 PS_ukb_r11_initial_case_control_MAF0_005_INFO7.logistic.hybrid_merged |\
 awk '{print $1,$2,$3,$4,$5,"INFO",$6,$7,$8,"PANEL"}' \
 > UKB_r11_final.summary2

sed 1d PS_ukb_r11_initial_case_control_MAF0_005_INFO7.logistic.hybrid_merged |\
 awk -F_ '{print $1,$2,$3}' |\
 awk -vc=${MARKERS} 'BEGIN{while(getline<c){a[$1]=$1}; i1=0; i2=0 }{
  if($3"_"$5"_"$4 in a) {
   i2++; print $1,$2,$3"_"$5"_"$4,$6,$7,".",$8,$9,$10,"."
  } else {
   i1++; print $1,$2,$3"_"$4"_"$5,$6,$7,".",$8,$9,$10,"."
  }} END {
     print "Num ref first",i1
     print "Num alt first",i2 }' \
 >> UKB_r11_final.summary2

grep ^Num UKB_r11_final.summary2 > countrefalt
grep -v ^Num UKB_r11_final.summary2 > UKB_r11_final.summary
rm UKB_r11_final.summary2
