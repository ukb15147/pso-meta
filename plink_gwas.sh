## Run using PLINK v2.00a3LM 64-bit Intel (11 May 2020)
## Run separately for each chromosome 1-22. Example code for chr3

plink2 \
--pgen [path_on_hpc]/ukb15147_imp_whiteBrit_unrel_MAF0_005_INFO7_chr3.pgen \
--psam [path_on_hpc]/ukb15147_imp_whiteBrit_unrel_MAF0_005_INFO7.psam \
--pvar [path_on_hpc]/ukb15147_imp_whiteBrit_unrel_MAF0_005_INFO7_chr3.pvar \
--remove [path_on_hpc]/withdrawals_20200204.exclude meta_cross_dataset_relateds.exclude \
--pheno r11_initial_case_control.pheno \
--covar [path_on_hpc]/samples_20PCs_array.covar \
--glm hide-covar cols=+a1freq,+a1freqcc,+machr2,+ax \
--ci 0.95 \
--out PS_ukb_r11_initial_case_control_MAF0_005_INFO7_chr3_v2
