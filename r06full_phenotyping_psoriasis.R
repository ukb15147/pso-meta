## NOTE THIS SHOULD BE A COPY OF r06 WITH THE LINE:
##   working.df <- df[1:1000, ]
## REPLACED BY THE FULL DATASET
##   working.df <- df
## AND THE OUTPUT FILENAMES SLIGHTLY ALTERED BY FINDING AND REPLACING r06_ WITH r06full_

## Script to pick out the different definitions of psoriasis

library(rhelpers)
library(reshape2)


## Load data and start with the first 1000 rows for efficiency

load("intermediate/initial_data.Rda")

key2 <- key
key2$suff <- paste0("_f", gsub("\\.", "_", gsub("-", "_", key2$field.html)))
key2$strip <- sapply(1:nrow(key2), function(i) gsub(key2[i, "suff"], "", key2[i, "col.name"]))

treatment.codes <- read.default("data_local/final_meds_list_pso.txt", as.char = TRUE)
treatment.codes <- subset(treatment.codes, INCLUDE_FINAL == "1")$code

hes.codes.excl.PsA <- c("L400", "L401", "L402", "L403", "L404", "L408", "L409")
hes.codes.incl.PsA <- c("L400", "L401", "L402", "L403", "L404", "L405", "L408", "L409")

primary.care.eids.auto <- read.default("data_local/psorisis_automapped_L40.eids", header = FALSE, as.char = TRUE)

working.df <- df

##----------- Later can set working.df to be the full dataset and rerun everything below -----------##



## Choose variables of interest

var_choose <- read.csv("excel_workbooks/uniq_vars_choose_20191218.csv",
                       stringsAsFactors = FALSE)
vars_needed <- subset(var_choose, diagnoses_all == 1)$Variable
vars <- subset(key2, strip %in% vars_needed)$col.name

working.df <- working.df[, vars]
rownames(working.df) <- NULL


##---------------------------##
## 1. Self-report definition ##
##---------------------------##

working.selfrep <- working.df[, grepp("noncancer_illness_code_selfreported", names(working.df))]
ps01.selfrep <- apply(working.selfrep, 1, function(x) {
  return(sum(x == "1453", na.rm = TRUE) > 0)
})
ps01.selfrep.ids <- working.df[ps01.selfrep, "eid"]

write.default(ps01.selfrep.ids, "output/r06full_ps01_selfreport.list", need.col.names = FALSE)
rm(ps01.selfrep.ids)


##-----------------------------##
## 2. Self-report visit 1 only ##
##-----------------------------##

working.selfrep.visit1 <- working.df[, grepp("noncancer_illness_code_selfreported_f20002_0", names(working.df))]
ps02.selfrep.visit1 <- apply(working.selfrep.visit1, 1, function(x) {
  return(sum(x == "1453", na.rm = TRUE) > 0)
})
ps02.selfrep.visit1.ids <- working.df[ps02.selfrep.visit1, "eid"]

write.default(ps02.selfrep.visit1.ids, "output/r06full_ps02_selfreport_visit1.list", need.col.names = FALSE)
rm(working.selfrep.visit1, ps02.selfrep.visit1.ids)


##------------------------------##
## 3. Self-report including PsA ##
##------------------------------##

ps03.selfrep.withPsA <- apply(working.selfrep, 1, function(x) {
  return(sum(x %in% c("1453", "1477"), na.rm = TRUE) > 0)
})
ps03.selfrep.withPsA.ids <- working.df[ps03.selfrep.withPsA, "eid"]

write.default(ps03.selfrep.withPsA.ids, "output/r06full_ps03_selfreport_withPsA.list", need.col.names = FALSE)
rm(working.selfrep, ps03.selfrep.withPsA.ids)


##------------------------##
## 4. Treatment code only ##
##------------------------##

## Represents having reported a medication that can be prescribed to treat psoriasis
working.srtreat <- working.df[, grepp("treatmentmedication_code", names(working.df))]
ps04.srtreat <- apply(working.srtreat, 1, function(x) {
  return(sum(x %in% treatment.codes, na.rm = TRUE) > 0)
})
ps04.srtreat.ids <- working.df[ps04.srtreat, "eid"]

## Most common medications reported
srtreats <- unlist(working.srtreat[ps04.srtreat, ])
srtreats <- srtreats[srtreats %in% treatment.codes]
srtreats <- sort(table(srtreats), decreasing = TRUE)
write.default(srtreats, "output/r06full_srtreat_counts.txt")
rm(srtreats)

write.default(ps04.srtreat.ids, "output/r06full_ps04_srtreat.list", need.col.names = FALSE)
rm(ps04.srtreat.ids)


##------------------------------------##
## 5. Self-report plus treatment code ##
##------------------------------------##

ps05.selfrep.srtreat <- ps01.selfrep & ps04.srtreat
ps05.selfrep.srtreat.ids <- working.df[ps05.selfrep.srtreat, "eid"]

## Most common medications reported in psoriasis patients
srtreats.ps <- unlist(working.srtreat[ps05.selfrep.srtreat, ])
srtreats.ps <- srtreats.ps[srtreats.ps %in% treatment.codes]
srtreats.ps <- sort(table(srtreats.ps), decreasing = TRUE)
write.default(srtreats.ps, "output/r06full_srtreat_counts_in_selfrep.txt")
rm(srtreats.ps)

write.default(ps05.selfrep.srtreat.ids, "output/r06full_ps05_selfrep_srtreat.list", need.col.names = FALSE)
rm(working.srtreat, ps05.selfrep.srtreat.ids)


##-----------------------##
## 6. HES main psoriasis ##
##-----------------------##

working.hes.main <- working.df[, grepp("diagnoses_main_icd10", names(working.df))]
working.hes.secondary <- working.df[, grepp("diagnoses_secondary_icd10", names(working.df))]
working.hes.speciality.recd <- working.df[, grepp("main_speciality_of_consultant_recoded", names(working.df))]

ps06.hes.main <- apply(working.hes.main, 1, function(x) {
  return(sum(x %in% hes.codes.excl.PsA, na.rm = TRUE) > 0)
})
ps06.hes.main.ids <- working.df[ps06.hes.main, "eid"]

write.default(ps06.hes.main.ids, "output/r06full_ps06_HESmain_noPsA.list", need.col.names = FALSE)
rm(ps06.hes.main.ids)


##-----------------------------------------------##
## 7. HES main, including arthropathic psoriasis ##
##-----------------------------------------------##

ps07.hes.main.withPsA <- apply(working.hes.main, 1, function(x) {
  return(sum(x %in% hes.codes.incl.PsA, na.rm = TRUE) > 0)
})
ps07.hes.main.withPsA.ids <- working.df[ps07.hes.main.withPsA, "eid"]

write.default(ps07.hes.main.withPsA.ids, "output/r06full_ps07_HESmain_withPsA.list", need.col.names = FALSE)
rm(working.hes.main, ps07.hes.main.withPsA.ids)


##----------------------------##
## 8. HES secondary psoriasis ##
##----------------------------##

ps08.hes.secondary <- apply(working.hes.secondary, 1, function(x) {
  return(sum(x %in% hes.codes.excl.PsA, na.rm = TRUE) > 0)
})
ps08.hes.secondary.ids <- working.df[ps08.hes.secondary, "eid"]

write.default(ps08.hes.secondary.ids, "output/r06full_ps08_HESsecondary_noPsA.list", need.col.names = FALSE)
rm(ps08.hes.secondary.ids)


##----------------------------------------------------##
## 9. HES secondary, including arthropathic psoriasis ##
##----------------------------------------------------##

ps09.hes.secondary.withPsA <- apply(working.hes.secondary, 1, function(x) {
  return(sum(x %in% hes.codes.incl.PsA, na.rm = TRUE) > 0)
})
ps09.hes.secondary.withPsA.ids <- working.df[ps09.hes.secondary.withPsA, "eid"]

write.default(ps09.hes.secondary.withPsA.ids, "output/r06full_ps09_HESsecondary_withPsA.list", need.col.names = FALSE)
rm(working.hes.secondary, ps09.hes.secondary.withPsA.ids)


##---------------------------##
## 10. HES main or secondary ##
##---------------------------##

ps10.hes.main.secondary <- ps06.hes.main | ps08.hes.secondary
ps10.hes.main.secondary.ids <- working.df[ps10.hes.main.secondary, "eid"]

write.default(ps10.hes.main.secondary.ids, "output/r06full_ps10_HESmain_or_secondary_noPsA.list", need.col.names = FALSE)
rm(ps10.hes.main.secondary.ids)


##-------------------------------------------------------------##
## 11. HES main or secondary, including arthropathic psoriasis ##
##-------------------------------------------------------------##

ps11.hes.main.secondary.withPsA <- ps07.hes.main.withPsA | ps09.hes.secondary.withPsA
ps11.hes.main.secondary.withPsA.ids <- working.df[ps11.hes.main.secondary.withPsA, "eid"]

write.default(ps11.hes.main.secondary.withPsA.ids, "output/r06full_ps11_HESmain_or_secondary_withPsA.list", need.col.names = FALSE)
rm(ps11.hes.main.secondary.withPsA.ids)


##--------------------------------------##
## 12. (Any) diagnosis by dermatologist ##
##--------------------------------------##

## Note this diagnosis may not necessarily refer to the diagnosis of psoriasis - but currently nothing better

ps12.hes.any.dermatologist <- apply(working.hes.speciality.recd, 1, function(x) {
  return(sum(x %in% c("1230"), na.rm = TRUE) > 0)
})
ps12.hes.any.dermatologist.ids <- working.df[ps12.hes.any.dermatologist, "eid"]

write.default(ps12.hes.any.dermatologist.ids, "output/r06full_ps12_HESanydermatologist.list", need.col.names = FALSE)
rm(working.hes.speciality.recd, ps12.hes.any.dermatologist.ids)


##------------------------------------------------------------------------##
## 13. HES main diagnosis of psoriasis and any diagnosis by dermatologist ##
##------------------------------------------------------------------------##

ps13.hes.main.any.dermatologist <- ps06.hes.main & ps12.hes.any.dermatologist
ps13.hes.main.any.dermatologist.ids <- working.df[ps13.hes.main.any.dermatologist, "eid"]

write.default(ps13.hes.main.any.dermatologist.ids, "output/r06full_ps13_HESmain_anydermatologist_noPsA.list", need.col.names = FALSE)
rm(ps13.hes.main.any.dermatologist.ids)


##---------------------------------------------------------------------------------------##
## 14. HES main diagnosis of psoriasis and any diagnosis by dermatologist, including PsA ##
##---------------------------------------------------------------------------------------##

ps14.hes.main.any.dermatologist.withPsA <- ps07.hes.main.withPsA & ps12.hes.any.dermatologist
ps14.hes.main.any.dermatologist.withPsA.ids <- working.df[ps14.hes.main.any.dermatologist.withPsA, "eid"]

write.default(ps14.hes.main.any.dermatologist.withPsA.ids, "output/r06full_ps14_HESmain_anydermatologist_withPsA.list", need.col.names = FALSE)
rm(ps14.hes.main.any.dermatologist.withPsA.ids)


##-------------------------------------------------------------------------------------##
## 15. HES main or secondary diagnosis of psoriasis and any diagnosis by dermatologist ##
##-------------------------------------------------------------------------------------##

ps15.hes.main.secondary.any.dermatologist <- ps10.hes.main.secondary & ps12.hes.any.dermatologist
ps15.hes.main.secondary.any.dermatologist.ids <- working.df[ps15.hes.main.secondary.any.dermatologist, "eid"]

write.default(ps15.hes.main.secondary.any.dermatologist.ids, "output/r06full_ps15_HESmain_secondary_anydermatologist_noPsA.list", need.col.names = FALSE)
rm(ps15.hes.main.secondary.any.dermatologist.ids)


##----------------------------------------------------------------------------------------------------##
## 16. HES main or secondary diagnosis of psoriasis and any diagnosis by dermatologist, including PsA ##
##----------------------------------------------------------------------------------------------------##

ps16.hes.main.secondary.any.dermatologist.withPsA <- ps11.hes.main.secondary.withPsA & ps12.hes.any.dermatologist
ps16.hes.main.secondary.any.dermatologist.withPsA.ids <- working.df[ps16.hes.main.secondary.any.dermatologist.withPsA, "eid"]

write.default(ps16.hes.main.secondary.any.dermatologist.withPsA.ids, "output/r06full_ps16_HESmain_secondary_anydermatologist_withPsA.list", need.col.names = FALSE)
rm(ps16.hes.main.secondary.any.dermatologist.withPsA.ids)


##------------------------------------------##
## 17. Any self report plus any HES, no PsA ##
##------------------------------------------##

ps17.selfrep.and.hes <- ps01.selfrep & ps10.hes.main.secondary
ps17.selfrep.and.hes.ids <- working.df[ps17.selfrep.and.hes, "eid"]

write.default(ps17.selfrep.and.hes.ids, "output/r06full_ps17_selfreport_AND_HESany.list", need.col.names = FALSE)
rm(ps17.selfrep.and.hes.ids)


##--------------------------------------------##
## 18. Any self report plus any HES, with PsA ##
##--------------------------------------------##

ps18.selfrep.and.hes.withPsA <- ps03.selfrep.withPsA & ps11.hes.main.secondary.withPsA
ps18.selfrep.and.hes.withPsA.ids <- working.df[ps18.selfrep.and.hes.withPsA, "eid"]

write.default(ps18.selfrep.and.hes.withPsA.ids, "output/r06full_ps18_selfreport_AND_HESany_withPsA.list", need.col.names = FALSE)
rm(ps18.selfrep.and.hes.withPsA.ids)


##----------------------------------------##
## 19. Any self report OR any HES, no PsA ##
##----------------------------------------##

ps19.selfrep.or.hes <- ps01.selfrep | ps10.hes.main.secondary
ps19.selfrep.or.hes.ids <- working.df[ps19.selfrep.or.hes, "eid"]

write.default(ps19.selfrep.or.hes.ids, "output/r06full_ps19_selfreport_OR_HESany.list", need.col.names = FALSE)
rm(ps19.selfrep.or.hes.ids)


##------------------------------------------##
## 20. Any self report OR any HES, with PsA ##
##------------------------------------------##

ps20.selfrep.or.hes.withPsA <- ps03.selfrep.withPsA | ps11.hes.main.secondary.withPsA
ps20.selfrep.or.hes.withPsA.ids <- working.df[ps20.selfrep.or.hes.withPsA, "eid"]

write.default(ps20.selfrep.or.hes.withPsA.ids, "output/r06full_ps20_selfreport_OR_HESany_withPsA.list", need.col.names = FALSE)
rm(ps20.selfrep.or.hes.withPsA.ids)


##-----------------------------------------------------------------##
## 21. At least two of self report, treatment, HES - excluding PsA ##
##-----------------------------------------------------------------##

ps21.two.of.selfrep.srtreat.hes <- (ps01.selfrep & ps04.srtreat) | (ps01.selfrep & ps10.hes.main.secondary) | (ps04.srtreat & ps10.hes.main.secondary)
ps21.two.of.selfrep.srtreat.hes.ids <- working.df[ps21.two.of.selfrep.srtreat.hes, "eid"]

write.default(ps21.two.of.selfrep.srtreat.hes.ids, "output/r06full_ps21_two_of_selfreport_srtreat_HES_noPsA.list", need.col.names = FALSE)
rm(ps21.two.of.selfrep.srtreat.hes.ids)


##-----------------------------------------------------------------##
## 22. At least two of self report, treatment, HES - including PsA ##
##-----------------------------------------------------------------##

ps22.two.of.selfrep.srtreat.hes.withPsA <- (ps03.selfrep.withPsA & ps04.srtreat) | (ps03.selfrep.withPsA & ps11.hes.main.secondary.withPsA) | (ps04.srtreat & ps11.hes.main.secondary.withPsA)
ps22.two.of.selfrep.srtreat.hes.withPsA.ids <- working.df[ps22.two.of.selfrep.srtreat.hes.withPsA, "eid"]

write.default(ps22.two.of.selfrep.srtreat.hes.withPsA.ids, "output/r06full_ps22_two_of_selfreport_srtreat_HES_withPsA.list", need.col.names = FALSE)
rm(ps22.two.of.selfrep.srtreat.hes.withPsA.ids)


##--------------------------------------------------##
## 23. Primary care - "auto" mapped - including PsA ##
##--------------------------------------------------##

## "auto" mapped means it uses the default mapping of read codes to ICD10 codes
## See Rosalind directory derm_ukb/primary_care_??

ps23.primary.care.auto <- working.df$eid %in% primary.care.eids.auto[, 1]
ps23.primary.care.auto.ids <- working.df[ps23.primary.care.auto, "eid"]

write.default(ps23.primary.care.auto.ids, "output/r06full_ps23_primary_care_auto.list", need.col.names = FALSE)
rm(ps23.primary.care.auto.ids)


##-------------------------------------------------------------------------------##
## 24. At least two of self report, treatment, HES, primary care - including PsA ##
##-------------------------------------------------------------------------------##

ps24.two.of.selfrep.srtreat.hes.primary.withPsA <-
  (ps03.selfrep.withPsA & ps04.srtreat) |
  (ps03.selfrep.withPsA & ps11.hes.main.secondary.withPsA) |
  (ps03.selfrep.withPsA & ps23.primary.care.auto) |
  (ps04.srtreat & ps11.hes.main.secondary.withPsA) |
  (ps04.srtreat & ps23.primary.care.auto) |
  (ps11.hes.main.secondary.withPsA & ps23.primary.care.auto)
ps24.two.of.selfrep.srtreat.hes.primary.withPsA.ids <- working.df[ps24.two.of.selfrep.srtreat.hes.primary.withPsA, "eid"]

write.default(ps24.two.of.selfrep.srtreat.hes.primary.withPsA.ids, "output/r06full_ps24_two_of_selfreport_srtreat_HES_primary_withPsA.list", need.col.names = FALSE)
rm(ps24.two.of.selfrep.srtreat.hes.primary.withPsA.ids)

