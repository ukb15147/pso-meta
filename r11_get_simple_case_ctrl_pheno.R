## Create an initial list of cases and controls for GWAS meta-analysis
## Picking a sensible definition of PS as of 12 May 2020; may revise in future

library(rhelpers)

fam <- read.default("data_local/ukb_step4_interim.fam", as.char = TRUE, header = FALSE, sep = " ")

case.list <- read.default("output/r06full_ps24_two_of_selfreport_srtreat_HES_primary_withPsA.list", as.char = TRUE, header = FALSE)


potential.exclusions.1 <- read.default("output/r06full_ps03_selfreport_withPsA.list", as.char = TRUE, header = FALSE)
potential.exclusions.2 <- read.default("output/r06full_ps04_srtreat.list", as.char = TRUE, header = FALSE)
potential.exclusions.3 <- read.default("output/r06full_ps11_HESmain_or_secondary_withPsA.list", as.char = TRUE, header = FALSE)
potential.exclusions.4 <- read.default("output/r06full_ps23_primary_care_auto.list", as.char = TRUE, header = FALSE)

## Leave out the treatments as exclusion criteria because broader than just pso
potential.exclusions <- unique(c(potential.exclusions.1[, 1],
                                 potential.exclusions.3[, 1],
                                 potential.exclusions.4[, 1]))
exclusions <- potential.exclusions[!potential.exclusions %in% case.list[, 1]]

## Write out
pheno <- fam[, c(1, 2, 6)]
pheno[, 3] <- 1
pheno[pheno[, 1] %in% case.list[, 1], 3] <- 2
pheno[pheno[, 1] %in% exclusions, 3] <- -9

write.default(pheno, "output/r11_initial_case_control.pheno", FALSE, FALSE)
