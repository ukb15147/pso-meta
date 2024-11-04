library(rhelpers)

## NOW PRIORITISE MIKIEL BEHIND ERLANGEN BEHIND KIEL

## Will all need redoing with revised kin1 files for Michigan group
## Plus make sure we have final est biobank files

## King kinship cutoffs: >0.354, [0.177, 0.354], [0.0884, 0.177] and [0.0442, 0.0884]

##------##
## Data ##
##------##

## Load data on non-Mich cohorts

kin1.9new <- read.default("data/all_pairs_9new.kin1", header = FALSE)
kin1.9new <- kin1.9new[, c(1, 3, 8)]
names(kin1.9new) <- c("Sample1", "Sample2", "Kinship")
kin1.9new$Cohort1 <- sapply(kin1.9new$Sample1, function(x) strsplit(x, "_")[[1]][1])
kin1.9new$Cohort2 <- sapply(kin1.9new$Sample2, function(x) strsplit(x, "_")[[1]][1])

## Load data on Mich cohorts

kin1.8mich <- read.default("data/all_pairs_8mich.kin1", header = FALSE)
kin1.8mich <- kin1.8mich[, c(1, 3, 8)]
names(kin1.8mich) <- c("Sample1", "Sample2", "Kinship")
kin1.8mich$Cohort1 <- sapply(kin1.8mich$Sample1, function(x) strsplit(x, "_")[[1]][1])
kin1.8mich$Cohort2 <- sapply(kin1.8mich$Sample2, function(x) strsplit(x, "_")[[1]][1])

## What are all of the cohorts?
shareable.cohorts <- unique(c(kin1.9new$Cohort1, kin1.9new$Cohort2, kin1.8mich$Cohort1, kin1.8mich$Cohort2))

## Add some manually that are not represented
shareable.cohorts <- sort(c(shareable.cohorts, "EST", "MiGENIZON"))

## Load HUNT data
hunt.vs.9new <- read.default("[filepath]/hunt_genotyped_common1_full_BetweenCohortRelated.kin0")
hunt.vs.8mich <- read.default("[filepath]/hunt_genotyped_common2_full_BetweenCohortRelated.kin0")
hunt <- rbind(hunt.vs.9new, hunt.vs.8mich)
rm(hunt.vs.8mich, hunt.vs.9new)
hunt <- hunt[, c(1, 3, 8)]
names(hunt) <- c("Sample1", "Sample2", "Kinship")
hunt$Cohort1 <- sapply(hunt$Sample1, function(x) strsplit(x, "_")[[1]][1])
hunt$Cohort2 <- "HUNT"

## Load EST_biobank data
est.vs.9new <- read.default("[filepath]/GSA_CLEANED_AUTOSOMAL_common1_full.kin0")
est.vs.9new <- est.vs.9new[, c("ID1", "ID2", "Kinship")]

est.vs.8mich <- read.default("[filepath]/GSA_CLEANED_AUTOSOMAL_common2_full.kin0")
est.vs.8mich <- est.vs.8mich[, c("ID1", "ID2", "Kinship")]

est <- rbind(est.vs.9new, est.vs.8mich)
rm(est.vs.9new, est.vs.8mich)
names(est) <- c("Sample1", "Sample2", "Kinship")
est$Cohort1 <- sapply(est$Sample1, function(x) strsplit(x, "_")[[1]][1])
est$Cohort1 <- ifelse(est$Cohort1 %in% shareable.cohorts, est$Cohort1, "ESTBIOBANK")
est$Cohort2 <- sapply(est$Sample2, function(x) strsplit(x, "_")[[1]][1])
est$Cohort2 <- ifelse(est$Cohort2 %in% shareable.cohorts, est$Cohort2, "ESTBIOBANK")


## Think within cohort has been sorted
within.8mich <- read.default("data/within_8mich.kin0", header = FALSE)


##--------------------------------##
## Order of priority - now edited ##
##--------------------------------##

order.of.priority.final <- data.frame(position=1:19, row.names = c("MiPSAGWAS", "GLVR", "LIAO", "SCOT", "KIEL19", "ERLANG", "MiKIEL",
                                                               "MiGSAMICH2", "MiGSATORONTO2", "MiWTCCC2",
                                                               "BSTOP", "MiGENIZON", "MiCASP", "EST", "MiEXMCHIP",
                                                               "MANC", "UKB", "ESTBIOBANK", "HUNT"))
order.of.priority.final["LIAO", "position"] <- order.of.priority.final["GLVR", "position"]
order.of.priority.final["SCOT", "position"] <- order.of.priority.final["GLVR", "position"]
order.of.priority.final["MiGSATORONTO2", "position"] <- order.of.priority.final["MiGSAMICH2", "position"]
order.of.priority.final["ESTBIOBANK", "position"] <- order.of.priority.final["UKB", "position"]
order.of.priority.final["HUNT", "position"] <- order.of.priority.final["UKB", "position"]


##----------##
## Analysis ##
##----------##

est.within <- subset(est, Cohort1 == "ESTBIOBANK" & Cohort2 == "ESTBIOBANK")

est.outside <- subset(est, !(Cohort1 == "ESTBIOBANK" & Cohort2 == "ESTBIOBANK"))
est.outside <- subset(est.outside, Cohort1 == "ESTBIOBANK" | Cohort2 == "ESTBIOBANK")

mich8.outside <- subset(kin1.8mich, Cohort1 != Cohort2)
new9.outside <- subset(kin1.9new, Cohort1 != Cohort2)

all.outside <- rbind(rbind(new9.outside, mich8.outside), rbind(hunt, est.outside))

all.outside$CohortA <- ifelse(all.outside$Cohort1 < all.outside$Cohort2,
                              all.outside$Cohort1, all.outside$Cohort2)
all.outside$CohortB <- ifelse(all.outside$Cohort1 > all.outside$Cohort2,
                              all.outside$Cohort1, all.outside$Cohort2)

cohort.levels <- sort(unique(c(all.outside$CohortA, all.outside$CohortB)))
all.outside$CohortA <- factor(all.outside$CohortA, levels = cohort.levels)
all.outside$CohortB <- factor(all.outside$CohortB, levels = cohort.levels)

#all.outside$Cohorts <- paste0(all.outside$CohortA, "_", all.outside$CohortB)
all.cohort.pairs <- table(all.outside[, c("CohortA", "CohortB")])
write.default(all.cohort.pairs, "output/FINAL_all_pair_counts.txt", TRUE, TRUE)

all.outside$relatedness <- ifelse(all.outside$Kinship > 0.354, "duplicate",
                                  ifelse(all.outside$Kinship > 0.177, "first-degree",
                                         ifelse(all.outside$Kinship > 0.10612,
                                                "second-degree-main",
                                                "second-degree-bottom-quintile")))

## Duplicates
all.outside.dups <- subset(all.outside, relatedness == "duplicate")
all.cohort.pairs.dups <- table(all.outside.dups[, c("CohortA", "CohortB")])
write.default(all.cohort.pairs.dups, "output/FINAL_all_pair_dups_counts.txt", TRUE, TRUE)

## First-degree
all.outside.1st <- subset(all.outside, relatedness == "first-degree")
all.cohort.pairs.1st <- table(all.outside.1st[, c("CohortA", "CohortB")])
write.default(all.cohort.pairs.1st, "output/FINAL_all_pair_1st_counts.txt", TRUE, TRUE)

## Second-degree (main)
all.outside.2nd.main <- subset(all.outside, relatedness == "second-degree-main")
all.cohort.pairs.2nd.main <- table(all.outside.2nd.main[, c("CohortA", "CohortB")])
write.default(all.cohort.pairs.2nd.main, "output/FINAL_all_pair_2nd_main_counts.txt", TRUE, TRUE)

## Second-degree (bottom quintile)
all.outside.2nd.low <- subset(all.outside, relatedness == "second-degree-bottom-quintile")
all.cohort.pairs.2nd.low <- table(all.outside.2nd.low[, c("CohortA", "CohortB")])
write.default(all.cohort.pairs.2nd.low, "output/FINAL_all_pair_2nd_bottom-quintile_counts.txt", TRUE, TRUE)


## Estonian Biobank appears to have inflated relatedness - ignore anything in lowest 20% of second-degree relative range
all.outside.2nd.low.2 <- subset(all.outside.2nd.low, !(Cohort1 == "ESTBIOBANK" | Cohort2 == "ESTBIOBANK"))
all.outside.2nd <- rbind(all.outside.2nd.main, all.outside.2nd.low.2)


## Remove duplicates

set.seed(1209983)

all.outside.dups$priority1 <- order.of.priority.final[all.outside.dups$Cohort1, "position"]
all.outside.dups$priority2 <- order.of.priority.final[all.outside.dups$Cohort2, "position"]
all.outside.dups$random <- runif(nrow(all.outside.dups))
all.outside.dups$exclude.sample <- with(all.outside.dups, ifelse(priority1 < priority2, Sample2,
                                                                 ifelse(priority1 > priority2, Sample1,
                                                                        ifelse(random < 0.5, Sample1, Sample2))))
all.outside.dups$exclude.cohort <- with(all.outside.dups, ifelse(priority1 < priority2, Cohort2,
                                                                 ifelse(priority1 > priority2, Cohort1,
                                                                        ifelse(random < 0.5, Cohort1, Cohort2))))
duplicate.exclusions <- unique(all.outside.dups[, c("exclude.sample", "exclude.cohort")])
duplicate.exclusions$reason <- "duplicate"


## Update other cohorts to exclude relationships with excluded dups

all.outside.1st.2 <- subset(all.outside.1st, !(Sample1 %in% duplicate.exclusions$exclude.sample |
                                                 Sample2 %in% duplicate.exclusions$exclude.sample))
all.outside.2nd.2 <- subset(all.outside.2nd, !(Sample1 %in% duplicate.exclusions$exclude.sample |
                                                 Sample2 %in% duplicate.exclusions$exclude.sample))


## Get rid of Biobanks

all.outside.1st.plus <- rbind(all.outside.1st.2, all.outside.2nd.2)
all.outside.1st.plus$exclude.sample <- with(all.outside.1st.plus,
                                            ifelse(Cohort1 == "UKB" | Cohort1 == "ESTBIOBANK" | Cohort1 == "HUNT",
                                                   Sample1,
                                                   ifelse(Cohort2 == "UKB" | Cohort2 == "ESTBIOBANK" | Cohort2 == "HUNT",
                                                          Sample2, NA)))
all.outside.1st.plus$exclude.cohort <- with(all.outside.1st.plus,
                                            ifelse(Cohort1 == "UKB" | Cohort1 == "ESTBIOBANK" | Cohort1 == "HUNT",
                                                   Cohort1,
                                                   ifelse(Cohort2 == "UKB" | Cohort2 == "ESTBIOBANK" | Cohort2 == "HUNT",
                                                          Cohort2, NA)))
biobank.exclusions <- na.omit(unique(all.outside.1st.plus[, c("exclude.sample", "exclude.cohort")]))
biobank.exclusions$reason <- "biobank"


all.outside.1st.plus.remaining <- subset(all.outside.1st.plus, is.na(exclude.sample))


## Now handle the first- and second-degree relationships accounting for multiple relationships

how.many.relationships <- as.data.frame(sort(table(c(all.outside.1st.plus.remaining$Sample1,
                                                     all.outside.1st.plus.remaining$Sample2)), decreasing = TRUE),
                                        stringsAsFactors = FALSE)
row.names(how.many.relationships) <- how.many.relationships[, 1]
how.many.relationships <- how.many.relationships[, -1, drop = FALSE]

single.relationships <- rownames(subset(how.many.relationships, Freq == 1))
multiple.relationships <- rownames(subset(how.many.relationships, Freq > 1))

all.outside.1st.plus.singles <- subset(all.outside.1st.plus.remaining,
                                       Sample1 %in% single.relationships &
                                         Sample2 %in% single.relationships)
all.outside.1st.plus.multiples <- subset(all.outside.1st.plus.remaining,
                                         Sample1 %in% multiple.relationships |
                                           Sample2 %in% multiple.relationships)


## Multiples first - keep samples that only appear once but are related to a sample appearing more than once

all.outside.1st.plus.multiples$count1 <- ifelse(all.outside.1st.plus.multiples$Sample1 %in% rownames(how.many.relationships),
                                                how.many.relationships[all.outside.1st.plus.multiples$Sample1, 1], NA)
all.outside.1st.plus.multiples$count2 <- ifelse(all.outside.1st.plus.multiples$Sample2 %in% rownames(how.many.relationships),
                                                how.many.relationships[all.outside.1st.plus.multiples$Sample2, 1], NA)

all.outside.1st.plus.multiples$exclude.sample <- with(all.outside.1st.plus.multiples,
                                                      ifelse(count1 == 1, Sample2,
                                                             ifelse(count2 == 1, Sample1, NA)))
all.outside.1st.plus.multiples$exclude.cohort <- with(all.outside.1st.plus.multiples,
                                                      ifelse(count1 == 1, Cohort2,
                                                             ifelse(count2 == 1, Cohort1, NA)))
multiples.exclusions.1 <- na.omit(unique(all.outside.1st.plus.multiples[, c("exclude.sample", "exclude.cohort")]))
multiples.exclusions.1$reason <- "multiple.1"


all.outside.1st.plus.multiples.2 <- subset(all.outside.1st.plus.multiples, !(Sample1 %in% multiples.exclusions.1$exclude.sample |
                                                                               Sample2 %in% multiples.exclusions.1$exclude.sample))

## The remaining "multiple" samples are presumable trios or more complex family setups.
## Just exclude the lowest coverage from each pair and we will be left with the best of each family

all.outside.1st.plus.multiples.2$priority1 <- order.of.priority.final[all.outside.1st.plus.multiples.2$Cohort1, "position"]
all.outside.1st.plus.multiples.2$priority2 <- order.of.priority.final[all.outside.1st.plus.multiples.2$Cohort2, "position"]
all.outside.1st.plus.multiples.2$random <- runif(nrow(all.outside.1st.plus.multiples.2))
all.outside.1st.plus.multiples.2$exclude.sample <- with(all.outside.1st.plus.multiples.2, ifelse(priority1 < priority2, Sample2,
                                                                                                 ifelse(priority1 > priority2, Sample1,
                                                                                                        ifelse(random < 0.5, Sample1, Sample2))))
all.outside.1st.plus.multiples.2$exclude.cohort <- with(all.outside.1st.plus.multiples.2, ifelse(priority1 < priority2, Cohort2,
                                                                                                 ifelse(priority1 > priority2, Cohort1,
                                                                                                        ifelse(random < 0.5, Cohort1, Cohort2))))
multiples.exclusions.2 <- unique(all.outside.1st.plus.multiples.2[, c("exclude.sample", "exclude.cohort")])
multiples.exclusions.2$reason <- "multiple.2"


## Exclude least coverage of single relationship pairs

all.outside.1st.plus.singles$priority1 <- order.of.priority.final[all.outside.1st.plus.singles$Cohort1, "position"]
all.outside.1st.plus.singles$priority2 <- order.of.priority.final[all.outside.1st.plus.singles$Cohort2, "position"]
all.outside.1st.plus.singles$random <- runif(nrow(all.outside.1st.plus.singles))
all.outside.1st.plus.singles$exclude.sample <- with(all.outside.1st.plus.singles, ifelse(priority1 < priority2, Sample2,
                                                                                         ifelse(priority1 > priority2, Sample1,
                                                                                                ifelse(random < 0.5, Sample1, Sample2))))
all.outside.1st.plus.singles$exclude.cohort <- with(all.outside.1st.plus.singles, ifelse(priority1 < priority2, Cohort2,
                                                                                         ifelse(priority1 > priority2, Cohort1,
                                                                                                ifelse(random < 0.5, Cohort1, Cohort2))))
singles.exclusions <- unique(all.outside.1st.plus.singles[, c("exclude.sample", "exclude.cohort")])
singles.exclusions$reason <- "single"


## EXCLUSION LISTS
##
## 1. duplicate.exclusions (n = 5797)
## 2. biobank.exclusions (n = 2637)
## 3. multiples.exclusions.1 (n = 33)
## 4. multiples.exclusions.2 (n = 12)
## 5. singles.exclusions (n = 372)


##---------------------------------------##
## Check that we break all relationships ##
##---------------------------------------##

all.exclusions <- do.call(rbind, list(duplicate.exclusions,
                                      biobank.exclusions,
                                      multiples.exclusions.1,
                                      multiples.exclusions.2,
                                      singles.exclusions))
sum(duplicated(all.exclusions))
sum(duplicated(all.exclusions$exclude.sample))

all.outside$is.sample1.excluded <- all.outside$Sample1 %in% all.exclusions$exclude.sample
all.outside$is.sample2.excluded <- all.outside$Sample2 %in% all.exclusions$exclude.sample
all.outside$are.samples.excluded <- paste(all.outside$is.sample1.excluded, all.outside$is.sample2.excluded, sep = "_")

check.false.false <- subset(all.outside, are.samples.excluded == "FALSE_FALSE")
tblNA(check.false.false$relatedness)
tblNA(check.false.false$Cohort2)
## It is all Estonian Biobank cases

check.true.true <- subset(all.outside, are.samples.excluded == "TRUE_TRUE")
tblNA(check.true.true[, c("Cohort1", "Cohort2")])
# 314 are MiKIEL and KIEL19 - could also be caught by Erlangen
# 300 are ESTBIOBANK and MiEXMCHIP - both of which had lots of overlap elsewhere
#   so could each have been lost as two different relationships
# 83 are MiEXMCHIP and GLVR - etc
# etc
check.true.true.kiel.kiel <- subset(check.true.true, Cohort1 == "KIEL19" & Cohort2 == "MiKIEL")
# Spot checks suggest reasonable

write.default(all.exclusions, "output/FINAL_final_exclusion_list_20191126.txt")

cohort.numbers <- as.data.frame(table(all.exclusions$exclude.cohort), stringsAsFactors = FALSE)
names(cohort.numbers) <- c("Cohort", "Number_to_exclude")
write.default(cohort.numbers, "output/FINAL_final_numbers_to_excluce_20191126.txt")
