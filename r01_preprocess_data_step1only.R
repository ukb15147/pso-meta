## Initial script to convert downloaded data to a useable dataframe using
## the ukbtools package (doi:10.32614/CRAN.package.ukbtools)

library(ukbtools)

## Make a data-frame
df <- ukb_df("ukb8308", path = "../rawdata")

## Make a key
key <- ukb_df_field("ukb8308", path = "../rawdata")

## Takes several minutes: save it! (Note saving takes several minutes too!)
save(df, key, file = "intermediate/initial_data.Rda")
