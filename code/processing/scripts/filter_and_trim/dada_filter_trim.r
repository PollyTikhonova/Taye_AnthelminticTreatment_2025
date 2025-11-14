# load libraries
library(yaml) 
library(dada2)
# load input data
args <- commandArgs(TRUE)

fnFs <- c(args[[2]])   # input raw reads for forward strand
fnRs <- c(args[[3]])   # input raw reads for reverse strand
filtFs <- c(args[[4]]) # output filtered reads for forward strand
filtRs <- c(args[[5]]) # output filtered reads for reverse strand
params <- yaml.load_file(args[[6]]) # load parameters for the filter_and_trim function

# perform filter_amd_trimming of the sample
res <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen = c(params$truncLenF, params$truncLenR), 
                     trimLeft = c(params$trimLeftF, params$trimLeftR),
                     maxN = params$maxN, 
                     maxEE = c(as.numeric(params$maxEEF), as.numeric(params$maxEER)), 
                     truncQ = params$truncQ, 
                     rm.phix = params$rm.phix, 
                     compress = TRUE, 
                     minLen = params$minLen, 
                     maxLen = as.numeric(params$maxLen))