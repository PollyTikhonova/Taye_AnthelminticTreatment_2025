# load libraries
library(yaml) 
library(dada2)
# load input data
args <- commandArgs(TRUE)

if (length(args) == 6){
    fnFs <- c(args[[1]])   # input raw reads for forward strand
    fnRs <- c(args[[2]])   # input raw reads for reverse strand
    params <- yaml.load_file(args[[3]]) # load parameters for the filter_and_trim function
    filtFs <- c(args[[4]]) # output filtered reads for forward strand
    filtRs <- c(args[[5]]) # output filtered reads for reverse strand
    log <- c(args[[6]]) # log file with filtering stats

    # perform filter_amd_trimming of the sample
    res <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                        truncLen = c(params$truncLen_f, params$truncLen_r), 
                        trimLeft = c(params$trimLeft_l, params$trimLeft_r),
                        maxN = params$maxN[1], 
                        maxEE = c(as.numeric(params$maxEE_f), as.numeric(params$maxEE_r)), 
                        truncQ = params$truncQ[1], 
                        rm.phix = params$rm_phix[1], 
                        compress = TRUE, 
                        minLen = params$minLen[1], 
                        maxLen = as.numeric(params$maxLen))
} else {
    fnFs <- c(args[[1]])   # input raw reads for forward strand
    params <- yaml.load_file(args[[2]]) # load parameters for the filter_and_trim function
    filtFs <- c(args[[3]]) # output filtered reads for forward strand
    log <- c(args[[4]]) # log file with filtering stats

    # perform filter_amd_trimming of the sample
    res <- filterAndTrim(fnFs, filtFs,
                        truncLen = params$truncLen_f,
                        trimLeft = params$trimLeft_l,
                        maxN = params$maxN, 
                        maxEE = as.numeric(params$maxEE_f),
                        truncQ = params$truncQ, 
                        rm.phix = params$rm_phix, 
                        compress = TRUE, 
                        minLen = params$minLen, 
                        maxLen = as.numeric(params$maxLen))
}
capture.output(print(res), file = log)