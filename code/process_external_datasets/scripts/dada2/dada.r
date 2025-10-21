##################################################################
# Step 4. dada: run the core dada2 algorithm of ASV identification.
#     Utilizes a error model to identify ASVs on dereplicated data.
#     Runs forward and reverse strands separately.
#     Creates one object per strand for all samples in the dataset.
# Output: a dada-class object in .rds format
##################################################################


# load input data
args <- commandArgs(TRUE)

# read and save the input
derep <- readRDS(args[[1]])  # read the file with dereplicated reads
err <- readRDS(args[[2]])    # read the error model
out <- c(args[[3]]) # a path for the output 
params <- args[[4]] # read a pool parameter for the dada pipeline

# Unlsee the pool parameter is pseudo, transform to the bool format
if (params != 'pseudo') {
    params <- as.logical(params)
}

# load libraries
library("dada2")
# run a dada algorithm
pool <- dada(derep, err=err, multithread=TRUE, pool=params, verbose=TRUE) 
# save the result to the rds format
saveRDS(pool, out)