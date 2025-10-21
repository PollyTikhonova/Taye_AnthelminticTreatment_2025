##################################################################
# Step 2. learn_errors: run a dada2 error model.
#    Runs a error model for all samples. In case of paired end sequencing, 
#    each strand is processed separately.
#    Creates one object per strand for all samples in the dataset.
# Output: a named list of error model in a .rds format.
##################################################################

# load input data
args <- commandArgs(TRUE)

# read and save the input
reads <- read.table(file = args[[1]], sep = '\t', header = TRUE) # read a sample table
strand <- args[[2]]  # read the strand the algorithm will be applied to
err_out <- args[[3]] # an output file path

# load libraries
library("dada2")
# ensure non-factors format
err_in <- as.character(reads[[strand]])
# run a error model
err <- learnErrors(err_in, multithread=TRUE, randomize=FALSE, MAX_CONSIST=10)
# save to the rds format
saveRDS(err, err_out)