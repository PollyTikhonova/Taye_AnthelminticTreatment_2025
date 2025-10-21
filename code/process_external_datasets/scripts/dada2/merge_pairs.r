##################################################################
# Step 5. merge_pairs: merge dada2 results, recieved from two strands.
#    Merges dada ASVs from two strands into one object.
# Output: a dada dataframe in a .rds format.
##################################################################

# load input data
args <- commandArgs(TRUE)

# read and save the input
dadaR1 <- readRDS(args[[1]])  # read dada algorithm results for the forward strand
dadaR2 <- readRDS(args[[2]])  # read dada algorithm results for the reverse strand
derepR1 <- readRDS(args[[3]]) # read dereplicated reads for the forward strand
derepR2 <- readRDS(args[[4]]) # read dereplicated reads for the reverse strand
output_file <- args[[5]] # an output file path

# load libraries
library("dada2")
library("parallel")
# merge forward and reverse strands
mergers <- mergePairs(dadaR1, derepR1, dadaR2, derepR2, verbose=TRUE)
# save to the rds format
saveRDS(mergers, output_file)