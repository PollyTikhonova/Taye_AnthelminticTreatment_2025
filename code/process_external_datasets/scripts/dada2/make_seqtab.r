##################################################################
# Step 6. make_seqtab: make an ASV table
#    Creates an ASV table from dada object (merged in case of paired-end reads).
# Output: an ASV table of read counts in a .rds format.
##################################################################

# load input data
args <- commandArgs(TRUE)

# read and save the input
mergers <- readRDS(args[[1]]) # read the resulting dada table (merged in case of paired-end sequences)
out <- c(args[[2]]) # an output file path

# load libraries
library("dada2")
# create a sequence table
seqtab <- makeSequenceTable(mergers)
# save to the rds format
saveRDS(seqtab, out)
