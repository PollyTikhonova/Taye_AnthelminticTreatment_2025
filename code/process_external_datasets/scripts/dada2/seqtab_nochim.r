##################################################################
# Step 7. seqtab_nochim: remove bimeras
#    Removes bimeras.
# Output: an ASV table in a .rds format.
##################################################################

# load input data
args <- commandArgs(TRUE)

# read and save the input
seqtab <- readRDS(args[[1]]) # read a sequence table
out <- c(args[[2]]) # an output file path

# load libraries
library("dada2")
# remove bimeras
seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=24)
# save to the rds format
saveRDS(seqtab_nochim, out)
