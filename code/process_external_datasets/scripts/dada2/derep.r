##################################################################
# Step 3. derep: dereplicate reads.
#    Dereplicates the reads for each sample in the dataset.
#    Creates one object per strand for all samples in the dataset.
# Output: a derep-class object in a .rds format.
##################################################################


# load input data
args <- commandArgs(TRUE)

# read and save the input
reads <- read.table(file = args[[1]], sep = '\t', header = TRUE) # read a sample table
strand <- args[[2]] # read the strand the algorithm will be applied to
output_path <- args[[3]] # an output file path

# load libraries
library("dada2")

# ensure non-factors format
files_in <- as.character(reads[[strand]])
samples <- as.character(reads$sample)

# dereplicate
derep <- derepFastq(files_in, verbose=TRUE)
# ensure the correct names
names(derep) <- samples
# save to the rds format
saveRDS(derep, output_path)