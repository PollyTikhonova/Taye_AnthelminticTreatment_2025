##################################################################
# Step 1. sample_table: create a .csv table that contains a list of all samples in the dataset.
#    Creates a sample table that has two columns in case of single-end data: sample, R1;
#    and three columns in case of paired-end data: sample, R1, R2.
#    "sample" column contains the names of the samples in the dataset, 
#    "R1" and "R2" columns contain paths to the correponding sample files, if existent.
# Output: .csv table
##################################################################

# load input data
args <- commandArgs(TRUE)

# read and save the input
reads_path <- args[[1]]        # a path to the reads folder
sample_table_path <- args[[2]] # a path to the output

# load libraries
library(utils)

# list all the file names in the reads directory
sample_names <- list.files(reads_path, pattern='', full.names = FALSE) 
# creates paths to the files for boths strands
fastq_r1 <- paste(reads_path, sample_names, paste0(sample_names, '_R1.fastq.gz'), sep='/')
fastq_r2 <- paste(reads_path, sample_names, paste0(sample_names, '_R2.fastq.gz'), sep='/')
# if all files for the second strand exist, create a table with three columns: sample, R1, R2;
# otherwise, just two: sample, R1.
# sample column contains sample names, R1/R2 columns contain paths to the correponding files
if (all(file.exists(fastq_r2) == TRUE)) {
    sample_table <- data.frame(sample = sample_names, R1 = fastq_r1, R2 = fastq_r2)
} else {
    sample_table <- data.frame(sample = sample_names, R1 = fastq_r1)
}
# saves the table to the tab-separated format
write.table(sample_table, 
            file = sample_table_path,
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE)