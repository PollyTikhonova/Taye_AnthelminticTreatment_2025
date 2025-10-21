##################################################################
# Step 8. assign_taxa: assign taxonomy
#    Assigns taxonomy using SILVA and SILVA species databases.
# Output: a taxonomy table in a .rds format with columns: Kingdom, Phylum, Class, Order, Family, Genus, Species.
##################################################################

# load input data
args <- commandArgs(TRUE)

# read and save the input
seqtab <- readRDS(args[[1]]) # reads sequence table
out <- c(args[[2]])          # a path to the output file
db_loc <- c(args[[3]])       # a path to the SILVA database
db2_loc <- c(args[[4]])      # a path to the SILVA species database

# load libraries
library("dada2")
# assign taxonomy
taxa.1 <- assignTaxonomy(seqtab, db_loc, multithread=TRUE, tryRC=TRUE)
# assign species
taxa <- addSpecies(taxa.1, db2_loc)
# save in rds format
saveRDS(taxa, out)