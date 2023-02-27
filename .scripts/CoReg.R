library(devtools)

list.of.packages <- c("CoReg")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install_github("LiLabAtVT/CoReg")

library(CoReg)
library(optparse)

option_list = list(
  make_option(c("-i", "--inputFile"), type="character", default=NULL, 
              help="Extended Network", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
parameters = parse_args(opt_parser);

print(parameters)
print("Loading data...")

# Loading only the TF and target columns

ExtendedNet <- read.table(parameters$inputFile, header = F, sep = "\t")[,c(2,3)]

# Converting into a iGraph our two columns format

ExtendedNetIGraph <- networkFromEdgeList(ExtendedNet)

# Getting and saving co-regulators results

print("Getting co-regulators")

ExtNetResults <- CoReg(ExtendedNetIGraph)

print("Saving files...")

outputFileCoRegulators <- gsub(".txt", "_Modules.txt", basename(parameters$inputFile))
outputFileMatrix <- paste("tmp/", gsub(".txt", "_Matrix.txt", basename(parameters$inputFile)), sep = "")
outputFileRange <-  paste("tmp/", gsub(".txt", "_Range.txt",  basename(parameters$inputFile)), sep = "")

# Matrix and Range will be saved in the tmp folder, Co-regulators in the current one

write.table(x = ExtNetResults$module, file = outputFileCoRegulators,
            quote = F, sep = "\t", col.names = T, row.names = F)
write.table(x = ExtNetResults$similarity_matrix, file = outputFileMatrix ,
            quote = F, sep = "\t", col.names = T, row.names = T)
write.table(x = ExtNetResults$rank, file = outputFileRange, 
            quote = F, sep = "\t", col.names = T, row.names = F)

print("All done")