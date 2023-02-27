library(easyPubMed)
library(optparse)

makeQueries <- function(query, fileName, retmax = 10){
  
  Message <- paste("Searching", query)
  print(Message)

  my_entrez_id <- get_pubmed_ids(query)
  
  if(as.integer(my_entrez_id$Count) < 1){
    warningMessage <- paste("   ", query, "returned zero results") 
    print(warningMessage)
    return(0)
  }
  
  my_art <- fetch_pubmed_data(my_entrez_id, format = "abstract", retmax = retmax)

  header_search = paste("\n#############", query ,"#############\n")
  
  write(x = header_search, file = fileName, append = TRUE)
  write(x = my_art, file = fileName, append = TRUE)
  
}

option_list = list(
  make_option(c("-i", "--inputFile"), type="character", default=NULL, 
              help="Rules or query to search on PubMed", metavar="character"),
  make_option(c("-o", "--outputFile"), type="character", default=NULL,
               help="Output file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
parameters = parse_args(opt_parser);

queries <- readLines(parameters$inputFile)
fileName <- parameters$outputFile

print("##### PARAMETERS #####")
print(paste("Input file name:", parameters$inputFile))
print(paste("Output file name:", parameters$outputFile))

# Empty file to save query results

write(x = "", file = fileName)

for(i in c(1:length(queries))) {
  makeQueries(query = queries[i], fileName = fileName)
}
