library(easyPubMed)
library(optparse)

MakeQueries <- function(query, file_name, retmax = 10) {

  message <- paste("Searching: ", query)
  print(message)

  my_entrez_id <- easyPubMed::get_pubmed_ids(query)
  number_matches <- as.integer(my_entrez_id$Count)

  if (number_matches < 1) {

    warning_message <- paste("   ", query, "returned zero results")
    print(warning_message)

    return(0)

  }

  message_matches <- paste("   ", query, " returned ",
                          number_matches, " results", sep = "")

  print(message_matches)

  my_art <- easyPubMed::fetch_pubmed_data(my_entrez_id, format = "abstract",
                              retmax = retmax)

  header_search <- paste("\n#############", query, "#############\n")

  write(x = header_search, file = file_name, append = TRUE)
  write(x = my_art, file = file_name, append = TRUE)

}

option_list <- list(
  make_option(c("-i", "--inputFile"), type = "character", default = NULL,
              help = "Rules or query to search on PubMed",
              metavar = "character"),
  make_option(c("-o", "--outputFile"), type = "character", default = NULL,
              help = "Output path and file name",
              metavar = "character")
);

opt_parser <- OptionParser(option_list = option_list);
parameters <- parse_args(opt_parser);

queries <- readLines(parameters$inputFile)

file_name <- paste(parameters$outputFile, parameters$inputFile,
                  "_pubmed", sep = "", collapse = "")

print("##### PARAMETERS #####")
print(paste("Input file name:", parameters$inputFile))
print(paste("Output path file name:", file_name))

# Empty file to save query results

write(x = "", file = file_name)

for (i in c(1:length(queries))) {
  MakeQueries(query = queries[i], file_name = file_name)
}
