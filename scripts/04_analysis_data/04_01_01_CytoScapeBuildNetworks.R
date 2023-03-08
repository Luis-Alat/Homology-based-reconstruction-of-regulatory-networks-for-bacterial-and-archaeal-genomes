library(RCy3)
library(optparse)
library(tools)

createCytoNets <- function(values = NULL, sleepTime = 2) {
 
  if (values$typeMethod == "new") {
    
    # Loading the extended network (only the values relevants: source, target and category(Known-New interaction))
    
    extenNetwork <- read.table(file = values$inputFile, sep ="\t", header = F)[,c(2,3,4)]

    # Parsing the extended network to get all the nodes and their attributes
    
    NodesID  <- c(extenNetwork[,1],extenNetwork[,2])
    NodesAttributes <- c(extenNetwork[,3],extenNetwork[,3])
    
    nodes <- data.frame(id=NodesID,
                        category=NodesAttributes)
    
    # Getting unique nodes

    NodesUniq <- unique(nodes)
    
    # Using this we get the nodes with "New" category and their names
    
    NodesNew <-  NodesUniq[NodesUniq[,2] == "New",]
    
    # Using this we get the nodes with "Known" category not included in "NodesNew" (No repeated)
    
    NodesKnown<- NodesUniq[!(NodesUniq[,1] %in% NodesNew[,1]),]
      
    # Merging results
    
    nodes <- rbind(NodesNew, NodesKnown)
    
    # Create a new network based in CytoScape
    
    edges <- data.frame(source=extenNetwork[,1],
                        target=extenNetwork[,2],
                        stringsAsFactors=FALSE)
    
    nameCollection <- gsub("_extended_network.txt","",basename(values$inputFile))
    createNetworkFromDataFrames(edges = edges, title = nameCollection, collection = nameCollection)
    
    # Sometimes cytoScape crashs due to R instructions are executed faster than cytoScape instructions
    Sys.sleep(sleepTime)
    
    # Load an additional column in the current table including "nodes" data ("category" data)
    
    loadTableData(nodes, data.key.column =  "id")
    
    #Set up an style for the new interactions
    
    style <- "newInteractionsStyle"
    defaults <- list(NODE_SHAPE="round_rectangle", NODE_BORDER_WIDTH=1.0)
    nodeLabels <- mapVisualProperty("node label", "name","p")
    nodeFills <- mapVisualProperty('node fill color','category','d',c("Known","New"), c("#89D0F5","#FA554A"))
    createVisualStyle(style.name =  style, defaults = defaults, mappings = list(nodeLabels, nodeFills))
    
    # Apply style
    
    print("Applying personlized style...")
    
    setVisualStyle(style.name = style)
    
    # Sometimes cytoScape crashs due to R instructions are executed faster than cytoScape instructions
    Sys.sleep(sleepTime)
    
    # Analyze network as directed graph, get summary and statistic table
    
    print("Analizying network...")
    
    SUMMARY <- analyzeNetwork(TRUE)
    
    # Sometimes cytoScape crashs due to R instructions are executed faster than cytoScape instructions
    Sys.sleep(sleepTime)
    
    SUMMARY <- as.matrix(SUMMARY)
    rownames(SUMMARY) <- paste("#", rownames(SUMMARY), sep = "")
    NodeStatistics <- getTableColumns(table = "node")
    NodeStatistics <- NodeStatistics[order(NodeStatistics$Outdegree, decreasing = T),]
    
  } else{
    
    stop("Invalid method\n")
    
  }

  OutPutFileName <- tools::file_path_sans_ext( basename(values$inputFile) )
  OutPutFileName <- paste(values$path, OutPutFileName, sep = "", collapse = "")
  
  # CytoScape has problems handling files with points (ONLY TO SAVE IMAGES!)
  
  #fileNameException <- gsub("\\.", "_", OutPutFileName)
  
  # Saving session, images and statistics

  print("Saving files...")

  # saveSession and exportImage are not working properly as used to; therefore, they are going to be saved by using CytoScape commands instead of Rcy3 sintax
  # exportImage(filename = OutPutFileName, type = 'SVG', overwriteFile = TRUE)
  
  #commandSaveSvg=paste("view export options=\"SVG (*.svg)\" view=CURRENT outputFile=", OutPutFileName, sep = "", collapse = "")
  #commandsRun(commandSaveSvg)
  RCy3::exportImage(filename = OutPutFileName, type = "SVG")
  
  #saveSession(OutPutFileName)
  #commandSaveSession=paste("session save as file=", OutPutFileName, sep = "", collapse = "")
  #commandsRun(commandSaveSession)  
  RCy3::saveSession(filename = OutPutFileName)
  
  # Saving network statistic
  write.table(x = SUMMARY, file = paste(OutPutFileName, "_metrics", sep = ""),
              quote = F,sep = "\t",row.names = T, col.names = F)
  write.table(x = NodeStatistics, file = paste(OutPutFileName, "_metrics", sep = ""),
              quote = F,sep = "\t",row.names = F, col.names = T, append = T)
  
  closeSession(save.before.closing = F)
  cytoscapeFreeMemory()  
  Sys.sleep(sleepTime)
  
  print("All done")

}


# Check cytoscape connection, user must launch it
cytoscapePing()

# Retrieve bash parameter

# Method: string describing if the networks will be processed as new interactions or orthologs
# Output files will be created using the basename of "--inputFile"

option_list = list(
  make_option(c("-i", "--inputFile"), type="character", default=NULL, 
              help="Extended Network or information about orthologs", metavar="character"),
  make_option(c("-t", "--typeMethod"), type = "character", default = NULL,
              help="Type of information to process: new or ortho", metavar="character"),
  make_option(c("-p", "--path"), type = "character", default = NULL,
              help="Path to save files", metavar = "character")
);

opt_parser = OptionParser(option_list=option_list);
parameters = parse_args(opt_parser);

print(parameters)


# Examples of values in "parameters" running in a "ortho" type method

#parameters$inputFile <- "/home/lromero/Documentos/IIMAS/IIMAS_V3/net/results/GCF_000009045.1_ASM904v1_B_subtilis_168_genomic_Orto_Summary.txt"
#parameters$modifiedNetwork <- "/home/lromero/Documentos/IIMAS/IIMAS_V3/modifiedNetworks/GCF_000009045.1_ASM904v1_B_subtilis_168_modified_net.txt"
#parameters$typeMethod <- "ortho"
#parameters$colors <- "#89D0F5,#9D6AC7,#97FA4A,#EFFA4A,#FAAD4A,#FA554A"
#parameters$path <- "/home/lromero/" # ending "/" is important
#values <- parameters

createCytoNets(values = parameters)
