library(igraph)
library(rtracklayer)
library(GenomicFeatures)
library(SplicingGraphs)
library(DEXSeq)
library(stringr)
library(parallel)

options(scipen=50)

from_bubble_to_exon_path <- function(g, bubblepath, node2pos) {
  path = c()
  for (i in 1:(length(bubblepath)-1)) {
    # weird way to get all the multi edges between two vertices
    edges = E(g)[node2pos[bubblepath[i]] %--% node2pos[bubblepath[i+1]]] 
    # keep only exon edges and not the exonic part edges. 
    # unfortunately, if we use get.edge.ids approach above it will always give exonic part edges instead
    edges = edges[edge_attr(g)$ex_or_in[as.numeric(edges)] != "None"]
    path = c(path, as.numeric(edges))	
  }
  return(path)
}



from_exon_path_to_exonic_part_path <- function(g, exonpath) {
  exon_path = c()
  exonic_part_path = c()
  vpaths_list = c()
  
  for (i in 1:length(exonpath)) {
    edge = exonpath[i]
    if (edge_attr(g)$ex_or_in[edge] == "ex") {
      start = ends(g, exonpath)[i, 1]
      end = ends(g, exonpath)[i, 2]
      for (x in 1:(as.numeric(V(g)[end]$id) - as.numeric(V(g)[start]$id))){
        vpaths_list <- c(vpaths_list, all_simple_paths(g, as.numeric(V(g)[start]$id) + x, as.numeric(V(g)[start]$id) + x + 1, cutoff=2))
      }
      for (x in 1:length(vpaths_list)){
        exonic_edge = E(g)[vpaths_list[[x]][1] %--% vpaths_list[[x]][2]]
        exonic_part_path = c(exonic_part_path, exonic_edge[edge_attr(g)$ex_or_in[exonic_edge] == "None"])
      }
    }
    else if (edge_attr(g)$ex_or_in[edge] == "in"){
      exon_path = c(exon_path, edge)
      exonic_part_path = c(exonic_part_path, edge)
    }
    else {
      print("should not come here. should only have ex or in in path.") 
    }
  }
  #stopifnot(all(exonpath == exon_path))
  return(exonic_part_path)
}



collapse_binary_bubbles <- function(complex_bubbles, binary_bubbles, graph, bubble_data, node2pos) {
  
  if (length(names(complex_bubbles)) == 0 | length(names(binary_bubbles)) == 0) {
    return (list(complex_bubbles, binary_bubbles, graph, bubble_data))
  }
  
  transcripts <- edge_attr_names(graph)[grepl("ENST", edge_attr_names(graph))]
  exonic_paths <- c()
  
  for (y in names(complex_bubbles)) {
    if (length(names(binary_bubbles[as.integer(sapply(strsplit(names(binary_bubbles), "-"), "[[", 1)) >= as.integer(unlist(strsplit(y, "-"))[1]) & 
                                    as.integer(sapply(strsplit(names(binary_bubbles), "-"), "[[", 2)) <= as.integer(unlist(strsplit(y, "-"))[2])])) == 0) {
      #print(paste0("Complex Bubble ", y, " contains no Binary Bubbles"))
      next
    }
    
    lowest_level = TRUE
    
    for (bubble_name in names(complex_bubbles[as.integer(sapply(strsplit(names(complex_bubbles), "-"), "[[", 1)) >= as.integer(unlist(strsplit(y, "-"))[1]) & 
                                              as.integer(sapply(strsplit(names(complex_bubbles), "-"), "[[", 2)) <= as.integer(unlist(strsplit(y, "-"))[2]) & 
                                              y != names(complex_bubbles)])) {
      
      if(length(names(binary_bubbles[as.integer(sapply(strsplit(names(binary_bubbles), "-"), "[[", 1)) <= as.integer(unlist(strsplit(bubble_name, "-"))[1]) & 
                                     as.integer(sapply(strsplit(names(binary_bubbles), "-"), "[[", 2)) >= as.integer(unlist(strsplit(bubble_name, "-"))[2])])) == 0) {
        lowest_level = FALSE
        break
      }
    }
    
    if (!lowest_level) {next}
    
    for (bubble_name in names(binary_bubbles[as.integer(sapply(strsplit(names(binary_bubbles), "-"), "[[", 1)) >= as.integer(unlist(strsplit(y, "-"))[1]) & 
                                             as.integer(sapply(strsplit(names(binary_bubbles), "-"), "[[", 2)) <= as.integer(unlist(strsplit(y, "-"))[2])])) {
      # cat('\n')
      # print(paste0("Complex Bubble ", y, " contains Binary Bubble ", bubble_name))
      bubbles <- strsplit(binary_bubbles[[bubble_name]]$Bubbles, " ")[[1]]

      for (i in 1:length(bubbles)) {
        nodes = strsplit(bubbles[i], ",")[[1]]
        # cat('\n')
        # print(paste0("Edges to Delete in Path ", bubbles[i],": "))
        
        for (node in 1:(length(nodes) - 1)) {
          edge <- E(graph)[node2pos[nodes[node]] %--% node2pos[nodes[node+1]]]
          edge <- edge[edge_attr(graph, "ex_or_in", edge) != "None" & is.na(edge_attr(graph, "exonic_weight", edge))]
          
          # for inner nodes (not first/source and last/sink), delete binary nodes from complex bubbles to collapse
          if (node > 1) {
            bubbles <- gsub(paste0(nodes[node], ','), '', bubbles)
          }
          
          if (length(edge) == 0) {next}
          
          #print(paste0(nodes[node], " -> ", nodes[node+1], ": edge ", as.numeric(edge)," with edge attribute ", edge_attr(graph, "ex_or_in", edge)))
            
          for (exonic_node in as.numeric(nodes[node]):(as.numeric(nodes[node+1]) - 1)) {
            exonic_edge <- E(graph)[node2pos[as.character(exonic_node)] %--% node2pos[as.character(exonic_node + 1)]]
            exonic_edge <- exonic_edge[!is.na(edge_attr(graph, "exonic_weight", exonic_edge))]
            if (length(exonic_edge) == 0){
              graph <- graph + edge(node2pos[as.character(exonic_node)], node2pos[as.character(exonic_node + 1)], 
                                    exonic_weight = 0, ex_or_in = edge_attr(graph, "ex_or_in", as.numeric(edge)), dexseq_fragment = "")
              exonic_edge <- E(graph)[gsize(graph)]
              #print(paste0("Created edge ", exonic_edge, " with ends ", paste(V(graph)[ends(graph, exonic_edge)]$id, collapse = ",")))
            }
            else{
              #print(paste0("Found edge ", exonic_edge, " with ends ", paste(V(graph)[ends(graph, exonic_edge)]$id, collapse = ",")))
            }
            
            if (edge_attr(graph, "ex_or_in", as.numeric(edge)) == "ex") {
              edge_attr(graph, "exonic_weight", exonic_edge) <- edge_attr(graph, "exonic_weight", exonic_edge) + 1
              edge_attr(graph, "ex_or_in", exonic_edge) <- "ex"
            }
                
            # copy transcript information to exonic part (only need to worry about transcripts that are true in edge)
            for (transcript in transcripts) {
              edge_attr(graph, transcript, exonic_edge) <- edge_attr(graph, transcript, exonic_edge) | edge_attr(graph, transcript, as.numeric(edge))
            }
            
            # if (exonic_node == as.numeric(nodes[node+1]) - 1) {
            #   exonic_path <- append(exonic_path, as.character(exonic_node))
            #   exonic_path <- append(exonic_path, as.character(exonic_node+1))
            # }
            # else {
            #   exonic_path <- append(exonic_path, as.character(exonic_node))
            # }
          }
            
          # delete edge from graph, collapsing it
          #print(paste0("Deleting edge ", edge, " with ends ", paste(V(graph)[ends(graph, as.numeric(edge))]$id, collapse = ",")))
          graph <- delete_edges(graph, as.numeric(edge))
        }
      }
      
        
      # delete redundant entries in bubbles
      bubbles <- unique(bubbles)
      
      # delete collapsed binary bubbles
      if(length(bubbles) == 1) {
        # cat('\n')
        # print(paste0("Binary Bubble ", bubble_name, " has been collapsed to a Flattened Path"))
        binary_bubbles[[bubble_name]] <- NULL
      }
      
      #cat("\n")
      #print(paste0("Flattened Exonic Path: ", paste(exonic_path, collapse = " -> ")))
      exonic_paths <- append(exonic_paths, bubbles)
      
      # print("Updating Paths of All Complex Bubbles Affected...")
      # cat('\n')
      
      # search through relevant complex bubble paths and delete path segments that do not exist in the graph anymore
      for (c_bubble_name in names(complex_bubbles[as.integer(sapply(strsplit(names(complex_bubbles), "-"), "[[", 1)) <= as.integer(unlist(strsplit(bubble_name, "-"))[1]) & 
                                                  as.integer(sapply(strsplit(names(complex_bubbles), "-"), "[[", 2)) >= as.integer(unlist(strsplit(bubble_name, "-"))[1])])) {
        c_bubbles <- strsplit(complex_bubbles[[c_bubble_name]]$Bubbles, " ")[[1]]
        #cat('\n')
        #print(paste0("For Complex Bubble ", c_bubble_name, "..."))
        
        for (i in 1:length(c_bubbles)) {
          nodes = strsplit(c_bubbles[i], ",")[[1]]
          if (length(nodes) == 2) {next}
          
          #cat('\n')
          #print(paste0("Traversing Complex Path ", paste(nodes, collapse = ",")))
          for (node in 2:(length(nodes)-1)) {
            counter = 0
            #print(paste0("Edges for node ", nodes[node], ": ", nodes[node-1], " -> ", nodes[node], " and ", nodes[node], " -> ", nodes[node+1]))
            for (j in 0:1) {
              edge <- E(graph)[node2pos[nodes[node-1 + j]] %--% node2pos[nodes[node + j]]]
              edge <- edge[edge_attr(graph, "ex_or_in", edge) != "None" & is.na(edge_attr(graph, "exonic_weight", edge))]
              if (length(edge) == 0){
                #print(paste0("Edge ", j+1, " does not exist"))
                counter = counter + 1
              }
            }
            
            if (counter == 2){
              #print(paste0("We can delete ", nodes[node], " from ", c_bubbles[i]))
              c_bubbles[i] <- sub(paste0(as.character(nodes[node]), ","), "", c_bubbles[i])
            }
            
            #print(paste0("Outcome: ", c_bubbles[i]))
          }
        }
        
        # delete redundant entries in c_bubbles
        c_bubbles <- unique(c_bubbles)
        
        paths_to_delete <- c()
        
        for (i in 1:length(c_bubbles)) {
          if (length(strsplit(c_bubbles[i], ",")[[1]]) == 2 & length(c_bubbles) > 2) {
            #cat('\n')
            #print(paste0("Complex Bubble ", c_bubble_name, " has > 2 paths and Path ", c_bubbles[i], " only contains the start and end nodes, so we can delete this entry"))
            paths_to_delete <- append(paths_to_delete, c_bubbles[i])
          }  
        }
        
        c_bubbles <- c_bubbles[!(c_bubbles %in% paths_to_delete)]
        
        # update complex bubble
        complex_bubbles[[c_bubble_name]]$Bubbles <- paste(unlist(c_bubbles), collapse = " ")
        #cat('\n')
        #print(paste0("Complex Bubble ", c_bubble_name," is now ", complex_bubbles[[c_bubble_name]]$Bubbles))
        if (length(strsplit(complex_bubbles[[c_bubble_name]]$Bubbles, " ")[[1]]) == 1) {
          exonic_paths <- append(exonic_paths, complex_bubbles[[c_bubble_name]]$Bubbles)
        }
      }
      
      # cat('\n')
      # print("Done Updating Paths of Affected Bubbles...")
    }
    
    #cat("\n\n")
  }
  
  
  exonic_paths <- unique(exonic_paths)
  
  # replace sections of each relevant complex bubble that matches with any of the exonic paths
  # the flattened path can be used for finding focal exons
  for (exonic_path in exonic_paths) {
    # cat('\n')
    # print(paste0("For Flattened Path ", exonic_path, "..."))
    start_node <- strsplit(exonic_path, ",")[[1]][1]
    end_node <- strsplit(exonic_path, ",")[[1]][2]
    
    for (c_bubble_name in names(complex_bubbles[as.integer(sapply(strsplit(names(complex_bubbles), "-"), "[[", 1)) <= as.integer(start_node) & 
                                                as.integer(sapply(strsplit(names(complex_bubbles), "-"), "[[", 2)) >= as.integer(start_node)])) {
      
      c_bubbles <- strsplit(complex_bubbles[[c_bubble_name]]$Bubbles, " ")[[1]]
      # cat('\n')
      # print(paste0("For Complex Bubble ", c_bubble_name, " - ", c_bubbles))
      
      for (i in 1:length(c_bubbles)) {
        flattened_path = as.character(start_node:end_node)
        binary_nodes = strsplit(exonic_path, ",")[[1]]
        nodes = strsplit(c_bubbles[i], ",")[[1]]
        
        # cat('\n')
        # print(paste0("Searching Complex Path ", paste(nodes, collapse = ",")))
        for (node in 1:(length(nodes)-1)) {
          if (nodes[node] == binary_nodes[1] & nodes[node+1] == binary_nodes[2]) {
            # print(paste0("Complex Nodes ", nodes[node], ",", nodes[node+1],
            #             " match the nodes from the Binary Path ", paste(binary_nodes, collapse = ",")))
            # print(paste0("Complex Nodes will be replaced with Flattened Exonic Path ", paste(flattened_path, collapse = ",")))
            c_bubbles[i] <- sub(paste0(nodes[node], ",", nodes[node+1]), paste(flattened_path, collapse = ","), c_bubbles[i])
          }
          
          if (length(E(graph)[node2pos[nodes[node]] %--% node2pos[nodes[node+1]]]) == 0) {
            # print(paste0("Edge ", nodes[node], "-", nodes[node+1], 
            #              " does not exist, so it will be replaced with a Flattened Path - ", 
            #              paste(nodes[node]:nodes[node+1], collapse = ",")))
            c_bubbles[i] <- sub(paste0(nodes[node], ",", nodes[node+1]), 
                                paste(nodes[node]:nodes[node+1], collapse = ","), c_bubbles[i])
          }
          
          #print(paste0("Result: ", c_bubbles[i]))
        }
      }
      
      # update complex bubble
      # cat("\n")
      # print(paste0("Complex Bubble ", c_bubble_name, " is now ", paste(c_bubbles, collapse = " ")))
      
      # if collapsed into flattened path, delete the complex bubble
      if (length(c_bubbles) == 1) {
        #print(paste0("Complex Bubble ", c_bubble_name, " has been collapsed to a Flattened Path"))
        complex_bubbles[[c_bubble_name]] <- NULL
      }
      # if collapsed into binary, copy complex bubble into binary_bubbles and delete
      else if (length(c_bubbles) == 2) {
        #print(paste0("Complex Bubble ", c_bubble_name, " has been collapsed to a Binary Bubble"))
        binary_bubbles[[c_bubble_name]] <- complex_bubbles[[c_bubble_name]]
        complex_bubbles[[c_bubble_name]] <- NULL
        
        #find focal exons of the new binary bubble and add it to bubble data
        
        bubble1 <- c_bubbles[1]
        bubble1 <- unlist(strsplit(bubble1, ","))
        
        bubble2 <- c_bubbles[2]
        bubble2 <- unlist(strsplit(bubble2, ","))
        
        bubble_data$BubbleType <- append(bubble_data$BubbleType, "Complex -> Flattened Binary")
        bubble_data$Bubbles <- append(bubble_data$Bubbles, binary_bubbles[[c_bubble_name]]$Bubbles)
        
        # print(paste0("Finding Focal Exons for new Binary Path"))
        # print(paste0("Bubble1: ", paste(bubble1, collapse = ","), " and Bubble2: ", paste(bubble2, collapse = ",")))
        
        path1 = from_bubble_to_exon_path(graph, bubble1, node2pos)
        path2 = from_bubble_to_exon_path(graph, bubble2, node2pos)
        
        # delete any exons that may be double counted as exonic parts, but are part of a different transcript path
        paths1_to_remove = c()
        paths2_to_remove = c()
        
        for (path in path1) {
          if (abs(as.numeric(V(graph)[ends(graph, path)]$id[1]) - as.numeric(V(graph)[ends(graph, path)]$id)[2]) == 1) {
            if (is.na(edge_attr(graph, "exonic_weight", path))) {
              paths1_to_remove <- append(paths1_to_remove, path)
            }
          }
        }
        
        for (path in path2) {
          if (abs(as.numeric(V(graph)[ends(graph, path)]$id[1]) - as.numeric(V(graph)[ends(graph, path)]$id)[2]) == 1) {
            if (is.na(edge_attr(graph, "exonic_weight", path))) {
              paths2_to_remove <- append(paths2_to_remove, path)
            }
          }
        }
        
        path1 <- path1[!(path1 %in% paths1_to_remove)]
        path2 <- path2[!(path2 %in% paths2_to_remove)]
        
        # convert exon paths into exonic part paths of the graph
        e_path1 = from_exon_path_to_exonic_part_path(graph, path1)
        e_path2 = from_exon_path_to_exonic_part_path(graph, path2)
        
        # find difference in paths
        focalexons = union(setdiff (e_path1, e_path2), setdiff (e_path2, e_path1))
        # exclude all introns
        focalexons = focalexons[edge_attr(graph)$ex_or_in[focalexons]  != "in"]
        
        bubble_data$FocalExons <- append(bubble_data$FocalExons, paste0("E", paste(unlist(edge_attr(graph)$dexseq_fragment[focalexons]), collapse=",E")))
        #print(paste0("Focal Exons Found: ", paste0("E", paste(unlist(edge_attr(graph)$dexseq_fragment[focalexons]), collapse=",E"))))
      }
    }
  }
  
  
  results_list <- collapse_binary_bubbles(complex_bubbles, binary_bubbles, graph, bubble_data, node2pos) 
  complex_bubbles <- results_list[[1]]
  binary_bubbles <- results_list[[2]]
  graph <- results_list[[3]]
  bubble_data <- results_list[[4]]
  
  return (list(complex_bubbles, binary_bubbles, graph, bubble_data))
}










# parallelize the script
numCores <- detectCores()

genes <- list.dirs("/home/dwito/BonnalData/Bnaive_CD8naive/grase/grase_results/gene_files", recursive=F, full.names=T)
genes <- genes[17]

#use mclapply for parallel processing of all genes in the dataset
mclapply(genes, function(x) {
  gene = sapply(strsplit(x, split='/')[1], tail, 1)
  print(gene)
  gene_graphml = paste0(x, "/output/", gene, ".graphml")
  gtf = paste0(x, "/", gene, ".gtf")
  
  
  g <- read_graph(gene_graphml, format = c("graphml"))
  
  vertices <- V(g) #vertices
  edges <- E(g) #edges
  vertex_labels <- vertex_attr(g) #get the vertices of graph
  
  gffFile = gtf
  gr <- import(gffFile)
  gr <- gr[!(mcols(gr)$type %in% c("start_codon", "stop_codon"))]
  transcript_id <- mcols(gr)$transcript_id
  txdb <- makeTxDbFromGRanges(gr)
  sg <- SplicingGraphs(txdb,  min.ntx=1)
  edges_by_gene <- sgedgesByGene(sg)
  
  # this whole code is to translate between node IDs and names (positions). 
  if (strand(sg) == '+'){
    node2pos1 = cbind(edges_by_gene[[gene]]$to, end(edges_by_gene[[gene]])+1)
    node2pos2 = cbind(edges_by_gene[[gene]]$from, start(edges_by_gene[[gene]]))
    node2pos = rbind(node2pos1, node2pos2)
    node2pos <- unique(node2pos[order(node2pos[,1]),])
    rownames(node2pos) <- node2pos[,1]
    node2pos <- node2pos[,2]
  }
  
  if (strand(sg) == '-'){
    node2pos1 = cbind(edges_by_gene[[gene]]$from, end(edges_by_gene[[gene]])+1)
    node2pos2 = cbind(edges_by_gene[[gene]]$to, start(edges_by_gene[[gene]]))
    node2pos = rbind(node2pos1, node2pos2)
    node2pos <- unique(node2pos[order(node2pos[,1]),])
    rownames(node2pos) <- node2pos[,1]
    node2pos <- node2pos[,2]
  }
  
  bubbles_df = bubbles(sg[gene])
  
  
  bubble_data = list() # holds tabular data for all bubbles in the gene
  binary_bubbles <- c() # named list that will contain information for each binary bubble in the gene
  binary_bubbles_index <- c()
  complex_bubbles <- c() # named list that will contain information for each complex bubble in the gene
  complex_bubble_index <- c()
  
  # if the gene has bubble paths, process it
  if (length(bubbles_df$paths) > 0){
    for (bubble in 1:length(bubbles_df$paths)){
      # process binary cases separately from more complex cases
      if (length(bubbles_df$paths[[bubble]]) == 2 & 
          bubbles_df$source[bubble] != 'R' & bubbles_df$source[bubble] != 'L' & 
          bubbles_df$sink[bubble] != 'R' & bubbles_df$sink[bubble] != 'L'){
        # binary bubbles
        bubble1 = unlist(c(bubbles_df$source[bubble], strsplit(gsub('[{}]','',bubbles_df$paths[[bubble]])[1], ","), bubbles_df$sink[bubble]))
        bubble2 = unlist(c(bubbles_df$source[bubble], strsplit(gsub('[{}]','',bubbles_df$paths[[bubble]])[2], ","), bubbles_df$sink[bubble]))

        bubble_data$BubbleType <- append(bubble_data$BubbleType, "Binary")
        binary_bubble_path <- paste0(paste(unlist(bubble1), collapse=","), " ",  paste(unlist(bubble2), collapse=","))
        bubble_data$Bubbles <- append(bubble_data$Bubbles, binary_bubble_path)
        
        # convert bubbles into exon paths of the graph
        path1 = from_bubble_to_exon_path(g, bubble1, node2pos)
        E(g)[path1]
        
        path2 = from_bubble_to_exon_path(g, bubble2, node2pos)
        E(g)[path2]
        
        # convert exon paths into exonic part paths of the graph
        e_path1 = from_exon_path_to_exonic_part_path(g, path1)
        e_path2 = from_exon_path_to_exonic_part_path(g, path2)
        
        # find difference in paths
        focalexons = union(setdiff (e_path1, e_path2), setdiff (e_path2, e_path1))
        # exclude all introns
        focalexons = focalexons[edge_attr(g)$ex_or_in[focalexons]  != "in"]
        bubble_data$FocalExons <- append(bubble_data$FocalExons, paste0("E", paste(unlist(edge_attr(g)$dexseq_fragment[focalexons]), collapse=",E")))
        
        # store bubble information into binary_bubbles data structure
        binary_bubbles_index <- append(binary_bubbles_index, bubble)
        binary_bubble_name <- paste0(bubbles_df$source[[bubble]], "-", bubbles_df$sink[[bubble]])
        
        binary_bubbles[[binary_bubble_name]]$Bubbles <- append(binary_bubbles[[bubbles_df$source[[bubble]]]]$Bubbles, 
                                                                        binary_bubble_path)
        # binary_bubbles[[bubbles_df$source[[bubble]]]]$Transcripts <- append(binary_bubbles[[bubbles_df$source[[bubble]]]]$Transcripts, 
        #                                                                     gsub("[{}]", "", paste(unlist(bubbles_df$partitions[[bubble]]), collapse=" ")))
        # binary_bubbles[[binary_bubble_name]]$Sink <- append(binary_bubbles[[bubbles_df$source[[bubble]]]]$Sink, 
        #                                                               bubbles_df$sink[[bubble]])
        # binary_bubbles[[bubbles_df$source[[bubble]]]]$FocalExons <- append(binary_bubbles[[bubbles_df$source[[bubble]]]]$FocalExons, 
        #                                                                    paste0("E", paste(unlist(edge_attr(g)$dexseq_fragment[focalexons]), collapse=",E")))
        
        # # Obtain exon that is left of focal exon
        # in_node <- incident(g, ends(g, focalexons)[1], mode = "in")
        # # Go to the next vertex if it doesn't have an exonic part
        # if(all(edge_attr(g)$dexseq_fragment[in_node] == "")){
        #   in_node <- incident(g, ends(g, in_node)[1], mode = "in")
        # }
        # # Obtain exon that is right of focal exon
        # out_node <- incident(g, ends(g, focalexons)[2], mode = "out")
        # # go to the next vertex if it doesn't have an exonic part
        # if(all(edge_attr(g)$dexseq_fragment[out_node] == "")){
        #   out_node <- incident(g, ends(g, out_node)[2], mode = "out")
        # }
        # 
        # adjexon1 <- in_node[edge_attr(g)$ex_or_in[in_node]  == "None"] 
        # adjexon2 <- out_node[edge_attr(g)$ex_or_in[out_node]  == "None"]
      }
      
      # process complex cases (> 2 bubbles)
      else if (length(bubbles_df$paths[[bubble]]) > 2 & 
               bubbles_df$source[bubble] != 'R' & bubbles_df$source[bubble] != 'L' & 
               bubbles_df$sink[bubble] != 'R' & bubbles_df$sink[bubble] != 'L'){
        
        bubble_data$BubbleType <- append(bubble_data$BubbleType, "Complex")
        
        # process the complex bubble and get all paths
        complex_bubble_path <- unlist(c(bubbles_df$source[bubble], strsplit(gsub('[{}]','',bubbles_df$paths[[bubble]])[1], ","), bubbles_df$sink[bubble]))
        complex_bubble_path <- paste0(complex_bubble_path, collapse=",")
        for (i in 2:length(bubbles_df$paths[[bubble]])){
          bubble1 = unlist(c(bubbles_df$source[bubble], strsplit(gsub('[{}]','',bubbles_df$paths[[bubble]])[i], ","), bubbles_df$sink[bubble]))
          complex_bubble_path <- paste0(complex_bubble_path, " ", paste(unlist(bubble1), collapse=","))
        }
        

        bubble_data$Bubbles <- append(bubble_data$Bubbles, complex_bubble_path)
        bubble_data$FocalExons <- append(bubble_data$FocalExons, "N/A")
        
        # store bubble data into complex_bubbles data structure
        complex_bubble_index <- append(complex_bubble_index, bubble)
        complex_bubble_name <- paste0(bubbles_df$source[[bubble]], "-", bubbles_df$sink[[bubble]])
          
        complex_bubbles[[complex_bubble_name]]$Bubbles <- append(complex_bubbles[[bubbles_df$source[[bubble]]]]$Bubbles, 
                                                                         complex_bubble_path)
        # complex_bubbles[[bubbles_df$source[[bubble]]]]$Transcripts <- append(complex_bubbles[[bubbles_df$source[[bubble]]]]$Transcripts, 
        #                                                                      gsub("[{}]", "", paste(unlist(bubbles_df$partitions[[bubble]]), collapse="  ")))
        # complex_bubbles[[complex_bubble_name]]$Sink <- append(complex_bubbles[[bubbles_df$source[[bubble]]]]$Sink, 
        #                                                                 bubbles_df$sink[[bubble]])
      }
    }
    
    if (length(bubble_data) > 0) {
      # begin processing to collapse binary bubbles
      # in processing binary bubbles, some complex bubbles will collapse into binary cases
      # collapse these new binary cases until they are exhausted
      collapsed_graph <- read_graph(gene_graphml, format = "graphml") # copy of graph that will be collapsed
      edge_attr(collapsed_graph, "exonic_weight", ) <- NA # stores weights of collapsed exons into their constituent exonic parts
      
      # collapse binary bubbles, then output updated data structures and graph
      results_list <- collapse_binary_bubbles(complex_bubbles, binary_bubbles, collapsed_graph, bubble_data, node2pos)
      complex_bubbles <- results_list[[1]]
      binary_bubbles <- results_list[[2]]
      collapsed_graph <- results_list[[3]]
      bubble_data <- results_list[[4]]
      
      table <- as.data.frame(bubble_data)
      write.table(table, file=paste0("/home/dwito/BonnalData/Bnaive_CD8naive/grase/grase_results/gene_files/",
                                     gene, "/output/", gene, ".FocalExons.v3.txt"),
                  sep='\t', quote=F, row.names=F)
      write_graph(collapsed_graph, file = paste0("/home/dwito/BonnalData/Bnaive_CD8naive/grase/grase_results/gene_files/",
                                                 gene, "/output/", gene, ".collapsed.graphml"))
    }
  }
}, mc.cores = 1, mc.preschedule = T)

