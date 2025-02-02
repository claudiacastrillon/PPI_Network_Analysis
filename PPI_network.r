library(httr)
library(igraph)
library(readr)
library(stats) # For k-means clustering

# Function to read the proteins from the CSV file
get_proteins_from_csv <- function(file_path) {
  data <- read_csv(file_path)
  print(names(data))  # Print the column names for debugging
  proteins <- data$protname  # Assuming 'protname' is the correct column name for protein symbols
  return(proteins)
}

# Function to fetch PPI data from STRING API and create network
create_ppi_network <- function(proteins) {
  protein_query <- paste(proteins, collapse = "\n")
  encoded_query <- URLencode(protein_query, reserved = TRUE)
  url <- paste0("https://string-db.org/api/tsv/network?identifiers=", encoded_query, "&species=9606")
  
  # Print URL for debugging
  print(paste("Fetching data from URL:", url))
  
  response <- GET(url)
  
  if (status_code(response) == 200) {
    content <- content(response, "text")
    ppi_data <- read.csv(text = content, sep = "\t", header = TRUE)
    return(ppi_data)
  } else {
    stop("Failed to fetch data from STRING API.")
  }
}

# Function to extract features for clustering
extract_features <- function(ppi_data) {
  interaction_counts <- table(c(ppi_data$preferredName_A, ppi_data$preferredName_B))
  feature_matrix <- as.data.frame(as.table(interaction_counts))
  colnames(feature_matrix) <- c("protein", "interaction_count")
  return(feature_matrix)
}

# Function to perform k-means clustering
perform_kmeans_clustering <- function(features, k) {
  set.seed(42)  # For reproducibility
  kmeans_result <- kmeans(features$interaction_count, centers = k)
  features$cluster <- as.factor(kmeans_result$cluster)
  return(features)
}

# Function to plot PPI network with clustering
plot_ppi_network <- function(ppi_data, features) {
  g <- graph_from_data_frame(ppi_data[, c("preferredName_A", "preferredName_B")], directed = FALSE)
  
  # Map clusters to proteins in the graph
  V(g)$cluster <- features$cluster[match(V(g)$name, features$protein)]
  
  # Define colors for clusters
  cluster_colors <- rainbow(length(unique(V(g)$cluster)))
  names(cluster_colors) <- unique(V(g)$cluster)
  
  plot(g, layout = layout_with_fr(g), 
       vertex.color = cluster_colors[V(g)$cluster], 
       vertex.size = 10, 
       vertex.frame.color = "gray", 
       edge.color = "gray50", 
       main = "Protein-Protein Interaction Network with Clustering")
  
  plot(g, layout = layout_with_fr(g), 
       vertex.color = cluster_colors[V(g)$cluster], 
       vertex.size = 15, 
       vertex.label = V(g)$name, 
       vertex.label.color = "black", 
       vertex.label.cex = 0.8, 
       vertex.frame.color = "gray", 
       edge.color = "gray50", 
       main = "Protein-Protein Interaction Network with Clustering and Labels")
}

# Main function to perform analysis
perform_analysis <- function(file_path, k) {
  proteins <- get_proteins_from_csv(file_path)
  if (length(proteins) > 0) {
    ppi_data <- create_ppi_network(proteins)
    if (nrow(ppi_data) > 0) {
      features <- extract_features(ppi_data)
      features <- perform_kmeans_clustering(features, k)
      plot_ppi_network(ppi_data, features)
    } else {
      print("No data available to plot.")
    }
  } else {
    stop("No proteins found in the CSV file.")
  }
}

# Replace with the path to your CSV file and specify the number of clusters
file_path <- "/Users/claudiacastrillonalvarez/Downloads/reactome_R-HSA-913531_proteins.csv"
k <- 3  # Number of clusters
perform_analysis(file_path, k)
