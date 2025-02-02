# 🧬 Protein-Protein Interaction Network Analysis
📌 Overview

This repository contains an R script for constructing and analyzing Protein-Protein Interaction (PPI) networks using data from PubChem and the STRING database. The script retrieves interaction data, constructs a network, and applies k-means clustering for protein grouping.

🔬 Methodology

Data Acquisition

Retrieves protein lists from PubChem based on Reactome IDs.

Queries the STRING database using httr in RStudio.

Loads interaction data in .tsv format.

Network Construction

Uses igraph to construct an undirected network.

Nodes represent proteins, edges represent interactions.

Weights edges by interaction confidence scores.

Clustering Analysis (K-Means)

Uses stats::kmeans to cluster proteins.

Positions nodes using Spring Layout (force-directed graph drawing).

Assigns cluster-based colors for visualization.

📦 Dependencies

Ensure you have the following R libraries installed:

install.packages(c("httr", "igraph", "readr", "stats"))

🛠️ Usage

1️⃣ Clone this repository:

git clone https://github.com/claudiacastrillon/PPI_Network_Analysis.git

2️⃣ Navigate to the project folder:

cd PPI_Network_Analysis

3️⃣ Run the R script:

source("PPI_network.r")

📊 Output

PPI Network Visualization with nodes and edges.

Clustering results with color-coded groups.

.CSV File Analysis showing extracted protein interactions.

🤝 Contributions

Feel free to contribute by submitting pull requests or reporting issues!

📜 License

This project is open-source. See LICENSE for details.

📩 Contact

For inquiries, contact claudiacastrillon via GitHub. 💡


