#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
import networkx as nx
import argparse as arg
import sys

if __name__ == "__main__":

    argp = arg.ArgumentParser(description="Script to get general topological metrics of networks")
    argp.add_argument("--inputFile", dest="inputFile", help="Network input file to read")
    argp.add_argument("--outputFile",dest="outputFile", help="Output file path and name to save results")
    args = argp.parse_args()

    print(f"Input file: {args.inputFile}")
    print(f"Output file: {args.outputFile} \n")

    print("Reading file...")
    
    try:
        # Reading only TF (0) y TG (1)
        network = pd.read_csv(args.inputFile,sep="\t",header=None)[[0,1]]
        network = network.rename({0:"TF",1:"TG"}, axis=1)
    except Exception as err:
        sys.stderr.write(str(err))
        exit(1)

    print("Creating network...")
    g = nx.DiGraph()
    g.add_edges_from(network.to_numpy())



    print("Creating output file...")
    with open(args.outputFile, "w") as oFile:
        # Number of nodes
        oFile.write("#Number of nodes: " + str(len(g.nodes)) + "\n")
        # Number of edges
        oFile.write("#Number of edges: " + str(len(g.edges)) + "\n")
        # Clustering coefficient
        oFile.write("#Clustering coefficient: " + str(nx.cluster.average_clustering(g)) + "\n")
        # Network density
        oFile.write("#Density: " + str(nx.density(g)) + "\n")
        # Connected components
        oFile.write("#Weakly connected components: " + str(nx.algorithms.number_weakly_connected_components(g)) + "\n")
        # Number of self loops
        oFile.write("#Self loops: " + str(nx.number_of_selfloops(g)) + "\n")


    print("Getting metrics...")
    # Getting metrics by nodes

    clustering_c = pd.Series(nx.algorithms.cluster.clustering(g))
    clustering_c = clustering_c.rename("Clustering_coefficient")

    closeness = pd.Series(nx.algorithms.closeness_centrality(g))
    closeness = closeness.rename("Closeness_centrality")

    edge_count = {}
    for i in g:
        edge_count[i] = len(g.edges(i))

    edge_count = pd.Series(edge_count)
    edge_count = edge_count.rename("Edge_count")

    in_degree = pd.Series(dict(g.in_degree()))
    in_degree = in_degree.rename("In_degree")

    out_degree = pd.Series(dict(g.out_degree()))
    out_degree = out_degree.rename("Out_degree")

    betweenness = pd.Series(nx.algorithms.betweenness_centrality(g))
    betweenness = betweenness.rename("Betweenness_centrality")

    # Defining dataFrame to save (Mergin/saving to file)
    print("Saving final file...")
    Statistics = pd.DataFrame(index=g.nodes)
    Statistics = Statistics.join([clustering_c,closeness,edge_count, in_degree, out_degree,betweenness])

    Statistics.to_csv(args.outputFile, sep="\t", header=True, index=True,index_label="#Nodes", mode="a")

    print("All done\n\n")



