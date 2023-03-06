import sys
import networkx as nx

Network = sys.argv[1]
Columns = sys.argv[2]

ColumnLis = Columns.strip().replace(" ","").split(",")
ColumnLis = [int(num) for num in ColumnLis]
ColumnLis[0] = int(ColumnLis[0]) - 1

MatrixNet = []
Nodes = []

with  open(Network) as iFile:
	for line in iFile:
		line = line.strip().split("\t")
		MatrixNet.append(line[ColumnLis[0]:ColumnLis[1]])
		Nodes.append(line[ColumnLis[0]])
		Nodes.append(line[ColumnLis[1] - 1])

Nodes = list(set(Nodes))

Graph = nx.DiGraph()
Graph.add_nodes_from(Nodes)
Graph.add_edges_from(MatrixNet)

MatrixNet = None
Nodes = None

hubs , authorities  = nx.hits(Graph)

hubs = [[key,value] for key,value in hubs.items()]
authorities = [[key,value] for key,value in authorities.items()]

[print(h[0],"\t",h[1],"\t",a[1]) for h,a in zip(hubs,authorities)]
