# %%
import pandas as pd
import networkx as nx
import argparse
from ast import literal_eval

# %%
def Arguments() -> argparse.ArgumentParser:

    '''
    This functions parses the command line arguments supplied by the user
    '''

    argp = argparse.ArgumentParser()
    
    argp.add_argument("--input_network", dest="input_network", help="Input file to calculate hubs of the network", type=str)
    argp.add_argument("--arguments_read_net", dest="arguments_read_net", help="Arguments to be used while reading the network by the pandas method read_csv")
    argp.add_argument("--output_path", dest="output_path", help="Path and name to place the files", type=str)

    args = argp.parse_args()

    args.arguments_read_net = literal_eval(args.arguments_read_net)

    return args

def ShowArguments(arguments:argparse.ArgumentParser):

    print("\n##### PARAMETERS #####\n")
    
    for arg in vars(arguments):
        print(f"{arg}: {getattr(arguments, arg)}")
    
    print("\n######################\n")

def CreateScoresFrame(nodes:list, hub_scores:list, authorities_scores:list):

    df = pd.DataFrame(columns=["Node", "Hub_Score"," Aut_Score"])
    df["Node"] = nodes
    df["Hub_Score"] = hub_scores
    df["Aut_Score"] = authorities_scores
    df = df.sort_values(by="Hub_Score", ascending=False)

    return df

# %%
if __name__ == "__main__":

    args = Arguments()
    ShowArguments(args)

    network_df = pd.read_csv(args.input_network, **args.arguments_read_net)
    network_df = network_df.to_numpy().tolist()

    graph = nx.DiGraph(network_df)
    hubs_values, authorities_values = nx.hits(graph)

    scores_net_df = CreateScoresFrame(nodes=hubs_values.keys(), 
                                      hub_scores=hubs_values.values(), 
                                      authorities_scores=authorities_values.values())
    
    scores_net_df.to_csv(args.output_path, sep="\t", header=True, index=False)


