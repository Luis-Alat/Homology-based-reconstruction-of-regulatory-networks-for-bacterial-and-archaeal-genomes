# %%
import scipy.stats as stats
import os
import pandas as pd
import numpy as np
from time import time
import argparse
from ast import literal_eval

# %%
def GetExpectedProbability(counts:pd.DataFrame, ax:int) -> pd.Series:
    
    '''
    This functions computes the expected values shown as probabilities given a contingency table
    
    Parameters
    ----------
    
    counts: pd.DataFrame
        Pandas DataFrame containing the counts of each group
    
    ax: int
        Integer 0 or 1 indicating along which axis will be used to compute the expected value
    
    Return
    ------
    
    pd.Series: pandas Series indicating the probabilities for every group
    
    '''
    
    total_sums = counts.sum(axis=ax)
    
    return total_sums / total_sums.sum()

# %%
def BinaryMetrics(A:set,B:set) -> pd.DataFrame:

    '''
    
    This function computes the True positives (TP), False Negatives (FN) and False Positives (FP)
    from the comparation between the interactions found in two networks
    
    Parameters
    ----------

    A: set
        Reference set to calculate True Positives (TP) and False Negatives (FN) 
    
    B: set
        Complement set to calculate True Positives (TP) and False Positives (FP)
        
    Returns
    -------
    
        pd.DataFrame: A data frame containing the TP, FP, FN of the networks compared

    '''

    TP = len(A.intersection(B))
    FP = len(B.difference(A))
    FN = len(A.difference(B))

    df = pd.DataFrame([TP,FP,FN]).rename({0:"Counts"}, axis=1)
    df.index = ["TP", "FP", "FN"]
    
    return df

# %%
def GetContingencyTable(reference_to_random:pd.DataFrame, reference_to_compare, n:int=10000) -> pd.DataFrame:
    
    '''
    
    This function get the contingency table (TP, FP, FN) from a random network obtained
    
    Parameters
    ----------
    
    reference_to_random: pandas DataFrame
        Network containing columns "TF" (Transcription factors) and TG (Target gene). Every pair of values 
        TF-TG will describe an interaction and this frame will be shuffled
    
    reference_to_compare: set
        Set containing pairs TF-TG for each index. This set is not modified and used to compared againts the shuffled frame

    n: int; default 10000
        Number of random networks to retrieve (experiments)
        
    Returns
    -------
    
    pd.DataFrame: A pandas data frame containing the counts of TP, FP, FN for every experiment
    
    '''

    print(f"Getting contingency table with n={n}")
    
    random_contingency_table = pd.DataFrame(dtype=int, index=["TP","FP","FN"])
    number_of_interactions = reference_to_random.shape[0]

    start = time()

    for experiment in range(n):

        print(f"\tExperiment: {experiment + 1}", end="\r")
        experiment = "rand_" + str(experiment)
        random_df = reference_to_random.copy()

        # If we do not transform into a list, then the data frame matches the rows by index
        # and the shuffle effect is lost 
    
        random_df["TF"] = random_df["TF"].sample(number_of_interactions).to_list()

        Interactions_sample = set(random_df[["TF","TG"]].apply(" | ".join, axis=1))
        metrics = BinaryMetrics(reference_to_compare, Interactions_sample)["Counts"]
        random_contingency_table = pd.concat([random_contingency_table, metrics.rename(experiment)], axis=1)

    end = time()
    print(f"\tExperiment: {n}")
    print(f"\tThis functions took {end - start} secs to run")
    
    return random_contingency_table

# %%
def FilterNetworks(graph:pd.DataFrame) -> "(pd.DataFrame, pd.DataFrame)":
    
    '''
    
    This function filters the "Known" and "New" (a.k.a inferred) interactions between two networks given
    a specific format
    
    Parameters
    
    graph: pandas dataFrame
        A pandas dataframe containing the whole network
        
    Returns
    -------
    
    (known = pd.DataFrame, new = pd.DataFrame): A two dimention tuple containing the known interactions
        and the new interactions filtered
    
    '''
    
    print("Filtering networks")
    
    know_net = graph.query("STATUS == 'Known'")
    mask_new_cond_1 = graph["STATUS"] == "New"
    mask_new_cond_2 = (graph["STATUS"] == "Known") & ((graph["ORG_REF"] != "NOT_REFR_ORG") | (graph["TU_UNIT"] != "NOT_TU_REFER")) 

    new_net = graph[mask_new_cond_1 | mask_new_cond_2]

    return (know_net, new_net)

# %%
def ListDirectory(path:str, filter:bool=False) -> "list":
    
    '''
    This function list the file names found inside a directory
    
    Parameters
    ----------
    
    path: str
        String describing the path of the target directory
    
    filter: str or RegEx, default: False
        String or RegEx instruction matching with the target file names
    
    Returns
    -------
    
    list: list containing the file names
    
    '''
    
    print(f"Path directory: {path}")
    print(f"Filter: {filter}\n")
    
    file_names = pd.Series(os.listdir(path))
    
    if filter:
        file_names = file_names[file_names.str.contains(f"{filter}")]
    
    return list(file_names)

# %%
def LoadParserNet(file:str, **kwargs) -> "pd.DataFrame":
    
    '''
    This function load into a pandas DataFrame the specified network
    
    Parameters
    ----------
    
    file: str
        String describing the file name
        
    **kwargs:
        Aditional named keys for the pd.read_csv object
    
    
    Returns
    -------
    
    pd.DataFrame: Pandas DataFrame with the network loaded
        
    '''
    
    print(f"Loading file: {file}")
    
    return pd.read_csv(file, **kwargs)

# %%
def Arguments() -> argparse.ArgumentParser:

    '''
    This function processes the arguments supplied from the command line

    Returns
    -------

    argparse.ArgumentParser: argparse object containing the arguments parsed

    '''

    argp = argparse.ArgumentParser()
    
    argp.add_argument("--input_path", dest="input_path",
                      help="Path to search all the networks", type=str)
    argp.add_argument("--filter_files", dest="filter_files",
                      help="RegEx to filter the files found in --input_path", type=str)
    argp.add_argument("--output_file", dest="output_file",
                      help="Path and fiel name to place results", type=str)
    argp.add_argument("--arguments_load_files", dest="arguments_load_files",
                      help="String containing a dictionary-like syntax with the arguments to pass to the pandas read_csv method")

    args = argp.parse_args()
    args.arguments_load_files = literal_eval(args.arguments_load_files)

    return args

# %%
def ShowArguments(arguments:argparse.ArgumentParser):

    '''
    This function shows on screen the arguments passed by command line
    '''

    print("\n##### PARAMETERS #####\n")
    
    for arg in vars(arguments):
        print(f"{arg}: {getattr(arguments, arg)}")
    
    print("\n######################\n")
    

# %%
if __name__ == "__main__":

    args = Arguments()
    ShowArguments(args)

    # Listing and saving file names of the networks
    network_file_names = ListDirectory(args.input_path, filter=args.filter_files)

    # Creating empty file
    open(args.output_file, mode="w").close()

    for net in network_file_names:
    
        real_path = os.path.join(args.input_path, net)
        plain_graph = LoadParserNet(real_path, **args.arguments_load_files)
    
        # Getting the inferred interactions and the known interactions
        New_network, Known_network = FilterNetworks(plain_graph)
    
        # Getting the contingency table of the observed values of the inferred network and the network used
        # as a reference
        Interactions_new = set((New_network[["TF","TG"]].apply(" | ".join, axis=1)))
        Interactions_known = set((Known_network[["TF","TG"]].apply(" | ".join, axis=1)))
    
        Known_versus_inferred_metrics = BinaryMetrics(Interactions_known, Interactions_new)
    
        # Getting contingency table and expected values of random networks based on the inferred network
        contingency = GetContingencyTable(New_network, Interactions_known, 10000)
        probabilities = GetExpectedProbability(contingency, 1)
    
        # Doing the test G for goodness of fit
    
        Observed = Known_versus_inferred_metrics.values.reshape(-1,)
        Expected = (probabilities * Observed.sum()).values.reshape(-1,)
        
        # G-test directly
        G_d, p_value_G_d = stats.power_divergence(Observed, Expected, lambda_ = 0)
        
        print("Tests:")
        print(f"\tGd:{G_d}; p:{p_value_G_d}")
        
        print("")

        # Creating output to save
        probabilities.name = "Expected_probability"
        observed_prob = Known_versus_inferred_metrics / Known_versus_inferred_metrics.sum()
        observed_prob.columns = ["Observed_probability"]
        
        contingency_frame = Known_versus_inferred_metrics.merge(observed_prob, left_index=True, right_index=True)
        contingency_frame = contingency_frame.merge(probabilities, left_index=True, right_index=True)

        with open(args.output_file, mode="a") as oFile:
            oFile.write("\n" + real_path + "\n\n")
            oFile.write("G-test" + "\n")
            oFile.write("score: " + str(G_d) + "\n")
            oFile.write("P-value: " + str(p_value_G_d) + "\n")
            oFile.write("df: " + str(int(contingency_frame.shape[0] - 1)) + "\n\n")
            oFile.write("Matrix\n")

        contingency_frame.to_csv(args.output_file, index=True, header=True, mode="a")


