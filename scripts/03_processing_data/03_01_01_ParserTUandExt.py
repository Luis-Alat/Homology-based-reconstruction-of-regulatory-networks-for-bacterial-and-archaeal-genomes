# %%
import pandas as pd
import argparse
import os
import ast
from time import time

# %%
def ParserArguments() -> argparse.ArgumentParser:
    
    '''
    This functions parser the command line arguments supplied by the user
    
    '''

    argp = argparse.ArgumentParser(description="Cross modified network and transcription unit (TU) files")

    argp.add_argument("--input_net", dest="input_net", help="Input file of the extended network", type=str)
    argp.add_argument("--arguments_read_net", dest="arguments_read_net", help="Arguments to be passed to the pandas read csv method while reading --inputNet")
    argp.add_argument("--input_tus", dest="input_tus", help="Input file of transcriptional units", type=str)
    argp.add_argument("--arguments_read_tus", dest="arguments_read_tus", help="Arguments to be passed to the pandas read csv method while reading --inputTus")
    argp.add_argument("--use_first_match", dest="use_first_match", help="True or False", type=bool, default=False)
    argp.add_argument("--output_path_name", dest="output_path_name", help="Path to save files", type=str)

    args = argp.parse_args()

    # Parsing arguments for pandas
    args.arguments_read_net = ast.literal_eval(args.arguments_read_net)
    args.arguments_read_tus = ast.literal_eval(args.arguments_read_tus)

    return args

def ShowArguments(arguments:argparse.ArgumentParser):

    print("\n##### PARAMETERS #####")
    
    for arg in vars(arguments):
        print(f"{arg}: {getattr(arguments, arg)}")
    
    print("######################","\n")

# %%
if __name__ == "__main__":

    start = time()

    args = ParserArguments()
    ShowArguments(args)
    
    network_df = pd.read_csv(args.input_net, **args.arguments_read_net)
    network_df.rename({args.arguments_read_net["usecols"][0]: "TF",
                       args.arguments_read_net["usecols"][1]: "TG"}, axis=1, inplace=True)
    network_df = network_df.astype(str)
    
    tus_df = pd.read_csv(args.input_tus, **args.arguments_read_tus)
    tus_df.rename({0:"id", 1:"TU", 2:"strand"}, axis=1, inplace=True)
    tus_df = tus_df.astype(str)
    
    print("Cross network-tus info...")
    
    final_output = []
    for TF, TG in zip(network_df["TF"], network_df["TG"]):
        for index_tus, TU in enumerate(tus_df["TU"]):
            
            if args.use_first_match:
            
                tus_as_list = TU.replace(" ","").split(",")
                if tus_as_list[0].find(TG) > -1:
                
                    # Iterating over the genes inside the TU
                    for gene in tus_as_list:
                    
                        # Giving a format to the output: STRAIN-TF-TG-TU_id
                        final_output.append( [tus_df.iloc[index_tus]["strand"], TF, gene, tus_df.iloc[index_tus]["id"] + "_tu"] )
            else:
            
                tus_as_str = TU.replace(" ", "")
                tus_as_list = TU.replace(" ","").split(",")
    
                if tus_as_str.find(TG) > -1:
                
                    # Iterating over the genes inside the TU
                    for gene in tus_as_list:
                    
                        # Giving a format to the output: STRAIN-TF-TG-TU_id
                        final_output.append( [tus_df.iloc[index_tus]["strand"], TF, gene, tus_df.iloc[index_tus]["id"] + "_tu"] )
    
    final_output = pd.DataFrame(final_output).drop_duplicates()

    print(f"Making {args.output_path_name}")
    final_output.to_csv(args.output_path_name, header=False, index=False, sep="\t")
    print(f"All donde in {time() - start}")


