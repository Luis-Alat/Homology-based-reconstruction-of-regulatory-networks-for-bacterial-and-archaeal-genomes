# Homology based reconstruction of regulatory networks on bacteria and archaeal genomes

In this work, we have inferred the gene regulatory network of 12230 bacterial and 649 archaeal genomes, by considering as reference six organisms with well-known regulatory interactions. Here, the repository of results/programs is presented, as well as a brief description of the folders associated with the project.

More information is avilable at [doi.org/10.3389](https://doi.org/10.3389/fmicb.2022.923105)

The full project can be found on [drive](https://drive.google.com/drive/folders/1n8zZq0-HJg0ULoMiqSr_IhApmJSTgn4h?usp=sharing)

UPDATE:

* Mapping from TUs it's already based on presence/absent of TG orthologues rather than the first match
* _B.subtillis_ network was replaced by the network provided by SubtiWiki

## Getting started

A large part of this repository (results) can be obtained by executing some files as shown below and those are already available here but, because of the limited size to storage all the files created in this project (mainly for the bacteria and archaeal genomes), most of them are not avilable here, but can be downloaded _**[here](insertart drive)**_

### Prerequisites

This work was built using:

* [python](python.org) - Programming language
* [perl](https://www.perl.org/) - Programming language
* [R](https://www.r-project.org/) - Programming language
* [Bash Unix shell](https://www.gnu.org/software/bash/) - Unix shell and command language
* [CytoScape](https://cytoscape.org/) - Visualizer and network analysis
* [proteinortho](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-124) - Tool to detect orthologous genes

Additionally, various libraries were used such as [pandas](https://pandas.pydata.org/) or [networkx](https://networkx.org/documentation/stable/#) for python or [Rcy3](https://bioconductor.org/packages/devel/bioc/vignettes/RCy3/inst/doc/Overview-of-RCy3.html) for R. But don't worry!, most of the requirements shown above can be installed using the [conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment provided here (enviroment.yml). To do this, the following line of code can be executed:

```
conda env create -f environment.yml
```

And activated as: 

```
conda activate homology_project
```

If you wish to know more about managing environments with [conda](https://docs.conda.io/projects/conda/en/latest/index.html) you can do it [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)

On the other hand, you can download [CytoScape](https://cytoscape.org/) by clicking on the respective link and following the instructions and that's it!

### Running

Once activated the conda enviroment, __<span style="color:green">you will need to launch CytoScape manually</span>__ and to be placed in the **scripts** folder. 

In order to recreate the exact and most of the folders presented here, you can run:

```
bash main.sh -g Fasta_files_path.txt -n Nets_files_path.txt -l Labels_organism.txt -t ../tus/moreno_models --extended_nets_output ../network/predicted_nets/models/ --proteinortho_output ../proteinortho/models/ --tables_output ../analysis/models/tables --cytoscape_output ../analysis/models/cytoscape --coreg_output ../analysis/models/coreg/ --networkx_output ../analysis/models/hits/ --g_test_output ../analysis/models/gtest/ --literature_output ../analysis/models/pubmed/ --literature_input Ecoli_Rules_PubMed.txt
```

I know, those are a lot of arguments, but allows more flexibility in the code. Here, this command will preprocess the networks of the reference organisms until it analyzes them and creates their extended networks based on the orthology relationships between species. While if we execute

For a deeper explanation about the arguments of the you can run one of the following commands

```
bash main.sh
```

or

```
bash main.sh --help
```

On the other hand, *main_BacArc.sh* will do something similar, but this time for the 12000 bacterial genomes and 650 archaea approximately.

To run the script with the bacteria, you can run the following:

```
bash main_BacArc.sh -g Fasta_files_path.txt -t Targets_bacteria.txt -n Nets_files_path.txt -l Labels_organism.txt -u Tus_bacteria.txt --proteinortho_output ../proteinortho/bacteria --extended_nets_output ../network/predicted_nets/bacteria
```

In order to execute the above script, multiprocessing is available and, therefore, multiple processes can be launched at the same time in batches. Therefore, an argument called _--batches_ can be set, by default, it is not multiprocessing or _--bactches 1_ .

For a device with 8GB RAM and 8 cores, I recommend _--batches 90_ . Although be aware that this could leave your device unavailable for a few hours.

Finally, some steps of the pipeline can be skip if you want but, for that, you'll need to comment the respective line of code of the scripts. For example, if you want to avoid running proteinortho, you'll need to comment the line calling the implementation of proteinortho in the script (*RunProteinortho*), If you do that, just take into account some steps of the pipeline recieves as input the output of the steps before. 

## Folder description

A lot of different folders were created to organize this work, a description of the content of the main one are shown below:

* **analysis**: Analysis for the 6 model organisms (pubmed, hubs, cytoscape, etc) and bacteria or archaeal genomes (networkx) 

* **genomes**: Genomes used in this work in fasta format for models, bacterial and archaeal

* **network**: Raw networks and modified networks of the six organism models. In addition, networks predicted in the model, bacteria and archaeal organisms

* **others**: Others files

* **proteinortho**: Proteinortho outputs for model, bacteria and archaeal organisms

* **scripts**: python, perl, bash and R scripts to recreate this project

* **tus**: Transcription units for model, bacteria and archaeal organisms
