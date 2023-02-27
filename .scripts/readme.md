# Notes

This short markdown has the purpose to describe some "particularities" of some scripts.

## General comments

The following scripts already have default file outputs and, when executed, __<span style="color:#D55D44">will NOT check for existing files</span>__ (at least not most of them). The same goes for input files, those will be located in a specific folder and expected to be found there (e.g. fasta files in _genomes/_ and modified networks in _modifiedNetworks/_). Those are just a few considerations to take into account.

On the other hand, the following files work through a series of steps found as bash functions, so if for some reason you wish to omit a step (e.g., no longer run __proteinortho__ because you already have the files) the corresponding function can be commented "#"

## Running main files

<br>

### pipeline_V4 file

<br>

_pipeline_V4.sh_ is the script that goes from preprocessing the networks of the six organism models, to extending and analyzing them. From the top directory, it can be executed as

```{bash}
bash scripts/pipeline_V4.sh
```

_pipeline_V4.sh_ calls and executes multiple scripts found here using python, R or perl one-liners commands (an enviroment.yml file was supplied in order to run this properly)

<br>

### predict_bact_arch directory

<br>

Inside of this folder, we can find _predict_network_without_reference.sh_. This the script developed to run __proteinortho__ and predict networks based on orthology relationships between model organisms and "target" organisms (~12000 bacteria and ~650). When processing a large number of files/genomes, __<span style="color:#D55D44">this script will run</span>__ several processes __<span style="color:#D55D44">in parallel (background) and use the available cores</span>__ to speed up this process, so your computer may be a bit unavailable. Additionally, this processing will take a bit of time (a long time to be honest).

_predict_network_without_reference.sh_ works as a tool type whose parameters can be viewed with either one of the following two lines

```{bash}
bash predict_network_without_reference.sh
bash predict_network_without_reference.sh -h
```

__<span style="color:green">As a reminder, _predice_network_without_reference.sh_ assumes some file names and paths to run</span>__, for example, it assumes _Escherichia coli_ with a nickname of _Ecoli_. This in order to identify and keep track of which files and at what time each of them are being processed, to change this, this must be modified directly in the file code (see __### Defining default variables__, __### Names of model organism__ and __OPTIONAL arguments__
sections inside the script).

Beyond this, this file was not run directly (although it could be) but was called by another script found as _pipeline_bact_arch.sh_,which, in turn, gets the metrics of the network topology

From the parent folder, that file can be run as easy as...

```{bash}
bash scripts/pipeline_bact_arch.sh
```
