# UCA Estimation

**If you use this method in published research, please cite the paper below.**
> [Jensen C, Kleinstein S, Hoehn K (2026). UCA PAPER. bioRxiv https://www.biorxiv.org/content/link](link)

Dowser can estimate the unmutated common ancestor (UCA) or full length naive B cell sequence of a B cell lineage. To do so, a few files and executables will need to supplied along with having Python and a few packages installed. 

## Preparing for UCA Estimation

### Obtain the OLGA Model Folder
We use `OLGA` for part of our UCA estimation. Please clone their github (you will need this later).

``` bash
git clone https://github.com/statbiophys/OLGA.git
```


### Have Python and some Packages Installed
Before proceeding, verify that Python is installed on your system.

Open a terminal (Mac/Linux) or Command Prompt (Windows) and run:

``` bash
python --version
# or, on some systems:
python3 --version
```

If Python is installed, you'll see something like `Python 3.11.4`. If you get an error or "command not found", you'll need to install it.

To install Python, go to [https://www.python.org/downloads/](https://www.python.org/downloads/) and download the latest stable version for your operating system (Windows, macOS, or Linux).

Once you have Python installed you can use the following function to install the packages needed for UCA estimation. 

``` r
py_path <- Sys.which("python3") # it may be python
dowser::installPythonDependencies(python = py_path)
```

### Install IgPhyML

IgPhyML is required for UCA estimation. If you already have a copy installed before January of 2026, you will need to update it. Follow [these instructions](https://igphyml.readthedocs.io/en/latest/install.html#compiling-from-source-recommended-for-linux) on how to install IgPhyML.

## Data Pre-processing 

Before this can occur, B cell sequences must be separated into clonal groups, had their clonal germline sequences reconstructed, and in the case of paired chain estimation, had the light subgroup identified. At this point we assume that clonal identification has been completed. This example will show how to find the UCA for a paired lineage, although it can be run on just one chain.


``` r
library(dowser)

# load example data 
data("ExampleMixedDb")

# if you have paired data and want the UCA for both you'll need to do resolveLightChains
# if only want heavy chains or have bulk data you can skip this step
ExampleMixedDb <- resolveLightChains(ExampleMixedDb)

# create the germline that is used during tree reconstruction 
# if you plan on only using heavy chains you can leave clone as clone_id
references <- readIMGT("imgt/human/vdj")
ExampleMixedDb <- createGermlines(ExampleMixedDb, references, trim_lengths = TRUE, 
                                  clone = 'clone_subgroup_id')

# format the data in to a formatClones object
# if you're only using heavy chains set chain = "H"
clones <- formatClones(ExampleMixedDb, chain = "HL")
```


## Estimating the UCA

You should now have a clone object with 4 clones. To do reconstruct the trees and UCAs, use the `getTreesAndUCAs` function in R. The base inputs of this function are:
- clones = the clones object
- data = your input airr tsv (ExampleMixedDb here)
- references = references read in using readIMGT
- exec = the location of the IgPhyML executable
- model_folder = the folder that contains Olga's heavy chain model parameters

Additionally, for paired estimation you'll want to include:
- model_folder_igk = the folder that contains Olga's light chain kappa model parameters  -- not required for heavy chain only inference
- model_folder_igl = the folder that contains Olga's light chain lambda model parameters -- not required for heavy chain only inference
- chain = "HL", tells the function to estimate both chains 
- partition = "hl", tells IgPhyML to partition tree reconstruction by chain


``` r
clones_UCA <- getTreesAndUCAs(clones = clones, data = ExampleMixedDb,
                               references = references,
                               exec = "igphyml/src/igphyml",
                               model_folder = "OLGA/olga/default_models/human_B_heavy",
                               model_folder_igk = "OLGA/olga/default_models/human_B_kappa",
                               model_folder_igl = "OLGA/olga/default_models/human_B_lambda",
                               nproc = 1, chain = "HL", partition = "hl")
```

The output of this function is your clones object with a new column `UCA` that contains the inferred UCA. This can be found in the data frame containing sequence information for each clone, or in the germline node of the tree. Additionally, you will have a `trees` column with the IgPhyML inferred lineage trees. 
