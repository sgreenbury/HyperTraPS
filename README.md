# HyperTraPS

This software provides the computational tools to perform statistical inference using the HyperTraPS framework.

The methodology is described and used in the following papers:

* [HyperTraPS: Inferring probabilistic patterns of trait acquisition in evolutionary and disease progression pathways](https://doi.org/10.1016/j.cels.2019.10.009)\
([Cell Systems](https://doi.org/10.1016/j.cels.2019.10.009) or [arXiv](https://arxiv.org/abs/1912.00762))
* [Evolutionary Inference across Eukaryotes Identifies Specific Pressures Favoring Mitochondrial Gene Retention](https://doi.org/10.1016/j.cels.2016.01.013)
* [Precision identification and prediction of high mortality phenotypes and disease progression pathways in severe malaria without requiring longitudinal data](https://doi.org/10.1101/425132)

## Set-up

Retrieve a local copy with:
```bash
    git clone https://github.com/sgreenbury/HyperTraPS.git
```
To compile the C++ source code, a `g++` compiler with `-std=c++11` support is required. If this is the case, run
```bash
    ./configure.sh
```
which will create the necessary local folders and compile the C++ source to binaries.

A python3 environment is required with the following packages for full functionality provided in `requirements.txt`.

If you are working with anaconda, you can install these into a "HyperTraPS" environment with
```bash
    conda create --name HyperTraPS -y python=3.7
    conda activate HyperTraPS
    pip install -r requirements.txt
```

## Tutorial

More detailed information on the full usage of the package is outlined in the sections below. However, a Jupyter Notebook is also available in:
```bash
    notebooks/HyperTraPS_tutorial.ipynb
```
with an example case performing HyperTraPS inference on a sythetic cross-sectional dataset to illustrate functionality.

## Dataset format

In order to perform inference with HyperTraPS, all that is required is a comma-separated values (csv) file in the format of _N_ samples (rows) x _L_ features (columns). If the data is longitudinal, the samples from one object should have the same index and be time ordered. If the data is phylogenetic, the samples should have the same index as the leaves contained in the corresponding phylogeny.

The csv file requires a header with the name of each feature and an index column. The table below illustrates the format for a dataset described by matrix __D__:

|                   |         Feature 1 | ... |         Feature L |
| ----------------- |-------------------| ----|-------------------|
|   Sample index  1 |  _D<sub>1,1</sub>_| ... |  _D<sub>1,L</sub>_|
|               ... |               ... | ... |               ... |
|   Sample index _N_|  _D<sub>N,1</sub>_| ... |  _D<sub>N,L</sub>_|


The values taken by _D<sub>i,j</sub>_ may be 0 (feature is absent from sample) or 1 (feature is present in sample). If the dataset is cross-sectional, a value of 2 can also be present if the value is missing. Note that the inference is currently not compatible with missing data for longitudinal and phylogenetic datasets but is an active area of research.

## Generic conversion of dataset to transitions

To perform inference of your dataset, the standard procedure is to create a subfolder and copy the scripts from `scripts/` into this subfolder. This can be achieved with the following general format
```bash
    mkdir <subfolder>
    cd <subfolder>
    cp ../scripts/* ./
    ./run_convert_data_to_transitions.sh\
      <csv file location>\
      <dataset type>\
      <optional: phylogeny file location>
```
All the default scripts required for performing inference and generating output post-processing are located in the subfolder.

## Example cross-sectional, longitudinal and phylogenetic datasets

An example of a cross-sectional dataset that can be analysed is contained within `datasets/synthetic/` along with other example datasets. A specific example can be set-up with
```bash
    mkdir test_data_2/
    cd test_data_2/
    cp ../scripts/* ./
    ./run_convert_data_to_transitions.sh\
      ../datasets/synthetic/test-data-2-7.csv\
      cross-sectional
```
An example of a synthetic longitudinal dataset that can be analysed is also contained within `datasets/synthetic/`. This example can be set-up with
```bash
    mkdir test_data_2_longitudinal/
    cd test_data_2_longitudinal/
    cp ../scripts/* ./
    ./run_convert_data_to_transitions.sh\
      ../datasets/synthetic/test-data-2-7-longitudinal.csv\
      longitudinal
```
An example of a cross-sectional dataset relating to genetic alterations in ovarian cancer is contained within `datasets/ovarian-cgh/`. The example can be set-up with
```bash
    mkdir ovarian_test/
    cd ovarian_test/
    cp ../scripts/* ./
    ./run_convert_data_to_transitions.sh\
      ../datasets/ovarian-cgh/ovarian-cgh-header.csv\
      cross-sectional
```
The ovarian alteration dataset is publicly available as part of the R package [Oncotree](https://rdrr.io/rforge/Oncotree/man/ov.cgh.html).

An example of a phylogenetic dataset relating to the acquisition of drug-resistance genetic polymorphisms is contained within `datasets/tuberculosis/` and can be set-up with
```bash
    mkdir tuberculosis_test/
    cd tuberculosis_test/
    cp ../scripts/* ./
    ./run_convert_data_to_transitions.sh\
      ../datasets/tuberculosis/tuberculosis-v5-header-7.csv\
      phylogenetic\
      ../datasets/tuberculosis/ng.2878-S2.txt
```
The tuberculosis dataset and phylogeny are made available with the associated research [publication](https://doi.org/10.1038/ng.2878) of Casali _et_ _al._ 2014.

Phylogenies are required to be provided in Newick file format for the python ete3 package to process.

Once `run_convert_data_to_transitions.sh` is run, several files named `transitions.*` are output that are utilised for subsequent procedures.

## Running the MCMC sampler

Once the conversion of the raw dataset to the correct transition dataset has been performed, the MCMC sampler may be run. To see options associated with the sampler, perform
```bash
    ../bin/RUN_MCMC_SAMPLER -h
```
Or to simply run the sampler with default options the script
```
    ./run_mcmc_sampler.sh
```
Depending on the size and dimensionality of your dataset, this could take minutes or days to run. The files output from the sampler are called `stats.csv` and `forwards.txt`. `stats.csv` contains the likelihood and acceptance ratios at each iteration. `forwards.txt` contains a sample from the MCMC chain every 100 iterations. Once a "burn-in" period is considered to have been reached, samples from "forwards.txt" may be treated as independent samples from the posterior distribution.

To examine traces of the likelihood from the sampler, you can run
```bash
    ./run_plot_mc_stats.sh
```
which will output a `.png` file displaying the MCMC trace alongside acceptance statistics for a basic examination of the mixing performance.

## Simulating random walks with posterior samples

Once it has been established the number of iterations for the "burn-in" period, random walks using samples `forwards.txt` can be performed. The random walks may be complete performing transitions from the zero state to the one state or "match the data" by performing transitions between all source and target states in the original dataset. To perform this simulation, the following can be run with more information on setting options for the simulated random walks
```bash
    ../bin/RUN_PW -h
```
In particular, the burn-in period from which to choose samples from can be set with the option `-b`. A ready made script to achieve the random walk simulations is also available by running
```bash
    ./run_pw.sh
```
Once the program has run, several files with the format `forwards\_*.txt` are output that can be analyses with plotting scripts described in the next section.


## Plotting outputs

Simulated random walks can be examined with ready made plots of the following types:
* Hypercube transition graph
* Ordering histograms
* Feature graph

Each can be run with default parameters with the following scripts
```bash
    ./run_plot_hypercube_graph.sh
    ./run_plot_ordering_histograms.sh
    ./run_plot_feature_graph.sh
```
The contents of the individual scripts can be customised in many ways to give the user control over the nature of the plots.

## Performing regularisation

The samples output from the MCMC sampler make use of the full space of parameters when choosing a first- or second- order model. As introduced and described [here](https://doi.org/10.1101/409656), a regularisation approach can be utilised to obtain more parsimonious models. To run the regularisation procedure, the following program is included
```bash
    ../bin/RUN_INFORMATION_CRITERION -h
```
and a ready made script can be run with
```bash
    ./run_information_criterion.sh
```
The regularisation procedure outputs `forwards\_information_criterion*` files capturing the regularisation process which may be plotted with the below script. In addition, new samples are stored `forwards\_regularised.txt` with the top 100 (as default) samples that optimise the Akaike Information Criterion (AIC) during the search.

## Plotting regularisation

To plot the output of the regularisation procedure, the following ready made script generates a plot of how the AIC changes with respect to the number of parameters. This can be run with
```bash
    ./run_plot_information_criterion.sh
```
and subsequently generates the file `forwards-information-criterion.pdf` with the plotted output. A comparison with the regularisation procedure applied to a comparison first order model can also be achieved with this script by providing the relative path of the folder where this model has been run and the regularisation procedure also performed.

The following script
```bash
    ./run_pw_regularised.sh
```
simulates random walks with the new regularised parameterisations and can be examined with the same plotting scripts as above.
