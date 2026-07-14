# fastBEST
fastBE supplementary tools (fastBEST) is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow that bundles the generic workflow from variant calling with Mutect2, to clonal denvolution with [fastBE](https://github.com/raphael-group/fastBE). A clonal tree was inferred with fasBE with it's search method, and the variants were clustered with fastBE's cluster method, and the optimal number of clones were inferred with the [kneedle algorithim](https://kneed.readthedocs.io/en/stable/).

------------

# Overview

1. [Setup](#setup)
2. [Overview](#overview)
3. [Usage](#usage)
4. [Outputs](#workflow-output)
5. [Workflow DAG](#workflow-DAG)
6. [Citation](#citation)

--------------

# Dependencies

* [conda](https://github.com/conda-forge/miniforge), version >24.1.2
* [Snakemake](https://snakemake.readthedocs.io/en/stable/), version >=7.32.4

---------------

# Setup

This pipeline would require that both [conda](https://github.com/conda-forge/miniforge) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) be installed; ensure bioconda and conda-forge channels are added. Below are the steps to ensure those requirements are met.

1. Install `conda` through [miniforge]((https://github.com/conda-forge/miniforge#install)).
2. Ensure the appropiate conda channels are added;


```
conda config --add channels bioconda
conda config --add channels conda-forge
```

1. Install appropiate [Snakemake](https://snakemake.readthedocs.io/en/stable/):

```
conda create -c conda-forge -c bioconda --name snakemake snakemake
```

---------------

# Overview

Our Snakemake pipeline's main function is to automate variant calling with Mutect2, phylogenetic reconstruction with fastBE, and subsequent analysis with custom-made python scripts to generate cancer cell fraction (CCF) stacked bar chart of the clonal composition of a sample, a phylogenetic reconstructed tree from fastBE outputs, and a fishplot to visualize the evolutionary dynamics of a cohort overtime.

One major advantage Snakemake worklows, along with other workflow management systems is the utilization of a configuration file to generalize a pipeline. We split the workflow to use two configuration files: one config file for samples. This is in the format;

1. `Samples` - the samples dictionary where you will place the path to your input BAM files.
2. `sample_name` - The name of the sample you are analyzing. This will be used to label newly created folders where your outputs will be placed.
3. `sample rows` - First row contains path to `tumor bam files`. Second row holds the `name of normal sample`. Final row is the path to `normal bam file`.

> [!IMPORTANT]
> This workflow is intended for use with tumor-normal samples. If you wish to use a panel of normals, or tumor only: modify rules/mutect2.smk to alter the variant calling shell command.



---------------



# Citation

For fastbe, please cite the original authors:

```
@article{schmidt2024regression,
  title={A regression based approach to phylogenetic reconstruction from multi-sample bulk DNA sequencing of tumors},
  author={Schmidt, Henri and Raphael, Benjamin J},
  journal={PLOS Computational Biology},
  volume={20},
  number={12},
  pages={e1012631},
  year={2024},
  publisher={Public Library of Science San Francisco, CA USA}
}
```
