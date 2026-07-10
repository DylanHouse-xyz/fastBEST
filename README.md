# fastBEST
fastBE supplementary tools (fastBEST) is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow that bundles the generic workflow from variant calling with Mutect2, to clonal denvolution with [fastBE](https://github.com/raphael-group/fastBE). A clonal tree was inferred with fasBE with it's search method, and the variants were clustered with fastBE's cluster method, and the optimal number of clones were inferred with the [kneedle algorithim](https://kneed.readthedocs.io/en/stable/).

------------

# Overview

1. [Setup](#setup)
2. [Usage](#usage)
3. [Outputs](#workflow-output)
4. [Workflow DAG](#workflow-DAG)
5. [Citation](#citation)

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
