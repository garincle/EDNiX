image:: add images to display the pipeline

EasyMRIbrain
=======

EasyMRIbrain enables a flexible analysis of MRI images of different mammlian species (anatomical, fMRI, and more soon!), PET, from raw data to statstical analysis. EasyMRIbrain has been design to adapt to the various constrain of each species as well as to the different MRI, softwares used in the various laboratories.
This pipeline may not cover all issues encounter during processing, but has been tested on the largest sample of species possible so far. Any issues encounter during processing can be reported in the issues section to improve the flexibility of this pipeline.
EasyMRIbrain enables longitudinal co-coregistration (when a individual is scanned multiple time).
Signal extration can rely on a library of cross-species atlases (https://www.pnas.org/doi/abs/10.1073/pnas.2202491119), provided along with the pipeline.

Download 
=======

    git clone https://github.com/garincle/EasyMRIbrain_sing

Dependencies
=======

Python and conda
--------------

**1. Virtual environment**

We recommend that you install and use ``conda`` (https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
To create and activate a new python environment:

    conda create -n emb python=3.9
    conda activate emb
    
**2. Python dependencies**

In this enviroment you can install the following dependencies

    cd path/to/your/EasyMRIbrain_sing
    python install_packages.py


Singularity
--------------
    code to install and download singularity files
 
Explanations how to do it
 
Roadmap
=======

**1. Setup a virtual environment**
