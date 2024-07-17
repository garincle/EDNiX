image:: add images to display the pipeline

EasyMRIbrain
=======

EasyMRIbrain is a brain MRI images treatment pipeline, coded in python, that provides important flexiblility and adaptation to different mammlian species (anatomical, fMRI, and more soon!), PET, from raw data to statstical analysis. 
EasyMRIbrain has been design to adapt to the various constrain of each species as well as to the different MRI providers, sequences, softwares used in the various laboratories.
This pipeline has been tested on the largest sample of species possible so far (XXpaper). Any issues encounter during processing can be reported in the issues section to improve the flexibility of this pipeline.
EasyMRIbrain enables longitudinal co-coregistration (when a individual is scanned multiple time).
Signal extration can rely on a library of cross-species atlases (https://www.pnas.org/doi/abs/10.1073/pnas.2202491119), provided along with the pipeline, but specific altases may also be used.


Download 
=======

    git clone https://github.com/garincle/EasyMRIbrain_sing

Dependencies
=======

Python and conda
--------------

**1. Virtual environment**

We recommend that you install and use ``conda`` (https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
Then, to create and activate a new python environment:

    conda create -n emb python=3.9
    conda activate emb
    
**2. Python dependencies**

In this enviroment you will need to install the following dependencies

    cd path/to/your/EasyMRIbrain_sing
    python install_packages.py


Singularity
--------------
    code to install and download singularity files
 
Explanations how to do it
 
Roadmap
=======

**1. From raw to BIDS**

**2. Defining variables adapted to you objectives and you species**

**3. Checking the quality**

**4. Statistical analysis**

