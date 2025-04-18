image:: add images to display the pipeline

EDNiX (EvoDevo NeuroImaging Explorer)
=======

EDNiX is a comprehensive brain MRI image processing pipeline, designed for flexibility and adaptability to various mammalian species. 
The pipeline supports anatomical MRI, functional MRI (fMRI), and PET scans images, facilitating data processing from raw images to statistical analysis. 
It accommodates as much as possible to the various constraints of different species, imaging setups, and acquisition sequences. 
EDNiX is developed to enhance cross-species and developmental neuroimaging research and has been tested on a diverse set of mammalian species from various laboratory thanks to open sources datasets. 
Along with EDNiX we also release two MRI datasets acquired in awake macaques and anesthetized mouse lemurs.
Signal extration can rely on a library of cross-species atlases (https://www.pnas.org/doi/abs/10.1073/pnas.2202491119), provided along with the pipeline, but study specific atlases may also be used.


Download 
=======

    git clone https://github.com/garincle/EDNiX

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

    cd path/to/your/EDNiX
    python install_packages.py


Singularity
--------------
    code to install and download singularity files
 
Explanations how to do it
 
Roadmap
=======

# Key Features

- Adaptable to multiple mammalian species, including humans and laboratory animals.
- Supports anatomical MRI, fMRI, and PET scans.
- Facilitates longitudinal co-coregistration for repeated individual scans.
- Includes cross-species atlas libraries for signal extraction as well as cross-species comparisons
- Tested extensively on diverse imaging setups and species.
- Encourages community feedback for ongoing improvements.
- 
# Setup Instructions

## 1.	Clone or download the EDNiX Repository

Clone or download the EDNiX repository from the official source.

`git clone https://github.com/garincle/EDNiX`

## 2.	Install Python Dependencies

Ensure Python 3.9 is installed on your system. Run the install_packages.py script to install necessary dependencies:
`conda create -n EDNiX python=3.9`
`conda activate EDNiX
`
`cd path/to/your/EDNiX`

`python install_packages.py `

## 3.	Explore the Launcher_anat/ Launcher_func/ folder

Review and copy past in a new text sheet the launcher the most adapted to your species and imaging modalities (e.g., anatomical MRI, fMRI)

## 4.	Review and adapt your new Launcher_ to your modality/species 

Change the paths reading the singularity folders, altas folders, as well as the various variables when necessary

### Usage Guide

To process MRI data, use the launcher scripts provided in the Launcher_**/ folder. These scripts are tailored for various species and data types (e.g., anatomical MRI, fMRI). For Usage:
- Use Launcher_anat_Human.py for processing human anatomical MRI data.
- Use Launcher_func_Mouse.py for processing mouse fMRI data.
Each script contains parameters and settings that you can customize to match your dataset.

### Folder Structure
1.	Launcher_**/: Contains example scripts for launching modality and species-specific data processing pipelines as well as statistical analysis (Launcher_analysis/).
Those launchers, allows to start functions contained within in one of the following directories:
a.	anatomical/: Preprocessing scripts for anatomical MRI, including orientation correction, bias field correction, skull stripping and co-registration and more.
b.	functions/: Preprocessing scripts for functional MRI (fMRI), such as orientation correction, motion correction, artefact correction, temporal filtering, masking, co-registration and more.
c.	PET/: Preprocessing scripts for functional MRI (fMRI), such as orientation correction, motion correction, artefact correction, temporal filtering, masking, co-registration and more.
d.	analyses/: Includes statistical analysis scripts for MRI data.
