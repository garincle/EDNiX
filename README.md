EasyMRIbrain
=======

EasyMRIbrain enables a flexible analysis of MRI images of different mammlian species (anatomical, fMRI, and more soon!), PET, from raw data to statstical analysis. EasyMRIbrain has been design to adapt to the various constrain of each species as well as to the different tools used in laboratories.
This pipeline may not cover all issues encounter during processing, but has been test on the largest sample of species possible so far. Any issues encounter during processing can be reported in the issues section to improve the flexibility of this pipeline.
EasyMRIbrain also for longitudinal co-coregistration (when a individual is scanned multiple time). Signal extration can rely on a library of cross-species atlases, provided along with the pipeline.


Install
=======


Dependencies
=======

Latest release
--------------

**1. Setup a virtual environment**

We recommend that you install ``nilearn`` in a virtual Python environment,
either managed with the standard library ``venv`` or with ``conda``
(see `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ for instance).
Either way, create and activate a new python environment.

With ``venv``:

.. code-block:: bash

    python3 -m venv /<path_to_new_env>
    source /<path_to_new_env>/bin/activate

Windows users should change the last line to ``\<path_to_new_env>\Scripts\activate.bat`` in order to activate their virtual environment.

With ``conda``:

.. code-block:: bash

    conda create -n nilearn python=3.9
    conda activate nilearn

**2. Install nilearn with pip**

Execute the following command in the command prompt / terminal
in the proper python environment:

.. code-block:: bash

    python -m pip install -U nilearn

Development version
-------------------

Please find all development setup instructions in the
`contribution guide <https://nilearn.github.io/stable/development.html#setting-up-your-environment>`_.

Check installation
------------------

Try importing nilearn in a python / iPython session:

.. code-block:: python

    import nilearn

If no error is raised, you have installed nilearn correctly.
