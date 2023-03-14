# BHMC
## Introduction
A workflow to generate the structures of cluster system by Basin-Hopping Monte Carlo algorithm
## Install necessary packages
It is suggested to create a virtual environment via conda for your project.  After your conda environment is activated, run the following command to install necessary packages.
、、、
    pip install dpdata
    pip install pyxtal
    pip install ase
    pip install dpdispatcher
、、、
## Basic usage
You have to prepare an input.json file, which can refer to template input files in the input-example folders. Here, we only offer the example which use cp2k or Deep-SE model to run geometry optimizations. Other calculators will be given in the future. Now you can run the workflow through nohup commands.
、、、
    nohup python run/BHMC.py 1>log 2>err &
、、、
