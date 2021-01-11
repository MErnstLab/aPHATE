# aPHATE

This package contains the script for performing the modified implementation of PHATE intended primarily for analysis of computational docking results.  
For a more in depth explanation of the method please see the Frontiers Neuroscience publication.  

## Quick start  
The data should be packed into a list, where:  
- the entry `Circular` contains a matrix of circular coordinates of atoms, where each row represents a single molecule pose,  
- the entry `Samples` contains a factor with sample information for every column in `Circular`,  
- the entry `Scoring` contains a factor with scores for every column in `Circular`.
