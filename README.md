# FeatureFinder

General features are first identified from a given structure file (PDB). 
This includes basic information such as whether the system contains proteins, nucleic acid chains, small molecules, etc. 
It narrows down this information further to identify likely analyses which will be performed later.
As an example, if a protein and small molecule are both identified in the starting structure, FeatureFinder will recommend binding energy calculations.
In situations where an unrecognized small molecule is present, it will further analyze the structure to determine features or substructures such as carbohydrate scaffolds or modified nucleotides.

### AGIMUS Integration

When installed as part of the AGIMUS suite, any unrecognized structures such as amino acids that have undergone post-translational modifications or small molecule drug candidates are identified, isolated, and passed to AutoParams.
