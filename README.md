# stafix


**scaleSTAFIX_amber.py** and **scaleSTAFIX_gromacs.py** are scripts to stafix-scale **Amber** and **Gromacs** topologies. epsilon for Lennard-Jones interactions of specifc RNA atoms is scaled by user-specified factor. More details can be found here: https://doi.org/10.1101/2022.07.22.501120 

requirements:

python 3 (remove input/output type definitions in is_to_be_scaled() function to use this also with python 2)

parmed - for scaleSTAFIX_amber.py (scaling is done via parmed tool of Ambertools)

usage:

python scaleSTAFIX_amber.py parm_file scaling_factor [mask]
	
[mask] is to be given in format resid1,resid2,resid3,... (e.g. 1,2,3,4). If no mask is specified, all RNA residues will be scaled

example:

	python scaleSTAFIX_amber.py complex.parm7 0.5 1,2,3,4,5,6

	python scaleSTAFIX_gromacs.py complex.top 0.75

do not use scaleSTAFIX_amber.py for scaling topologies that are to be converted to gromacs, use scaleSTAFIX_gromacs.py on gromacs topology instead!

flows:

scaleSTAFIX_amber.py:

run parmed printDetails for original topology file, then read the output and find atoms that are to be scaled (based on atomname and residuename) and get their IDs and LJ parameters
        
run parmed again this time with addLJType to make new type for atoms that we scale and then changeLJpair to scale epsilon for the selected atom pairs

scaleSTAFIX_gromacs.py:

read original topology, get definitions of all present atomtypes from there and write a temporary topology that will have atomtypes of atoms that are to be scaled changed by appending "Y" to their atomtypename

pick atomtypes that are to be scaled

write new topology where [ nonbond_params ] section is added after [ atomtypes ] and contains definition of LJ params for scaled atomtypes (those with "Y" appended)


**scaleSTAFIX.py** is script to stafix-scale **Amber** topology with parmed tool of Ambertools

requirements:

parmed and python 3 (remove input/output type definitions in is_to_be_scaled function to use this also with python 2)

usage:

python scaleSTAFIX.py parm_file scaling_factor [mask]
	
[mask] is to be given in format resid1,resid2,resid3,... (e.g. 1,2,3,4). If no mask is specified, all RNA residues will be scaled

flow:

run parmed printDetails for original topology file, then read the output and find atoms that are to be scaled (based on atomname and residuename) and get their IDs and LJ parameters
	
run parmed again this time with addLJType to make new type for atoms that we scale and then changeLJpair to scale epsilon for the selected atom pairs
  
example:

	python scaleSTAFIX.py complex.parm7 0.5 1,2,3,4,5,6

do not use this for scaling topologies that are to be converted to gromacs!

