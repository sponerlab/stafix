# stafix

**scaleSTAFIX.py** is script to stafix-scale **Amber** topology with parmed tool of Ambertools\n
requirements:
	parmed and python 3 (remove input/output type definitions for is_to_be_scaled function to use this also with python 2)
usage:
	python scaleSTAFIX.py parm_file scaling_factor [mask]
		[mask] is to be given in format resid1,resid2,resid3,... (e.g. 1,2,3,4). If no mask is specified, all RNA residues will be scaled
flow:
	run parmed printDetails for original topology file, then read the output and find atoms that are to be scaled (based on atomname and residuename) and get their IDs and LJ parameters
	run parmed again this time with addLJType to make new type for atoms that we scale and then changeLJpair to scale epsilon for the selected atom pairs
  
do not use this for scaling topologies that are to be converted to gromacs!

example: python scaleSTAFIX.py complex.parm7 0.5 1,2,3,4,5,6

