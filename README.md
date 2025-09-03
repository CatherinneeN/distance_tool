# distance_tool
The tool is a ChimeraX extension that maps UniProt residues onto PDB structures using local alignment, then creates pseudobonds between specified residues. It calculates and logs Euclidean distances between two atoms, providing both visual and quantitative insight into residue proximities.

# Features

Local sequence alignment of UniProt to PDB residues
Maps residues even when sequence numbering differs between UniProt and PDB
Creates pseudobonds between crosslinked residues
Calculates Cα–Cα Euclidean distances and logs results
Visualises distance constraints directly in ChimeraX

# Installation

Save the plugin code to a local directory.
Open ChimeraX and run:
dev install /path/to/saved/code
Replace /path/to/saved/code with the folder where you stored the plugin files.

The tool will then be available as a ChimeraX command.

# Usage

Example command in ChimeraX:

hello


Inputs: FASTA sequences in one file and XL-MS TSV file containing residue pairs.

Outputs: Pseudobonds visualised in ChimeraX and logged distances between them.

