# METAMAP - Metagenome analysis with Megan and Python

[![Build Status](https://travis-ci.com/lucass122/METAMAP.svg?branch=master)](https://travis-ci.com/lucass122/METAMAP)

Collection of scripts for metagenomic downstream analysis of bioreactor data using the tool MEGAN and the programming language Python.

# Installation

NOTE: REQUIREMENTS NOT ADDED YET PLEASE DOWNLOAD SCRIPTS MANUALLY FOR NOW TO USE THEM - WILL UPDATE THIS ONCE PROJECT IS FINISHED

1. Clone this repository
2. Navigate to repository
3. pip install .

# Taxonomic analysis

```
taxonomy_plots.py
```

reads in multiple .daa files containing metagenomic reads mapped to NR database, extracts taxonomic information from them with MEGAN and then uses python's matplotlib to visualise the data

# Functional analysis


```
gc_assembly_map.py
```

takes a .daa file as input as well as a SEED ID and computes a list with organisms in the .daa file that have the gene described by the ID. It uses gene centric assembly included with MEGAN as well as NCBI online Blast


```
extract_strain_table.py
```

Takes BLAST results from gc_assembly_map and processes them. The output is a table with all enzymes that were targeted by the gc-assembly and the bacterial strains they are associated with in the .daa file
