# MiniProject_1_code
A non-exhaustive depository of code for the first mini project of my PhD on the MIBTP program hosted by the University of Leicester and undertaken at the University of Warwick

Files are a series of scripts used for various tasks in the project including initial processing and generation of data, data analysis, generating plots and stats testing.

HPC_files - a directory containing scripts ran on an HPC. Generally for data processing and generation and more intense analysis of outputs. Some scripts were run both locally and on an HPC.

  1_RepairPDB.sh - Generates repaired pdb files from alphafold structures
  2_BuildModel.sh - Generates resistant structures from the repaired wildtypes built with   1_RepairPDB.sh
  AutodownloaderCMD.R - Downloads AA fasta files of wildtypes and names them using the RefGene Catalog
  FoldXProcessingR.R - Processing FoldX outputs to get summary tables of unfolding energy
  MutateXAnalysis - Processing for MutateX outputs to get summary tables of unfolding energy
  MutateXSubsetAnalysis.R - Subsettng MutateX outputs and processing to get summary tables.
  Mutationeditor.R - Creates fasta files for gyrA resistant point mutants based on the wildtype fasta files and mutation coordinates obtained from the RefGene Catalog
  ProFitAnalysis.R - Generates summary tables from ProFit output, where ProFit was run on FoldX structures
  ProFitRandom.sh - Runs ProFit on MutateX structures
  ProFitRandomAnalysis_pbp.R - Generates summary tables from ProFit output, where ProFit was run on MutateX structures
  ProFitSubsetAnalysis.R - Generates summary tables from ProFit output, where ProFit was run on MutateX structures. Takes only a subset of the output, using structures from around the same region as the mean mutation position for that protein for that species 
  pbpmutationeditor.R - Creates fasta files for pbp resistant point mutants based on the wildtype fasta files and mutation coordinates obtained from the RefGene Catalog



Local_files - a directory containing scripts ran on my local machine. Typically for simpler data processing and generation, simple calculations, stats and plotting.

  Autodownloader.R - Downloads AA fasta files of wildtypes and names them using the RefGene Catalog
  More_plots.R - Generates plots for some outputs
  Residue_counter.R - Counts conservative, near conservative, and radical point mutations
  Scriptforsomegraphs.R - Generates plots for some outputs
  Sequence_Logos.R - Generates sequence logos and heatmaps for conservative/radical point mutations
  Stats_testing.R - Stats, overwritten multiple times so only the latest test is shown.
  mutantfilemaker.R - Produces mutant files needed for ProFit to run.
  pLDDT_scores.R - Finds and produces a dataframe of pLDDT scores for AlphaFold Structures.
