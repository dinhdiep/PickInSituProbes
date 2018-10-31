# DEPENDENCIES

bowtie1, R packages (dplyr, stringr, tidyr), MELTING5.1

# REFERENCE DATA

Need to download all these reference data files and put them into a 'Data' folder (add link to download here).

Multiple fasta file and its index (.fai) file:
*GRCh38_latest_genomic.fna
*GRCh38_latest_genomic.fna.fai

Bowtie1 index files:
*GRCh38.1.ebwt
*GRCh38.2.ebwt
*GRCh38.3.ebwt
*GRCh38.4.ebwt
*GRCh38.rev.1.ebwt
*GRCh38.rev.2.ebwt

Files for genes2targets.r:
*ncbi_refseq_GRCh38_cds.gtf.gz
*ncbi_refseq_GRCh38.gtf
*ref_conversion_table
*rsem_GRCh38.p2.gtf.gz

# USAGE:

## Step 1: generate target file 

First make a genes.txt file. Then in the directory with this file, run the script genes2targets.r using Rscript. Make sure that **data_dir** is set correctly to the correct path for your Data folder.

```{bash}
$ Rscript ./src/genes2targets.r
```

The outputs from this are (1) cds.txt, (2) shortest_isoforms.txt, and (3) constitutive.txt. Shortest_isoforms are generated for targets without cds regions (non-coding RNAs). Constitutive exons are generated for all the genes. Use either cds and shortest isoforms or constitutive to design probes.

## Step 2: modify the params file for your specific design

See params_example folder for examples of the parameters that can be set. The variables in the params file is written as a Perl script and used as input to the PickInSituProbes.pl program. Make sure that the correct path to your Data folder is set for the bowtie_cmd lines. Also, make sure that the correct paths are given to the MELTING5.1 and bowtie1 softwares.

## Step 3: run PickInSituProbes.pl

The first argument to this is the multiple fasta file (.fna or .fa) for the reference genome. The second argument is the job name. The third argument is the chosen params file. Then you also need to pipe in the genes target file that was generated in step 1. 

```{bash}
$ ./src/PickInSituProbes.pl Data/GRCh38_latest_genomic.fna MALAT1 params_example/our_optimal_params.pl < example.txt
```

The outputs from this that you'll want is the MALAT1.output (*job name*.output) file. The other files are for troubleshooting and evaluating the effects of different parameters on the outputs.


