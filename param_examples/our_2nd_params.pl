our $primerMaxTm = 75;
our $primerMinTm = 50; # 55 for 16mer and 17mer
our $primerMinLength = 18;
our $primerMaxLength = 24; #22 for 16mer and 17mer
our $gapSizeMin = 0;
our $gapSizeMax = 2;
our $deltaTm_cutoff = 15;

our $maxKmerOffTarget = 17;
our $minKmerOffTarget = 12;
our $minKmerOffTarget_max_occ = 800;

our $path_to_melting = "~/softwares/MELTING5.1.1/executable";
our $path_to_bowtie = "~/softwares/bowtie-1.1.1";

# Bowtie mappings
our $bowtie_filter_cmd = "$path_to_bowtie/bowtie -p 4 -f -a --strata --best Data/GRCh38 $job_name.kmer.fa > $job_name.kmer.bowtie.txt 2>> $job_name.kmer.log";
our $bowtie_find_off_target_cmd = "$path_to_bowtie/bowtie -p 4 -f -a -l 10 --best -n 2 Data/GRCh38 $job_name.oligo.fa > $job_name.bowtie.txt 2> $job_name.oligo.log";

# In M concentrations
our $melting_cmd = "java -cp $path_to_melting/melting5.jar melting.BatchMain -F 1 -E Na=0.33:formamide=2.5 -P 0.00000001 -H dnadna";
