#!/usr/bin/perl -w
use threads;

#Copyright 2018 The Regents of the University of California.  All Rights Reserved.

our %kmerTable;

our $genome_fa = $ARGV[0];
our $job_name = $ARGV[1];


our $primerMaxTm = 75;
our $primerMinTm = 55; # 55 for 16mer and 17mer
our $primerMinLength = 18;
our $primerMaxLength = 22; #22 for 16mer and 17mer
our $gapSizeMin = 1;
our $gapSizeMax = 2;
our $deltaTm_cutoff = 20;

our $maxKmerOffTarget = 16;
our $minKmerOffTarget = 12;
our $minKmerOffTarget_max_occ = 1000;

our $path_to_melting = "~/softwares/MELTING5.1.1/executable";
our $path_to_bowtie = "~/softwares/bowtie-1.1.1";

# Bowtie mappings
our $bowtie_filter_cmd = "$path_to_bowtie/bowtie -p 4 -f -a --strata --best GRCh38 $job_name.kmer.fa > $job_name.kmer.bowtie.txt 2>> $job_name.kmer.log";
our $bowtie_find_off_target_cmd = "$path_to_bowtie/bowtie -p 4 -f -a -l 10 --best -n 2 GRCh38 $job_name.oligo.fa > $job_name.bowtie.txt 2> $job_name.oligo.log";

# In M concentrations
our $melting_cmd = "java -cp $path_to_melting/melting5.jar melting.BatchMain -F 1 -E Na=0.33:formamide=2.5 -P 0.00000001 -H dnadna";


my $input_file = $ARGV[2];
eval `cat $input_file` or die "couldn't parse input file";

our %cmpTable = ("A", "T", "T", "A", "G", "C", "C", "G", "N", "N", "M", "K", "R", "Y", "W", "W", "S", "S", "Y", "R", "K", "M", "V", "B", "H", "D", "D", "H", "B", "V");
our %oligoTm;

sub main{

    my $start_time = time;

    my %geneTable;

    while(my $line = <STDIN>){
        chomp($line);
        my @fields = split "\t", $line;
        my ($geneName, $chrName, $chrStart, $chrEnd, $strand) = ($fields[0], $fields[1], $fields[2], $fields[3], $fields[4]);
        push(@{$geneTable{$geneName}}, $strand . "\t" . $chrName . "\t" . $chrStart . "\t" . $chrEnd); 
    }

    print "Found ", scalar(keys %geneTable), " genes.\n";
    my $oligo_id = 1;
    my @oligoSet;
    foreach my $gene (keys %geneTable){
        my @candidates = @{$geneTable{$gene}};
        my $num_candidate_probes = 0;
        for(my $i = 0; $i < scalar(@candidates); $i++){
            my ($strand, $chr, $startPos, $endPos) = split "\t", $candidates[$i];
            my $target_seq = getRNASeq($strand, $chr, $startPos, $endPos);
            for(my $j = 0; $j < length($target_seq) - $primerMinLength; $j++){
                for($k = $primerMinLength; $k <= $primerMaxLength; $k++){
                    my $oligo = substr($target_seq, $j, $k);
                    my $gc = $oligo =~ tr/gcGC//;
                    next if($oligo =~ m/GGGGGGG/ or $oligo =~ m/CCCCCCC/ or $oligo =~ m/AAAAAAA/ or $oligo =~ m/TTTTTTT/);
                    next if($gc/length($oligo) < 0.45 or $gc/length($oligo) > 0.70);
                    addKmer($oligo, $maxKmerOffTarget, $oligo_id);
                    push(@oligoSet, $j . "," . $k . "," . $oligo . "," . $candidates[$i]);
                    $oligo_id++;
                }
            }
        }
    }

    generateKmerFasta();
    system($bowtie_filter_cmd);
    my @filtered_max_kmer_off_target = removeCandidates(\@oligoSet, 1);

    # removes those with 10-mer greater than 1000 positions
    $oligo_id = 1;
    for(my $i = 0; $i < scalar(@filtered_max_kmer_off_target); $i++){
        my ($j, $k, $oligo, $candidate) = split "," , $filtered_max_kmer_off_target[$i];
        addKmer($oligo, $minKmerOffTarget, $oligo_id);
        $oligo_id++;
    }

    generateKmerFasta();
    system($bowtie_filter_cmd);
    my @filtered_max_kmer_and_10mer_off_target = removeCandidates(\@filtered_max_kmer_off_target, $minKmerOffTarget_max_occ);
    unlink("$job_name.kmer.fa");
    unlink("$job_name.kmer.bowtie.txt");
    unlink("$job_name.kmer.log");

    my %filtered_oligoSet;

    for(my $i = 0; $i < scalar(@filtered_max_kmer_and_10mer_off_target); $i++){
        my ($j, $k, $oligo, $candidate) = split "," , $filtered_max_kmer_and_10mer_off_target[$i];
        $oligoTm{$oligo} = 1;
        push(@{$filtered_oligoSet{$candidate}}, $j . "\t" . $k . "\t" . $oligo);
    }
    
    my @oligoList = keys %oligoTm;
    generateOligoSeq(\@oligoList);
    computeMeltingTemp("$job_name.1.seq", "$job_name.1.seq.results.csv");

    open(MELTING, "$job_name.1.seq.results.csv") or die("Error reading melting table\n");
    while(my $line = <MELTING>){
        chomp($line);
        my @fields = split "\t", $line;
        next if(!$oligoTm{$fields[0]});
        $oligoTm{$fields[0]} = $fields[4];
    }
    close(MELTING); 
    
    my @candidateProbes;
    foreach my $gene (keys %geneTable){
        my @candidates = @{$geneTable{$gene}};
        my $num_candidate_probes = 0;
        for(my $i = 0; $i < scalar(@candidates); $i++){
            next if(!$filtered_oligoSet{$candidates[$i]});
            my ($strand, $chr, $startPos, $endPos) = split "\t", $candidates[$i];
            my @oligoSet = @{$filtered_oligoSet{$candidates[$i]}};
            my ($Tm_1, $Tm2);
            for(my $m = 0; $m < scalar(@oligoSet) - 1; $m++){
                my ($j_1, $k_1, $a_oligo) = split "\t", $oligoSet[$m];
                $Tm_1 = $oligoTm{$a_oligo};
                next if($Tm_1 < $primerMinTm or $Tm_1 > $primerMaxTm);
                for(my $n = $m + 1; $n < scalar(@oligoSet); $n++){
                    my ($j_2, $k_2, $b_oligo) = split "\t", $oligoSet[$n];
                    my $candidate_start = $startPos + $j_1;
                    # 10-12-18 Dinh: fix candidate end pos "-1" offset
                    my $candidate_end = $startPos + $j_2 + $k_2 - 1;
                    my $gapSize = $candidate_end - $candidate_start + 1 - length($a_oligo) - length($b_oligo);
                    next if($gapSize > 2 or $gapSize < 1); # skip if gapsize is less than 1 or greater than 2
                    $Tm_2 = $oligoTm{$b_oligo};
                    next if($Tm_2 < $primerMinTm or $Tm_2 > $primerMaxTm);
                    push(@candidateProbes, "$gene,$chr,$startPos,$endPos,$candidate_start,$candidate_end,$gapSize,$Tm_1,$Tm_2,$a_oligo,$b_oligo," . join(",", getSNAIL($a_oligo, $b_oligo, "_DECODING_")) );
                    $num_candidate_probes++;
                }
            }
        } 
        #print "Found ", $num_candidate_probes, " candidate probes for $gene.\n";
    }

    # Rank probes using off target annealing Tm difference
    generateOligoFasta(\@candidateProbes);
    system($bowtie_find_off_target_cmd); 
    my @annotated_kmer_off_target_mappings = annotateCandidates(\@candidateProbes);
     
    open(OUT, ">$job_name.output") or die("Error writing $job_name.output");
    print OUT "GENE,REF_CHROM,REF_START,REF_END,HYB_START,HYB_END,GAPSIZE,PRIMER_TM,PADLOCK_TM,OLIGO_A,OLIGO_B,SNAIL_PRIMER,SNAIL_PADLOCK,DELTA_TM\n";
    my %targetSet;
    my $windowSize = 100;
    my $max_size = 4;

    foreach my $value (@annotated_kmer_off_target_mappings){
        my @probe_info = split ",", $value;
        my $maxDeltaTm = $probe_info[13];
        my $totalTm = $probe_info[7] + $probe_info[8];
        my $target = $probe_info[0] . "," . $probe_info[1] . "," . $probe_info[2] . "," . $probe_info[3];
        next if($maxDeltaTm < $deltaTm_cutoff);
        next if(abs($probe_info[7] - $probe_info[8]) > 5);
        my ($gene, $probe_start) = ($probe_info[0], $probe_info[4]);
        my $bin = int($probe_start/$windowSize);
        push(@{$targetSet{$gene}->{$bin}}, $value);

    }

    foreach my $gene (keys %targetSet){
        my %max_Tm_by_bin = ();
        my %best_Probe_by_bin = ();
        # only 1 probe maximum per bin (each bin is 20 bp)
        foreach my $bin (keys %{$targetSet{$gene}}){
            $max_Tm_by_bin{$bin} = 0;
            my @candidates = @{$targetSet{$gene}->{$bin}};
            for(my $i = 0; $i < scalar(@candidates); $i++){
                my $probe_line = $candidates[$i];
                my @probe_info = split ",", $probe_line;
                if($probe_info[13] > $max_Tm_by_bin{$bin}){
                    $max_Tm_by_bin{$bin} = $probe_info[13];
                    $best_Probe_by_bin{$bin} = $probe_line;
                }
            }
        }
        my $num = 0;
        foreach my $bin (sort { $max_Tm_by_bin{$a} <=> $max_Tm_by_bin{$b} } keys %max_Tm_by_bin){
            last if($num == $max_size);
            $num++;
            print OUT $best_Probe_by_bin{$bin}, "\n";
        }
    }

    close(OUT);

    my $end_time = time;
    print "Total time: ", $end_time - $start_time, "\n";
}

sub computeMeltingTemp{
    my $input = shift;
    my $output = shift;
    open(IN, "$input") or die("Error opening $input\n");
    open(TMP, ">$input.tmp") or die("Error writing $input.tmp\n");
    open(OUT, ">$output") or die("Error writing $output\n");
    close(OUT);
    my $count = 0;
    while(my $line = <IN>){
        if($count % 5000 == 1){
            close(TMP);
            system("$melting_cmd $input.tmp | sort -u >> melting.$input.log");
            system("cat $input.tmp.results.csv >> $output");
            open(TMP, ">$input.tmp") or die("Error writing $input.tmp\n");
        }
        print TMP $line;
        $count++;
    }
    close(IN);
    # process last lines
    close(TMP);
    system("$melting_cmd $input.tmp | sort -u >> melting.$input.log");
    system("cat $input.tmp.results.csv >> $output");
    unlink("$input");
    unlink("$input.tmp");
    unlink("$input.tmp.results.csv");
    unlink("melting.$input.log");
}


sub getRNASeq{
    my ($strand, $chr, $start, $end) = @_;
    my @seq_lines = `samtools faidx $genome_fa $chr:$start-$end | grep -v [0-9]`;
    my $target_seq = join('', map { chomp; $_ } @seq_lines);
    if($strand eq '-'){
        $target_seq = RevComp($target_seq);
    }
    return $target_seq;
}


sub getAntiSense{
    my $rna_seq = shift;
    my $antisense_seq = RevComp($rna_seq);
    return $antisense_seq;
}


sub getSNAIL{
    my $oligo_a = shift;
    my $oligo_b = shift;
    my $decoding = shift;
    my ($anti_a, $anti_b) = (getAntiSense($oligo_a), getAntiSense($oligo_b));
    my $primer = $anti_b . "TAATGTTATCTT";
    my $padlock = "ACATTA" . $anti_a . "ATTATTA" . $decoding . "AAGATA";
    return ($primer, $padlock);
}


sub annotateCandidates{
    my $candidates = shift;
    my @annotated_candidates;
    my %countHits;
    open(IN, "$job_name.bowtie.txt") or die("Error reading bowtie table\n");
    while(my $line = <IN>){
        chomp($line);
        my @fields = split "\t", $line;
        #2       +       NC_000001.11    33493342        AGCAGTGTCG      IIIIIIIIII      851
        my ($strand, $chr, $startPos, $endPos) = ($fields[1], $fields[2], $fields[3]+1, $fields[3]+length($fields[4])); 
        my $oligo = $fields[4];
        $oligo = RevComp($fields[4]) if($strand eq '-');
        if($countHits{$oligo}->{"counts"}){
            $countHits{$oligo}->{"counts"}++;
            if($fields[7]){
                my @mm_list = split ",", $fields[7];
                next if($countHits{$oligo}->{"minMM"} < scalar(@mm_list));
                my $off_target_seq = flipDNA( getRNASeq($strand, $chr, $startPos, $endPos) );
		# assumes 5C difference per base pair mismatch
                $countHits{$oligo}->{"off_target_Tm"} = $oligoTm{$oligo} - 5*scalar(@mm_list);
                my $nc = $off_target_seq =~ tr/Nn//;
                next if($nc > 0);
                push(@{$countHits{$oligo}->{"off_target_seq"}}, $off_target_seq);
                $countHits{$oligo}->{"minMM"} = scalar(@mm_list);
            }
        }else{
            # first hit always perfect
            $countHits{$oligo}->{"counts"} = 1;
            $countHits{$oligo}->{"minMM"} = 10; # set maximum mismatches
            $countHits{$oligo}->{"off_target_Tm"} = 0; # set maximum mismatches
        }
    }
    close(IN);

    open(OUT, ">$job_name.2.seq") or die("Error writing seq file \n");
    foreach my $value (keys %countHits){
        next if($countHits{$value}->{"counts"} == 1);
        next if($countHits{$value}->{"minMM"} == 0 or $countHits{$value}->{"minMM"} == 10);
        my @off_target_seqs = @{$countHits{$value}->{"off_target_seq"}};
        foreach my $seq (@off_target_seqs){
            print OUT $value, "\t", $seq, "\n";
        }
    }
    close(OUT);

    computeMeltingTemp("$job_name.2.seq", "$job_name.2.seq.results.csv");

    open(IN, "$job_name.2.seq.results.csv") or die("Error reading melting results \n");
    while(my $line = <IN>){
        chomp($line);
        my @fields = split "\t", $line;
        next if(!$countHits{$fields[0]});
        if($countHits{$fields[0]}->{"off_target_Tm"} < $fields[4]){
            $countHits{$fields[0]}->{"off_target_Tm"} = $fields[4];
        }
    }
    close(IN);

    for(my $i = 0; $i < scalar(@{ $candidates }); $i++){
        my @probe_info = split ",", $candidates->[$i];
        my $a_oligo = $probe_info[9];
        my $b_oligo = $probe_info[10];
        my $minDeltaTm = 0;
        my $deltaTm_a = 1000;
        $deltaTm_a = $probe_info[7] - $countHits{$a_oligo}->{"off_target_Tm"};
        my $deltaTm_b = 1000;
        $deltaTm_b = $probe_info[8] - $countHits{$b_oligo}->{"off_target_Tm"}; 
        $minDeltaTm = $deltaTm_a;
        if($deltaTm_a > $deltaTm_b){
            $minDeltaTm = $deltaTm_b;
        }

        next if($minDeltaTm < 5);
        push(@annotated_candidates, $candidates->[$i] . ",$minDeltaTm");
    }
    return @annotated_candidates;
}


sub generateOligoFasta{
    my $candidates = shift;
    my %oligoTable;
    open(OUT, ">$job_name.oligo.fa") or die("Error writing $job_name.oligo.fa file\n");
    for(my $i = 0; $i < scalar(@{ $candidates }); $i++){
        my @probe_info = split ",", $candidates->[$i];
        $oligoTable{$probe_info[9]} = 1;
        $oligoTable{$probe_info[10]}= 1;
    }
    my $cur_candidate = 1;
    foreach my $oligo (keys %oligoTable){
        print OUT ">$cur_candidate\n";
        print OUT $oligo, "\n";
        $cur_candidate++;
    }
    close(OUT);
}


sub generateOligoSeq{
    my $oligoList = shift;
    open(OUT, ">$job_name.1.seq") or die("Error writing to $job_name.1.seq");
    print OUT join("\n", @{ $oligoList }), "\n";
    close(OUT);
}


sub removeCandidates{
    my $candidates = shift;
    my $cutoff = shift;
    my @filtered_candidates;
    my %removeIds;
    my %countHits;
    open(IN, "$job_name.kmer.bowtie.txt") or die("Error reading kmer bowtie table");
    while(my $line = <IN>){
        chomp($line);
        my @fields = split "\t", $line;
        next if($fields[2] =~ m/NW_/);
        #2       +       NC_000001.11    33493342        AGCAGTGTCG      IIIIIIIIII      851
        my $kmer = uc($fields[4]);
        $kmer = RevComp($kmer) if($fields[1] eq '-');
        $countHits{$kmer}++;
    }
    close(IN);

    foreach my $kmer (keys %countHits){
        next if(!$kmerTable{$kmer});
        if($countHits{$kmer}){
            next if($countHits{$kmer} <= $cutoff);
            my @kmer_associated_candidates = split ",", $kmerTable{$kmer};
            foreach my $value (@kmer_associated_candidates){
                $removeIds{"candidate_" . $value} = 1;
            }
        }
        delete($kmerTable{$kmer});
    }
    
    %kmerTable = ();

    print scalar(@{ $candidates }), "=num candidates \n";
    print scalar(keys %removeIds), "=num removes \n";
    for(my $i = 0; $i < scalar(@{ $candidates }); $i++){
        my $id = $i+1;
        next if($removeIds{"candidate_" . $id});
        push(@filtered_candidates, $candidates->[$i]);
    }
    return @filtered_candidates;
}


sub generateKmerFasta{
    my $num = 1;
    open(OUT, ">$job_name.kmer.fa") or die("Error writing to $job_name.kmer.fa");
    foreach my $kmer (keys %kmerTable){
        print OUT ">", $num, "\n";
        print OUT $kmer, "\n";
        $num++;
    }
    close(OUT);
}


sub addKmer{
    my $oligo = shift;
    my $kmer_len = shift;
    my $candidate_id = shift;
    for(my $i = 0; $i < length($oligo) - $kmer_len; $i++){
        my $cur_kmer = uc(substr($oligo, $i, $kmer_len));
        if($kmerTable{$cur_kmer}){
            $kmerTable{$cur_kmer} = $kmerTable{$cur_kmer} . "," . $candidate_id;
        }else{
            $kmerTable{$cur_kmer} = $candidate_id;
        }
    }
}


sub RevComp{
    my $seq = shift;
	my $seqLen = length($seq);
	my $revcom = '';
	for(my $i = 0; $i < $seqLen; $i++){
		#For any base that is not A/T/G/C, such as a degenerate base, use "N".
		my $rcBase = $cmpTable{substr($seq, $i, 1)} ? $cmpTable{substr($seq, $i, 1)} :'N';
		$revcom = $rcBase . $revcom;
	}
    return $revcom;
}


sub flipDNA{
    my $oligo = shift;
    my $UColigo = uc($oligo);
    $UColigo =~ tr/G/c/;
    $UColigo =~ tr/C/g/;
    $UColigo =~ tr/A/t/;
    $UColigo =~ tr/T/a/;
    my $newoligo = uc($UColigo);
    return $newoligo;
}

main();
