#!/usr/bin/perl -w
use strict;

my %geneTable;
my %probeCounts;
my $windowSize = 100;
my $max_size = 4;
my $header;
while(my $line = <STDIN>){
    chomp($line);
    if($line =~ m/GENE/){
	$header= $line;
	next;
    }
    my @probe_info = split ",", $line;
    next if($probe_info[9] =~ m/GGGGGG/ or $probe_info[9] =~ m/CCCCCC/ or $probe_info[9] =~ m/TTTTTT/ or $probe_info[9] =~ m/AAAAAA/);
    next if($probe_info[10] =~ m/GGGGGG/ or $probe_info[10] =~ m/CCCCCC/ or $probe_info[10] =~ m/TTTTTT/ or $probe_info[10] =~ m/AAAAAA/);
    next if(abs($probe_info[7] - $probe_info[8]) > 5);
    my ($gene, $probe_start) = ($probe_info[0], $probe_info[4]);
    my $bin = int($probe_start/$windowSize);
    push(@{$geneTable{$gene}->{$bin}}, $line);
    $probeCounts{$gene}++;
}

print $header,"\n";
foreach my $gene (keys %geneTable){
    my %max_Tm_by_bin = ();
    my %best_Probe_by_bin = ();
    # only 1 probe maximum per bin (each bin is 20 bp)
    foreach my $bin (keys %{$geneTable{$gene}}){
        $max_Tm_by_bin{$bin} = 0;
        my @candidates = @{$geneTable{$gene}->{$bin}};
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
        print $best_Probe_by_bin{$bin}, "\n";
    }
}
