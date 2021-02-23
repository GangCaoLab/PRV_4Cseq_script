#!/usr/bin/perl -w
use strict;
use warnings;

opendir(DIR,"./")or die $!;
my @file=readdir DIR;
@file = grep (/(.fastq)$/gi, @file);
mkdir "./cutedread";
foreach my $primer(@file){
    open(IN, "./$primer") or die $!;
    my $n=0;
    my @line;
    my $en1;
    my @en1pos;
    my $en2;
    my @en2pos;
    while (<IN>) {
        chomp;
        $n++;
        push(@line,$_);
        if ($n%4==2) {
            if ($_=~m/GGATCC/g) {
                $en1++;
                push(@en1pos,pos($_));
            }elsif($_=~m/AGCT/g){
                $en2++;
                push(@en2pos,pos($_));
            }
        }
    }
    close IN;
    my $site=0;
    if ($en1>3*$en2) {
        $site=&getpos(\@en1pos);
        &splitread(\@line,$site);
    }elsif($en1<3*$en2){
        next;
    }else{
        print "data($primer) has same problem";
        next;
    }
    
    sub getpos{
        my $posarr=@_;
        my %hash;
        foreach my $posele(@$posarr){
            $hash{$posele}++;
        }
        my $n=0;
        my $pos=0;
        foreach my $key(sort {$hash{$b}<=>$hash{$a}} keys %hash){
            $n++;
            if ($n==1) {
                $pos=$key;
                last;
            }
        }
        return $pos;
    }
    
    sub splitread{
        my ($line,$pos)=@_;
        my @read=@$line;
        open(OUT, ">./cutedread/$primer.cuted.fastq") or die $!;
        for(my $i=0;$i<@read;$i++){
            if (($i+1)/2==0) {
                print OUT substr($read[$i],$pos);
            }else{
                print OUT $read[$i];
            }
        }
        close OUT;
    }
}


