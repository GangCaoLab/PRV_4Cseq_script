#!/usr/bin/perl -w
use strict;
use warnings;

my $exp=$ARGV[0];

my $cutoff=1;
my $fragLen = 250;
open(OUT, ">exp$exp.$fragLen-$cutoff.cluster") or die $!;
open(IN, "0.exp$exp.bed") or die $!;
my @data;
while (<IN>) {
    chomp;
    if ($_=~/^NC_010/) {
        my @term=split(/\t/,$_);
        my $mid=int(($term[2]+$term[1])/2);
        $term[0]=~s/NC_010//;#NC_010443.4
        $term[0]=int($term[0])-442;
        
        #print $term[0];
        push(@term,$mid);
        push(@data,\@term);
    }
}
close IN;

sub campare {
    $a->[0] <=> $b->[0] or $a->[4] <=> $b->[4] or $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2]
}
@data = sort campare @data;
#print OUT "CHR\t"."MID1\t"."MID2\t"."LEN\t"."COUNT\t"."UNIQ\t"."START\t"."END\t"."CHR_ID\t"."PRIMER\n";
my $count=0;
my $mid1;
my $mid2;
my $midx;
my $chr;
my $n;
my $id;
my %hash;
my %hashp;
my $s;
my $e;
foreach my $subdata (@data) {
    $n++;
    if ($n==1) {
        $count=1;
        $mid1=$subdata->[4];
        $chr=$subdata->[0];
        my $key=$subdata->[1]."_".$subdata->[2];
        $hash{$key}++;
        $s=$subdata->[1];
        $e=$subdata->[2];
        next;
    }
    $mid2=$subdata->[4];
    if (abs($mid2-$mid1)<=$fragLen && $chr==$subdata->[0]) {
        $count++;
        $midx=$mid2;
        my $key=$subdata->[1]."_".$subdata->[2];
        $hash{$key}++;
        $subdata->[3]=~/(\d+).\d/;
        my $primernum=$1;
        $hashp{$primernum}++;
        if ($subdata->[1]<=$s) {
            $s=$subdata->[1];
        }
        if ($subdata->[1]>=$s) {
            $e=$subdata->[2];
        }
    }else{
        if ($count>=$cutoff) {
            $id++;
            my $len=$midx-$mid1;
            my $size= keys %hash;
            if ($chr==1){ $chr="NC_010443.4";
            }elsif($chr==2){ $chr="NC_010444.3";
            }elsif($chr==3){ $chr="NC_010445.3";
            }elsif($chr==4){ $chr="NC_010446.4";
            }elsif($chr==5){ $chr="NC_010447.4";
            }elsif($chr==6){ $chr="NC_010448.3";
            }elsif($chr==7){ $chr="NC_010449.4";
            }elsif($chr==8){ $chr="NC_010450.3";
            }elsif($chr==9){ $chr="NC_010451.3";
            }elsif($chr==10){ $chr="NC_010452.3";
            }elsif($chr==11){ $chr="NC_010453.4";
            }elsif($chr==12){ $chr="NC_010454.3";
            }elsif($chr==13){ $chr="NC_010455.4";
            }elsif($chr==14){ $chr="NC_010456.4";
            }elsif($chr==15){ $chr="NC_010457.4";
            }elsif($chr==16){ $chr="NC_010458.3";
            }elsif($chr==17){ $chr="NC_010459.4";
            }elsif($chr==18){ $chr="NC_010460.3";
            }elsif($chr==19){ $chr="NC_010461.4";
            }elsif($chr==20){ $chr="NC_010462.2";
            }
            my $sizex= keys %hashp;
            print OUT "$chr\t$mid1\t$midx\t$len\t$count\t$size\t$s\t$e\t$subdata->[0]\t$sizex";
            #foreach my $key(sort {$a<=>$b} keys %hashp){
                #print OUT "\t$key";
            #}
            print OUT "\n";
        }
        %hash=();
        %hashp=();
        $s=$subdata->[1];
        $e=$subdata->[2];
        my $key=$subdata->[1]."_".$subdata->[2];
        $hash{$key}++;
        $mid1=$subdata->[4];
        $chr=$subdata->[0];
        $count=1;
    }
}
close OUT;
