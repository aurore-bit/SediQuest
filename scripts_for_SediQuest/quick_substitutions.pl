#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Basename;
use Math::CDF qw(:all);

my ($help, @ident, %data);

GetOptions ("help" => \$help);
&help if $help;

my @bams = @ARGV;
&help unless scalar @bams >= 1;

print join ("\t", "#ID", "5'CT", "3'CT", "5'CT_95CI", "3'CT_95CI", "5'#refC", "3'#refC",
                  "cond5'CT", "cond3'CT", "cond5'CT_95CI", "cond3'CT_95CI", "cond5'#refC", "cond3'#refC", "\n");

foreach my $bam (@bams) {
    %data = ();
    my $outname = basename($bam);
    $outname =~ s/\.bam$// or die "$bam is not a bam file\n";
    push (@ident, $outname);
    open READ, "samtools-0.1.18 view -X -F pu $bam |" or die "could not read file $bam\n";
    print STDERR "\nProcessing $bam\n";
    my ($counter, $interval);
    while (my $line = <READ>) {
        $counter++; $interval++;
        if ($interval == 10_000) {
            print STDERR "\r$counter";
            $interval = 0;
        }
        my @tmp = split /\t/, $line;
        my ($sequence, $locus, $aln_start, $mapqual, $orientation, $cigar, $header) = ("\U$tmp[9]", $tmp[2], $tmp[3], $tmp[4], $tmp[1], $tmp[5], $tmp[0]);
	my $md;
        foreach my $element (@tmp) {
            if ($element =~ /MD:[ZA]:(\S+)$/) {
                $md = $1;
            }
        }
        die "no MD field seen in\n $line\n" unless $md;
	$md =~ s/^0//; $md =~ s/(\D)0$/\1/; # remove leading or trailing zeros
	my $first_base = substr $sequence, 0, 1;
	next unless $first_base =~ /[ACGT]/;
	my $last_base = substr $sequence, -1;
	next unless $last_base =~ /[ACGT]/;
	my ($first_refbase, $last_refbase);
	my $first_md = substr $md, 0, 1;
	my $last_md = substr $md, -1;
	if ($first_md =~ /\d/) {
	    $first_refbase = $first_base;
	} elsif ($first_md =~ /([ACGT])/) {
	    $first_refbase = $1;
	} else {
	    next;
	    print "Problem: $md\n";
	}
	if ($last_md =~ /\d/) {
	    $last_refbase = $last_base;
	} elsif ($last_md =~ /([ACGT])/) {
	    $last_refbase = $1;
	} else {
	    next;
	    print "Problem: $md\n";
	}
	if ($orientation =~ /r/) {
	    $first_base =~ tr/ACGT/TGCA/;
	    $first_refbase =~ tr/ACGT/TGCA/;
	    $last_base =~ tr/ACGT/TGCA/;
	    $last_refbase =~ tr/ACGT/TGCA/;
	    ($first_base, $last_base) = ($last_base, $first_base);
	    ($first_refbase, $last_refbase) = ($last_refbase, $first_refbase);
	}

	my ($deam5, $deam3);
	if ($first_refbase eq "C") {
	    $data{"5refC"}++;
	    if ($first_base eq "T") {
		$data{"5seqT"}++;
		$deam5++;
	    }
	}
	if ($last_refbase eq "C") {
#	    print "$header\n";###debug
	    $data{"3refC"}++;
	    if ($last_base eq "T") {
		$data{"3seqT"}++;
		$deam3++;
	    }
	}
	if ($deam5 && $last_refbase eq "C") {
	    $data{"cond3refC"}++;
	    $data{"cond3seqT"}++ if $last_base eq "T";
	}
	if ($deam3 && $first_refbase eq "C") {
	    $data{"cond5refC"}++;
	    $data{"cond5seqT"}++ if $first_base eq "T";
	}
	
	#	print "$sequence\t$md\n$first_base\t$first_refbase\t$last_base\t$last_refbase\n\n";
    }
    my ($freq5, $freq3, $binom5, $binom3, $condfreq5, $condfreq3, $condbinom5, $condbinom3) = ("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA");
    my $ref5 = $data{"5refC"} || 0;
    my $ref3 = $data{"3refC"} || 0;
    my $seq5 = $data{"5seqT"} || 0;
    my $seq3 = $data{"3seqT"} || 0;
    my $condref5 = $data{"cond5refC"} || 0;
    my $condref3 = $data{"cond3refC"} || 0;
    my $condseq5 = $data{"cond5seqT"} || 0;
    my $condseq3 = $data{"cond3seqT"} || 0;
    $freq5 = sprintf ("%.1f", $seq5 / $ref5 *100) if $ref5;
    $freq3 = sprintf ("%.1f", $seq3 / $ref3 *100) if $ref3;
    $condfreq5 = sprintf ("%.1f", $condseq5 / $condref5 *100) if $condref5;
    $condfreq3 = sprintf ("%.1f", $condseq3 / $condref3 *100) if $condref3;
    if ($ref5) {my ($low, $high) = &binom($seq5, $ref5); $binom5 = "${low}-$high";}
    if ($ref3) {my ($low, $high) = &binom($seq3, $ref3); $binom3 = "${low}-$high";}
    if ($condref5) {my ($low, $high) = &binom($condseq5, $condref5); $condbinom5 = "${low}-$high";}
    if ($condref3) {my ($low, $high) = &binom($condseq3, $condref3); $condbinom3 = "${low}-$high";}
    print join ("\t", $outname, $freq5, $freq3, $binom5, $binom3, $ref5, $ref3, $condfreq5, $condfreq3, $condbinom5, $condbinom3, $condref5, $condref3), "\n";
}



sub report {
    foreach my $file (@ident) {
	print "$file\t";
	my ($freq5, $freq3, $ref5, $ref3) = ($data{$file}{"freq5"}, $data{$file}{"freq3"}, $data{$file}{"ref5"}, $data{$file}{"ref3"});
	if ($ref5) {
	    printf ("%.1f\t", $freq5);
	} else {
	    print "NA\t";
	}
	if ($ref3) {
	    printf ("%.1f\t", $freq3);
	} else {
	    print "NA\t";
	}
	print "$ref5\t";
	print "$ref3\t";
	if ($ref5 > 0) {
	    my $positive5 = sprintf("%.0f", $freq5/100*$ref5);
	    my $negative5 = $ref5 - $positive5;
	    my ($low5, $high5) = &binom($positive5, $ref5);
	    print "$low5", "-$high5\t";
	} else {
	    print "NA\t";
	}
	if ($ref3 > 0) {
	    my $positive3 = sprintf("%.0f", $freq3/100*$ref3);
	    my $negative3 = $ref3 - $positive3;
	    my ($low3, $high3) = &binom($positive3, $ref3);
	    print "$low3", "-$high3";
	} else {
	    print "NA\t";
	}
	print "\n";
    }
}

sub binom {
    my ($x, $n) = @_;
    return ("NA", "NA") unless $n > 0;
    my $alpha = 0.05;
    my $lower = qbeta($alpha / 2, $x, $n - $x + 1);
    my $upper = qbeta(1 - $alpha / 2, $x + 1, $n - $x);
    $lower = 0 if $x == 0;
    $upper = 1 if $x == $n;
    die "WTF is with $lower: $x/$n\n" unless defined $lower;
    die "WTF is with $upper: $x/$n\n" unless defined $upper;
    return (sprintf("%.1f", $lower * 100), sprintf("%.1f", $upper * 100));
}

sub help {
print "
This script quickly computes 5' and 3' substitution frequencies (incl. conditional substitutions) for alignment files in bam format. Paired and unmapped sequences are disregarded.

[usage]
./quick_substitutions.pl in1.bam in2.bam ...

[output]
-screen        regular and conditional substitution frequencies and 95% binomial confidence intervals

";
exit;
}
