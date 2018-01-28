#!/usr/bin/env perl

# (c) Victor Soria-Carrasco
# victor.soria.carrasco@gmail.com
# Last modified: 18/03/2016 15:19:40
# Based on code from Zach Gompert:
# This scripts converts a bcf file to a simpler format for downstream analysis. 
# I am calling this format genetoype likelihood (gl). The first line lists: the 
# number of individuals and loci. The next line has individual ids. This is followed
# by one line per SNP that gives the SNP id (scaffold, position) and the phred-scaled 
# genotype likelihoods, three per individual. 

# Changelog:
#  17/03/2016
#   - added support for bcftools 1.x
#  18/03/2016
#	- fixed bug when reading files
	
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Spec;

my $version='1.4.1-2016.03.18';

my $bcftools='bcftools';
# my $bcftools='/usr/local/extras/Genomics/apps/bcftools/1.2/bin/bcftools';

my $bcftools_version=`$bcftools |& grep Version | awk '{print \$2}'`;

my ($infile, $outfile, $outfileal);
my $entropy=0;
my $outal=0; # output alleles to file
GetOptions( 
    'i|I=s'  => \$infile, 
    'o|O=s'  => \$outfile,
    'a|A'    => \$outal,
    'e|E'    => \$entropy,
    'h|help' => \&usage 
)
or (print "\nERROR: Argument invalid\n" and &usage);

&usage if (!defined($infile));

$infile=File::Spec->rel2abs($infile);

# Determine if format is vcf or bcf
my $input=$infile;
if ($input =~ m/\.bcf$/i){
	if ($bcftools_version=~ m/^0\.1/){ # bcftools 0.1.x
		$input="$bcftools view -N $infile |";
	}
	else{ # bcftools v 1.2
		$input="$bcftools view -e 'REF=\"N\"' $infile |";
	}
}
elsif ($input =~ m/\.vcf\.gz$/i){
	$input="gzip -cd $infile |";
}
elsif ($input =~ m/\.vcf\.bz2$/i){
	$input="lbzip2 -cdk $infile |";
}
elsif ($input =~ m/\.vcf$/i){
	$input="cat $infile |";
}
else{
	die ("\nFormat of input file $infile not recognised\n\n");
}

my $nloc=`$input grep -v "^#" | grep -chv "[ATGC],[ATGC]"`;
chomp($nloc);

my $nind=`$input grep -m1 "^#CHROM" | sed 's/.*FORMAT[[:space:]]//' | awk '{print NF}'`;
chomp($nind);

print "\nNumber of loci: $nloc; number of individuals $nind\n";

if (!defined ($outfile)){
	$outfile=$infile;
	$outfile=~ s/\.[vb]cf(\.gz|\.bz2)?$/\.gl/;
}
$outfile=File::Spec->rel2abs($outfile);

if ($outal==1){
	$outfileal=$outfile;
	$outfileal=~ s/\.gl$/\.alleles\.txt/;
}

my @inds = ();
my $readgendata=0; # read genetic data
my @alleles=();
open (IN, "$input") or die "\nCan't read file $infile\n\n";
open (OUT, "> $outfile") or die "\nCan't write to $outfile\n\n";
	while (<IN>){
		chomp;
		## get individual ids
		if (m/^#CHROM/){
			my @aux=split(m/\s+/, $_);
			foreach my $i (9..$#aux){
				$aux[$i]=~ s/(\.sorted)?(\.[s|b]am)?//g;
				push (@inds, $aux[$i]);
				# $nind++;
			}
			print OUT "$nind $nloc\n";
			print OUT join (" ", @inds)."\n";
			print OUT join (" ", (('1') x scalar(@inds)))."\n" if ($entropy==0);
			
			$readgendata=1;
			next;
		}
		## read genetic data lines, write gl
		elsif ($readgendata && (!m/[AGCT],[AGCT]/)){ # second condition discards multiallelic snps
			my @aux=split(/\s+/,$_);
			my $id = "$aux[0]".":"."$aux[1]";
			print OUT "$id ";
			@aux = split(m/\s+/, $_);
			next if ($aux[4] =~ m/\,/); # discard multiallelic snps

			# get column index for genotype likelihoods from FORMAT field	
			my @aux2=split(/\:/,$aux[8]);
			my $plcol=0;
			while ($aux2[$plcol] ne "PL"){
				$plcol++;
			}
			# print "col $plcol - $aux2[$plcol]\n";
			foreach my $a (@aux[9..$#aux]){
				my @aux2=split(/\:/,$a);
				my $genotype=$aux2[0];
				my $genolhl=$aux2[$plcol];
				# print "$id - $a - genolhl $genolhl\n";
				if ($genolhl =~ /[0-9]+,[0-9]+,[0-9]+/){
					$genolhl =~ s/,/ /g;
					print OUT " $genolhl";
				}
				elsif ($genotype eq '0/0' && $genolhl eq '0'){ # invariant
					print OUT " 1 0 0";	
				}
				elsif ($genotype eq './.' && $genolhl eq '.'){ # missing data
					print OUT " 0 0 0";
				}
				else{
					die ("\nERROR: genotype ($genotype) and genotype likelihoods ($genolhl) format not recognised\n\n");
				}
			}
			print OUT "\n";
			push (@alleles, "$aux[3] $aux[4]") if ($outal==1);
		}	
	}
close (OUT);
close (IN);

print "\nFormat conversion finished. File saved as: $outfile\n\n";

if ($outal==1){
	open (OUTAL, ">$outfileal")
		or die ("\nCan't write to file $outfileal\n\n");
		print OUTAL join ("\n",@alleles)."\n";
	close (OUTAL);
	print "\n\tAlleles saved in: $outfileal\n\n";
}



# ==============================================================================
# ==============================================================================
# ============================== SUBROUTINES ===================================
# ==============================================================================
# ==============================================================================


# Show copyright
# ==============================================================================
sub author{
    print "\n";
    print "#########################################\n";
    print "  ".basename($0)."\n";
	print "  version $version\n";
    print "  (c) Victor Soria-Carrasco             \n";
    print "  victor.soria.carrasco\@gmail.com      \n";
    print "  based on code written by Zach Gompert \n";
    print "#########################################\n";
	print "\n";
}
# ==============================================================================

# Show usage
# ==============================================================================
sub usage{
    print "\n";
	print "  Usage:\n";
    print "    ".basename($0)."\n";
	print "      -i <input file (bcf file)>\n";
	print "      -o <output .gl file (optional)>\n";
	print "      -e <switch, output a slightly different format for entropy (optional,default=no)>\n";
	print "      -a <switch, output alleles to a file (optional,default=no)>\n";
	print "      -h this help\n";
    print "\n";
    exit;
}

