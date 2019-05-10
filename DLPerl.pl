#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
#Argument handling
my $input;
my $source = 'none';
my $mtf_file;
my $help = '';
GetOptions('source=s' => \$source, 'input=s' => \$input, 'motifs=s' => \$mtf_file, 'help' => \$help) or die "Usage: $0 --input NAME\n";
(my $fname = $input) =~ s{.*/}{};
$fname =~ s/\.[^.]+$//;
#Display Header
print <<HEADER;
 ____  _     ____           _ 
|  _ \\| |   |  _ \\ ___ _ __| |
| | | | |   | |_) / _ \\ '__| |
| |_| | |___|  __/  __/ |  | |
|____/|_____|_|   \\___|_|  |_|

HEADER
#Display Help if required
if($help ne '')
{
	print "DLPerl gives information, as gff and tsv files, about GC, ATGC, N content, motif presence and its location, from FASTA files.\n\nList of commands:\n\t--input FILENAME\tFASTA file to be analyzed\n\t--source\t\tOrigin of the FASTA file (default: none)\n\t--motifs FILENAME\tFile containing motifs in FASTA format\n\t--help\t\t\tDisplays this information\n";
	exit;
}
#Motif variable initialization
my @find;
my @id;
my @name;
my @chain = ("+", "-");
my @count;
#Motif file opening
if($mtf_file ne '')
{
	open MOTIF, $mtf_file or die "Error: could not open specified Motifs file: $mtf_file\n";
	while(<MOTIF>)
	{
		chomp $_;
		if($_ =~ /^>/)
		{
			$_ = substr($_, 1, length($_));
			push(@name,$_);push(@id,$_."F");push(@id,$_."R");push(@count,0);push(@count,0);
		}
		else
		{
			push(@find,$_);
			push(@find,revcomp($_));
		}
	}
	close MOTIF;
}
#FASTA file opening
open FASTA, $input or die "FASTA file error: $input is unespecified or empty\n";
my $seq = "";
my $seqid = "";
while(<FASTA>)
{
	chomp $_;
	if($_ =~ /^>/)
	{
		$_ = substr($_, 1, length($_));
		my @spl = split(' ', $_);
		$seqid = $spl[0];
		print "DNA sequence found:\n\tseqid: $seqid\n\tsource: $source\n\n";
	}
	else
	{
		$seq .= $_;
	}
}
close FASTA;
my $seq_len = length($seq);
print "Sequence was loaded into memory\nCalculating GC, ATGC and N content\n";
#GC content calculation
my $gc_content = `grep -v ">" ../cbeijerinckii/cbeijerinckii1.fasta | grep -E "G|C" -o | wc -l`;
chomp $gc_content;
my $atcg_content = `grep -v ">" ../cbeijerinckii/cbeijerinckii1.fasta | grep -E "A|T|G|C" -o | wc -l`;
chomp $atcg_content;
#Búsqueda de patrones
if($mtf_file ne '')
{
print "Looking for specified motifs\n";
open OUTPUT, ">$fname.gff3" or die "Error: Gff file could not be generated\n";
print OUTPUT "##gff-version 3\n";
my $find_len = @find;
	for (my $j = 0; $j < $find_len; $j++)
	{
		while ($seq =~ m/$find[$j]/ig)
		{
			my $end = $seq_len - length($');
			my $start = length($`) + 1;
			$count[$j]++;
			my $curr_chain = $chain[$j%2];
			my $curr_name = $name[$j/2];
			print OUTPUT "$seqid\tDLPerl\t$curr_name\t$start\t$end\t.\t$curr_chain\t.\tID=$id[$j]$count[$j]\n";
		}
	}
close OUTPUT;
print "Gff file generated as $fname.gff3\n";
}
open OUTPUT_SUMMARY, ">$fname\_summary.tsv" or die "Summary file could not be created\n";
my $other_content = $seq_len - $atcg_content;
my $gc_percentage  = ($gc_content/$seq_len)*100;
print OUTPUT_SUMMARY "#ID\t\%G+C\tN content\tTotal\t".join("\t",@id)."\n";
print OUTPUT_SUMMARY "$seqid\t$gc_percentage\t$other_content\t$seq_len\t".join("\t",@count)."\n";
close OUTPUT_SUMMARY;
print "$fname\_summary.tsv file was created\nBye!\n";
sub revcomp
{
	my $rc = reverse $_;
	$rc =~ tr/({})[]ATGCatgc/)}{(][TACGtacg/;
	return $rc;
}
exit;
