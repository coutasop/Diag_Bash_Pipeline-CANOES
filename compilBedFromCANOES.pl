use Getopt::Long;
use strict;

#Define boolean constant
use constant false => 0;
use constant true  => 1;

my $in;
my $out;

GetOptions(
	"in=s"  => \$in,
	"out=s" => \$out
);

open( INFILE,  "<$in" )  or die("can't open input file : $in\n");
open( OUTFILE, ">$out" ) or die("can't open output file : $out\n");
print OUTFILE 'Chr', "\t", 'Start', "\t", 'End', "\t", 'Sample', "\t", 'CNVType', "\t", 'Kb', "\t", 'CopyNbr', "\t", 'QSome', "\t", 'Gene',"\n";
my $lastGene;
my $fLine = <INFILE>;
my %storedline;
chomp($fLine);
my @fLines = split( "\t", $fLine );
$storedline{'chr'}   = $fLines[0];
$storedline{'start'} = $fLines[1];
$storedline{'end'}   = $fLines[2];
$storedline{'ind'}    = $fLines[3];
$storedline{'dd'}    = $fLines[4];
$storedline{'kb'}    = $fLines[5];
$storedline{'numTarg'}    = $fLines[6];
$storedline{'MLCN'}    = $fLines[7];
$storedline{'QSome'}    = $fLines[8];
$storedline{'gene'} = $fLines[13];
$lastGene = $fLines[13];

while ( my $line = <INFILE> ) {
	chomp($line);
	my @lines = split( "\t", $line );
	my $diff = $lines[1] - $storedline{'end'};
	if (   ( $lines[0] eq $storedline{'chr'} )
		&& ( $diff <= 10000 )
		&& ( $lines[4] eq $storedline{'dd'} )
		&& ( $lines[3] eq $storedline{'ind'} )
	 )
	{
		$storedline{'end'}=$lines[2];
		if ($lastGene ne $lines[13] ){
			$storedline{'gene'}="$storedline{'gene'},$lines[13]";
			$lastGene = $lines[13];
		}
	}
	else {
		my @list = split(",",$storedline{'gene'});
		my %genes;
		foreach my $gene (@list){
			$genes{$gene}=1;
		}
		#~ print OUTFILE $storedline{'chr'}, "\t", $storedline{'start'}, "\t",
		  #~ $storedline{'end'}, "\t", $storedline{'ind'}, "\t", $storedline{'dd'}, "\t",join(",",keys(%genes)),
		  #~ "\n";
		print OUTFILE $storedline{'chr'}, "\t", $storedline{'start'}, "\t",	$storedline{'end'}, "\t", $storedline{'ind'}, "\t", $storedline{'dd'}, "\t", $storedline{'kb'}, "\t", $storedline{'MLCN'}, "\t", $storedline{'QSome'}, "\t", join(",",keys(%genes)),"\n";
		$storedline{'chr'}   = $lines[0];
		$storedline{'start'} = $lines[1];
		$storedline{'end'}   = $lines[2];
		$storedline{'ind'}    = $lines[3];
		$storedline{'dd'}    = $lines[4];
		$storedline{'kb'}    = $lines[5];
		$storedline{'numTarg'}    = $lines[6];
		$storedline{'MLCN'}    = $lines[7];
		$storedline{'QSome'}    = $lines[8];
		$storedline{'gene'} = $lines[13];
		$lastGene = $lines[13];
	}
}
my @list = split(",",$storedline{'gene'});
my %genes;
foreach my $gene (@list){
	$genes{$gene}=1;
}
#~ print OUTFILE $storedline{'chr'}, "\t", $storedline{'start'}, "\t",
	#~ $storedline{'end'}, "\t", $storedline{'ind'}, "\t", $storedline{'dd'}, "\t",join(",",keys(%genes)),
        #~ "\n";
print OUTFILE $storedline{'chr'}, "\t", $storedline{'start'}, "\t",	$storedline{'end'}, "\t", $storedline{'ind'}, "\t", $storedline{'dd'}, "\t", $storedline{'kb'}, "\t", $storedline{'MLCN'}, "\t", $storedline{'QSome'}, "\t", join(",",keys(%genes)),"\n";

close(INFILE);
close(OUTFILE);
