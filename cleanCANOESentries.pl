use Getopt::Long;
use strict;
use List::Util qw[min max];

#Define boolean constant
use constant false => 0;
use constant true  => 1;

my $gc;
my $reads;
my $outPrefix;

GetOptions(
	"gc=s"        => \$gc,
	"reads=s"     => \$reads,
	"outPrefix=s" => \$outPrefix
);

open( GCFILE, "<$gc" ) or die("can't open gc file : $gc\n");
open( READSFILE, "<$reads" ) or die("can't open reads file : $reads\n");

open( GCOUTFILE, ">$outPrefix.gc.txt" ) or die("can't open gc output file : $outPrefix.gc.txt\n");
open( READSOUTFILE, ">$outPrefix.reads.txt" ) or die("can't open reads output file : $outPrefix.reads.txt\n");

my %remove;
my %chr;
$chr{"1"}=[];
while(my $line = <READSFILE>){
	my $keep = true;
	chomp($line);
	my @lines = split("\t",$line);
	for (my $i=2;$i<scalar(@lines);$i++){
		if($lines[$i]<=0){	#seuil profondeur reads.txt
			$keep = false;
		}
	}
	if($keep){
		push(@{$chr{$lines[0]}},[@lines]);
	}else{
		my $newStart=$lines[1]+1;
		my $pos="$lines[0]:$newStart-$lines[2]";
		$remove{$pos}=1;
	}
}
while (my $line = <GCFILE>){
	my @lines = split("\t",$line);
	if (!exists($remove{$lines[0]})){
		print GCOUTFILE $line;
	}
}

for (my $i = 1; $i <=24; $i++){
	if(exists($chr{$i})){
		my @list = sort { $a->[1] <=> $b->[1]} @{$chr{$i}};
		for my $aref (@list){
			print READSOUTFILE join("\t",@$aref),"\n";	
		}
	}
}

if(exists($chr{"X"})){
	print "chrX found \n";
	my @list = sort { $a->[1] <=>  $b->[1]} @{$chr{"X"}};
	for my $aref (@list) {
		print READSOUTFILE join("\t",@$aref),"\n";
	}
}

if(exists($chr{"Y"})){
	my @list = sort { $a->[1] <=>  $b->[1]} @{$chr{"Y"}};
	for my $aref (@list) {
		print READSOUTFILE join("\t",@$aref),"\n";
	}
}

print "Output: $outPrefix.gc.txt\n";
print "Output: $outPrefix.reads.txt\n";
