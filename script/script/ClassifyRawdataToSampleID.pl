#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::ArgParse;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Switch;
use JSON;
use YAML::XS;

my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
 
my $arguments = Getopt::ArgParse->new_parser(
		description => 'This is a program,atuo link rawdata with its real sampleID',
		epilog      => 'Support samples with two or more sequence data, which will exist in result txt with multiple lines ',
);

$arguments -> add_arg('--rawdir', '-r', required => 1, type => "Scalar", help => "directory store raw sequences data");
$arguments -> add_arg('--maplist', '-m', required => 1, type => "Scalar", help => "text file contain sampleID LibraryID,split with \\t");
$arguments -> add_arg('--out', '-o', required => 1, type => "Scalar", help => "directory to store classify result , output will stored in out/raw_fastq.txt");
$arguments -> add_arg('--outfmt', '-f', required => 0,choices => [ 'text', 'json','yaml' ], type => "Scalar", help => "output file format  , defualt is text , can change to json or yaml",default=>"text");

$arguments -> print_usage() if( @ARGV < 1); 

my $opt=$arguments -> parse_args();

mkdir $opt->out unless(-d $opt->out);

#######################################################################################
# use linux command find to list all fastq files 
# fastq files should end with gz , and  must be paired R1 R2 
#######################################################################################

my $rawdir=$opt->rawdir;
my @fastqlist=sort `find $rawdir -name "*gz"|xargs -i realpath \{\} `;

my @left=map{$fastqlist[$_]}grep { $_%2==0}(0..$#fastqlist);
my @right=map{$fastqlist[$_]}grep {$_%2==1}(0..$#fastqlist);
#######################################################################################
# check fastq files.  
# the file name of paired R1 R2 has only one different char  ,
# which is 1 in R1 and 2 in R2 , and this case R2 after R1 when use sort command 
# so I sort the fastq list ,check if eary two file has same length, and only one char difference
#######################################################################################
compareR1R2(\@left,\@right);


my %idmap;
my %library;
my %cleanData;
#maplist file can be yaml format
if ($opt->maplist =~/yaml/||$opt->maplist =~/yml/) {
	my $yaml=YAML::XS::LoadFile($opt->maplist);
	%idmap=%{$yaml->{"library2sample"}};
}else{
	open MAP,$opt->maplist or die $!;
	while (<MAP>) {
		chomp $_;
		next if $_!~/\S/;
		my @lines=split "\t",$_;
		$idmap{$lines[1]}=$lines[0];
		$library{$lines[0]}=$lines[1];
		$cleanData{$lines[0]}=$lines[2] if scalar(@lines) > 2;
	}
	close MAP;

}

if (scalar(keys %idmap) != scalar(keys %library)) {
	die "numbers of samples and library not match, make sure its one to one mapping !!\n";
}

#######################################################################################
# link fastq files to sampleID  
# for some libraryID is to simple, like LRA0001 and LRA0001-1 .if we match LRA0001 first ,LRA0001-1 will also be matched
# So first sort the sampleID by the reversed order of libraryID length ,than match the fastq file to libraryID one by one 
#######################################################################################

my %sample2fastq;

foreach my $lib (sort {length($idmap{$b}) <=> length($idmap{$a})} keys %idmap) {
	my $num=$#left;
	foreach my $i (0..$num) {
		my $read1=shift @left;
		my $read2=shift @right;
		chomp ($read1,$read2);
		if ($read1=~/$lib/) {
			push @{$sample2fastq{$idmap{$lib}}{"R1"}},$read1;
			push @{$sample2fastq{$idmap{$lib}}{"R2"}},$read2;
		}else{
			push @left,$read1;
			push @right,$read2;
		}
	}
	if (!defined $sample2fastq{$idmap{$lib}}) {
		print "No fastq files matched $lib for sample:",$idmap{$lib},"\n";
	}
}

if (@left) {
	print "Some fastq file didnt mapped to samples , you should check the ids , make sure all sample is given, no one left\n ";
	print join("\n",@left),"\n";
}

#######################################################################################
# out put match result
# I out result in text format by defualt ,and also can be json or yaml format ,for other developments
#######################################################################################


switch($opt->outfmt){                             
	case "text" {
		open OUT,">".$opt->out."/raw_fastq.txt";
		foreach my $sample (sort keys %sample2fastq) {
			for (my $i=0;$i < @{$sample2fastq{$sample}{"R1"}} ;$i++) {
				my @tables=($sample,$library{$sample},$sample2fastq{$sample}{"R1"}[$i],$sample2fastq{$sample}{"R2"}[$i]);
				push @tables,$cleanData{$sample} if(defined $cleanData{$sample});
				print OUT join("\t",@tables),"\n";
			}
		}
		close OUT;
	}        
	case "json" {
		open OUT,">".$opt->out."/raw_fastq.json";
		my $json_txt=encode_json(\%sample2fastq);
		print OUT $json_txt;
		close OUT;
	}            
	case "yaml" {
		open OUT,">".$opt->out."/raw_fastq.yaml";
		my $yaml_txt=YAML::XS::Dump(\%sample2fastq);
		print OUT $yaml_txt;
		close OUT;
	}            
}                                             

################################################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

sub compareR1R2{
	my @left=@{$_[0]};
	my @right=@{$_[1]};
	

	if (scalar(@left) != scalar(@right)) {
		print "Number of R1 and R2 is no equal , please check!!\nNumber of R1 :",scalar(@left),"\nNumber of R2 :",scalar(@right),"\n";
		die;
	}
	my $erro=0;	
	for (my $i=0;$i<@left ;$i++) {
		my $str1=$left[$i];
		my $str2=$right[$i];

		chomp($str1,$str2);
		
		my %diff=FindDiffernceBettweenString($str1,$str2);
		$erro++ if (scalar(keys %diff) > 1) ;
		print "As defualt set, the difference bettween filenames of R1 and R2 is only one chart, which is 1 and 2 \n",ShownDifference($str1,$str2,\%diff),"\n" if(scalar(keys %diff) > 1) ;
	}

	die "Error when compare R1 and R2,please check!!\n" if $erro; 

}

sub FindDiffernceBettweenString{
	my @chars1=split "",$_[0];
	my @chars2=split "",$_[1];
	
	die "strings with different length\n",join("\n",@_),"\n" if(scalar(@chars1) != scalar(@chars2));

	my %difference;
	
	for (my $i=0;$i<@chars1;$i++) {
		my $char1=$chars1[$i];
		my $char2=$chars2[$i];
		next if ($char1 eq  $char2);
		$difference{$i}=$char1." AND ".$char2;
	}
	return %difference;


}


sub ShownDifference{
	my $str=$_[0]."\n".$_[1]."\n";
	foreach my $i (0..length($_[0])-1) {
		my $char= defined $_[2]->{$i} ? "*":"-";
		$str.=$char;
	}
	return $str;
}


# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

