#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

=pod

=head1 Usage
	perl $0 [option] <bam> <outdir>
		-q	base quality [default 0]
		-Q	mapping quality [default 0]
		-s	sample name[default 'sample']
		-r 	region file(bed format, just for WES. WGS need not this file)
		-h	help		
		# Any bugs please report to zhangchao@novogene.cn 
=cut
		

my ($basethres,$mapQthres,$total_chr,$bedfile,$sample_name,$addnFile,$help);
GetOptions("q:i"=>\$basethres,"Q:i"=>\$mapQthres,"r:s"=>\$bedfile,"s:s"=>\$sample_name,"n:s"=>\$addnFile,"h"=>\$help);

#$total_chr ||= 2471805657;
$sample_name ||= 'sample';
$basethres ||= 0;
$mapQthres ||= 0;
#$addnFile="/PUBLIC/database/HEALTH/genome/human/b37_gatk/all.NBlock.larger1000bp.bed";
$addnFile="Nblock.bed";
die `pod2text $0` if(@ARGV<2 || $help);


my $bam=shift;
my $outdir=shift;

####################### init
my %depth=();
my $maxCov=0;
my $Average_sequencing_depth=0;
my $Average_sequencing_depth4=0;
my $Average_sequencing_depth10=0;
my $Average_sequencing_depth20=0;
my $Coverage=0;
my $Coverage4=0;
my $Coverage10=0;
my $Coverage20=0;
my $Coverage_bases=0;
my $Coverage_bases_4=0;
my $Coverage_bases_10=0;
my $Coverage_bases_20=0;

my $total_Coverage_bases=0;
my $total_Coverage_bases_4=0;
my $total_Coverage_bases_10=0;
my $total_Coverage_bases_20=0;



######################## end init

my %bychr_depth = ();
my %bychr_coverage = ();
my %length_chr_WES = ();

##calcualte genome length
my %length_chr = ();
open CHR,"samtools view -H $bam | " or die;
while (<CHR>)
{
	chomp;
	if($_ =~ /SQ\s+SN\:([\d|X|Y]+)\s+LN\:(\d+)/)
	{
		$length_chr{$1} = $2;
		
	}
	if($_ =~ /LN\:(\d+)/)
	{
		$total_chr += $1;
	}

}
close CHR;

my %length_nblock = ();

if (!defined($bedfile))
{

	open NBLOCK,"$addnFile" or die;
	while (<NBLOCK>)
	{
		my @words = split;
		$length_nblock{$words[0]} = $words[3];
		$total_chr -= $words[3];
	}
	close NBLOCK;


	`mkdir -p $outdir` unless -d $outdir;

	
	open DEPTH,"samtools depth -q $basethres -Q $mapQthres $bam | " or die;
	while(<DEPTH>)
	{
		chomp;
		my @arr = split;
		next if ($arr[2] == 0);
		$depth{$arr[2]}+=1;

		#for coverage bychr
		$bychr_depth{$arr[0]} += $arr[2];
		$bychr_coverage{$arr[0]} += 1;
	}
	close DEPTH;


	my @depth=sort {$a<=>$b} keys %depth;
	my %counts = ();

	open HIS,">$outdir/depth_frequency.xls" or die;
	open CUM,">$outdir/cumu.xls" or die;
	open CUN,">$outdir/$sample_name\.sample_cumulative_coverage_counts" or die;
	print CUM "Depth\tPercent\n0\t1\n";

	for my $i (0..500)
	{
		print CUN "\tget_$i";
	}
	print CUN "\n$sample_name";
	$counts{"0"} = $total_chr;

	foreach my $depth1 (@depth)
	{
		
		my $per=$depth{$depth1}/$total_chr;
		$total_Coverage_bases += $depth1*$depth{$depth1};
		$Coverage_bases += $depth{$depth1};

		if($depth1>=4)	
		{
			$total_Coverage_bases_4 += $depth1 * $depth{$depth1};
			$Coverage_bases_4 += $depth{$depth1};
		}
		if($depth1>=10)
		{
			$total_Coverage_bases_10 += $depth1 * $depth{$depth1};
			$Coverage_bases_10 += $depth{$depth1};
		}
		if($depth1>=20)
		{
			$total_Coverage_bases_20 += $depth1 * $depth{$depth1};
			$Coverage_bases_20 += $depth{$depth1};
		}

		
		$maxCov=$per if($maxCov<$per);
		my $tmp=0;
		print HIS "$depth1\t$per\t$depth{$depth1}\n";
		foreach my $depth2(@depth)
		{
			$tmp+=$depth{$depth2} if($depth2 >= $depth1);
		}
		$counts{$depth1} = $tmp;
		$tmp=$tmp/$total_chr;
		print CUM "$depth1\t$tmp\n";


	}

	for my $i (0..500)
	{
		print CUN "\t$counts{$i}";
	}
	close HIS;
	close CUM;
	close CUN;


	$Average_sequencing_depth=$total_Coverage_bases/$total_chr;
	$Coverage=$Coverage_bases/$total_chr;
	$Average_sequencing_depth4=$total_Coverage_bases_4/$total_chr;
	$Coverage4=$Coverage_bases_4/$total_chr;
	$Average_sequencing_depth10=$total_Coverage_bases_10/$total_chr;
	$Coverage10=$Coverage_bases_10/$total_chr;
	$Average_sequencing_depth20=$total_Coverage_bases_20/$total_chr;
	$Coverage20=$Coverage_bases_20/$total_chr;

	open STAT,">$outdir/information.xlsx" or die $!;
	print STAT "Average_sequencing_depth:\t",sprintf("%.2f",$Average_sequencing_depth),"\n";
	print STAT "Coverage:\t",sprintf("%.2f%%",100*$Coverage),"\n";
	print STAT "Coverage_at_least_4X:\t",sprintf("%.2f%%",100*$Coverage4),"\n";
	print STAT "Coverage_at_least_10X:\t",sprintf("%.2f%%",100*$Coverage10),"\n";
	print STAT "Coverage_at_least_20X:\t",sprintf("%.2f%%",100*$Coverage20),"\n";
	close STAT;

}



###for WES
my $Initial_bases_on_target=0;
my $Initial_bases_near_target=0;
my $Initial_bases_on_or_near_target=0;
my $Total_effective_reads=0;
my $Total_effective_yield=0;   # = $total_Coverage_bases
my $Average_read_length=0;
my $Effective_sequences_on_target=0;
my $Effective_sequences_near_target=0;
my $Effective_sequences_on_or_near_target=0;
my $Fraction_of_effective_bases_on_target=0;
my $Fraction_of_effective_bases_on_or_near_target=0;
my $Average_sequencing_depth_on_target=0;
my $Average_sequencing_depth_near_target=0;
my $Mismatch_rate_in_target_region=0;
my $Mismatch_rate_in_all_effective_sequence=0;
my $Base_covered_on_target=0;
my $Coverage_of_target_region=0;
my $Base_covered_near_target=0;
my $Coverage_of_flanking_region=0;
my $Fraction_of_target_covered_with_at_least_20x=0;
my $Fraction_of_target_covered_with_at_least_10x=0;
my $Fraction_of_target_covered_with_at_least_4x=0;
my $Fraction_of_flanking_region_covered_with_at_least_20x=0;
my $Fraction_of_flanking_region_covered_with_at_least_10x=0;
my $Fraction_of_flanking_region_covered_with_at_least_4x=0;	


##init
my %length_flank_chr_WES = ();
my $Mismatch_base_in_all_effective_sequence = 0;
my $Mismatch_base_in_target_region = 0;


if (defined($bedfile))
{
	

	open CHRWES,"$bedfile" or die;
	open TMP,">$outdir/tmp.region.bed";
	while(<CHRWES>)
	{
		chomp;
		my @info = split;
		$length_chr_WES{$info[0]} += $info[2] - $info[1];

		$Initial_bases_on_target += $info[2]-$info[1];

		my $pri_beg_pos=$info[1]-200;
		my $pri_end_pos=$info[1];
		my $next_beg_pos=$info[2];
		my $next_end_pos=$info[2]+200;
		if($pri_beg_pos<0) { $pri_beg_pos=0; }
		
		print TMP "$info[0]\t$pri_beg_pos\t$pri_end_pos\n";
		print TMP "$info[0]\t$next_beg_pos\t$next_end_pos\n";
	}
	close TMP;
	close CHRWES;


	`sortBed -i $outdir/tmp.region.bed | mergeBed -i stdin | subtractBed -a stdin -b $bedfile > $outdir/flank_region.bed`;
	`cat $outdir/flank_region.bed $bedfile | sortBed -i stdin | mergeBed -i stdin > $outdir/flank_target_region.bed`;
	
	`samtools depth -q $basethres -Q $mapQthres -b $outdir/flank_region.bed $bam > $outdir/flank_region.depth`;
	`samtools depth -q $basethres -Q $mapQthres -b $bedfile $bam > $outdir/target_region.depth`;

	$Total_effective_yield=`samtools depth -q $basethres -Q $mapQthres  $bam | awk '{total+=\$3};END{print total}'`;
	
	open FLANK,"$outdir/flank_region.bed" or die;
    while(<FLANK>)
    {
        chomp;
        my @info = split;
        #$length_flank_chr_WES{$info[0]} += $info[2] - $info[1];

        $Initial_bases_near_target += $info[2] - $info[1];
    }
    close FLANK;

	%bychr_depth = ();
	%bychr_coverage = ();

	open DEPTH,"$outdir/target_region.depth" or die;
	while(<DEPTH>)
	{
		chomp;
		my @arr = split;
		next if ($arr[2] == 0);

		#for coverage bychr
		$bychr_depth{$arr[0]} += $arr[2];
		$bychr_coverage{$arr[0]} += 1;
	}
	close DEPTH;
		
	$Initial_bases_on_or_near_target = $Initial_bases_on_target + $Initial_bases_near_target;

	my $tmp1=`awk '{total4++};{total+=\$3};\$3 >=20 {total1++};\$3 >=10 {total2++};\$3 >=4 {total3++};END{print total"\t"total1"\t"total2"\t"total3"\t"total4}' $outdir/target_region.depth`;
	chomp($tmp1);
	my @info1;
	@info1 = split /\t/, $tmp1;
	
	$Effective_sequences_on_target = $info1[0];
	if(defined($info1[1]) or $info1[1]=0 )
	{
		$Fraction_of_target_covered_with_at_least_20x = $info1[1]/$Initial_bases_on_target;
	}
	if(defined($info1[2]) or $info1[2]=0 )
	{
		$Fraction_of_target_covered_with_at_least_10x = $info1[2]/$Initial_bases_on_target;
	}
	if(defined($info1[3]) or $info1[3]=0 )
	{
		$Fraction_of_target_covered_with_at_least_4x = $info1[3]/$Initial_bases_on_target;
	}
	
	
	$Base_covered_on_target = $info1[4];
	

	#$Effective_sequences_near_target=`awk '{total+=\$3};END{print total}' $outdir/flank_region.depth`;
	#chomp($Effective_sequences_near_target);
	my $tmp2 = `awk '{total4++};{total+=\$3};\$3 >=20 {total1++};\$3 >=10 {total2++};\$3 >=4 {total3++};END{print total"\t"total1"\t"total2"\t"total3"\t"total4}' $outdir/flank_region.depth`;
	chomp($tmp2);
	my @info2;
	@info2 = split /\t/, $tmp2;

	$Effective_sequences_near_target = $info2[0];
	if(defined($info2[1]) or $info2[1]=0 )
	{
		$Fraction_of_flanking_region_covered_with_at_least_20x=$info2[1]/$Initial_bases_near_target;
	}
	if(defined($info2[2]) or $info2[2]=0 )
	{
		$Fraction_of_flanking_region_covered_with_at_least_10x=$info2[2]/$Initial_bases_near_target;
	}
	if(defined($info2[3]) or $info2[3]=0 )
	{
		$Fraction_of_flanking_region_covered_with_at_least_4x=$info2[3]/$Initial_bases_near_target;
	}

	$Base_covered_near_target = $info2[4];

	$Effective_sequences_on_or_near_target = $Effective_sequences_on_target + $Effective_sequences_near_target;



	#Total_effective_reads
	open BAM,"samtools view -F 0x0004 $bam | " or die $!;
	while(<BAM>)
	{
		chomp;
		my @_F=split /\t/;
		#if($_F[1]=~/d/) { next; }
		$Total_effective_reads++;
		if($_ =~ /MD:Z:([\d|A|T|C|G|N]+)\t/)
		{
			my $count=()=$1=~/\d[A|T|C|G|N]/g;
			$Mismatch_base_in_all_effective_sequence += $count;
    	}
	}
	close BAM;

	open BAM,"samtools view -F 0x0004 -L $bedfile $bam | " or die $!;
	while(<BAM>)
	{
		chomp;
		my @_F=split /\t/;
		if($_ =~ /MD:Z:([\d|A|T|C|G|N]+)\t/)
		{
			my $count=()=$1=~/\d[A|T|C|G|N]/g;
			$Mismatch_base_in_target_region += $count;
    	}
	}
	close BAM;

	#
	
	$Average_read_length = $Total_effective_yield/$Total_effective_reads;

	#Fraction_of_effective_bases_on_target
	$Fraction_of_effective_bases_on_target = $Effective_sequences_on_target/$Total_effective_yield;

	#Fraction_of_effective_bases_on_or_near_target
	$Fraction_of_effective_bases_on_or_near_target = $Effective_sequences_on_or_near_target/$Total_effective_yield;

	#Average_sequencing_depth_on_target
	$Average_sequencing_depth_on_target = $Effective_sequences_on_target/$Initial_bases_on_target;

	#Average_sequencing_depth_near_target
	$Average_sequencing_depth_near_target = $Effective_sequences_near_target/$Initial_bases_near_target;

	#Mismatch_rate_in_target_region
	$Mismatch_rate_in_target_region=$Mismatch_base_in_target_region/$Effective_sequences_on_target;

	#Mismatch_rate_in_all_effective_sequence
	$Mismatch_rate_in_all_effective_sequence=$Mismatch_base_in_all_effective_sequence/$Total_effective_yield;

	#Coverage_of_target_region
	$Coverage_of_target_region=$Base_covered_on_target/$Initial_bases_on_target;

	#Coverage_of_flanking_region
	$Coverage_of_flanking_region=$Base_covered_near_target/$Initial_bases_near_target;




	#output
	open STAT,">$outdir/information.xlsx" or die $!;
	print STAT "Initial_bases_on_target:\t$Initial_bases_on_target\n";
	print STAT "Initial_bases_near_target:\t$Initial_bases_near_target\n";
	print STAT "Initial_bases_on_or_near_target:\t$Initial_bases_on_or_near_target\n";
	print STAT "Total_effective_reads:\t$Total_effective_reads\n";
	printf STAT "Total_effective_yield(Mb):\t%.2f\n",$Total_effective_yield/1000000;
	printf STAT "Average_read_length(bp):\t%.2f\n",$Average_read_length;
	printf STAT "Effective_sequences_on_target(Mb):\t%.2f\n",$Effective_sequences_on_target/1000000;
	printf STAT "Effective_sequences_near_target(Mb):\t%.2f\n",$Effective_sequences_near_target/1000000;
	printf STAT "Effective_sequences_on_or_near_target(Mb):\t%.2f\n",$Effective_sequences_on_or_near_target/1000000;
	printf STAT "Fraction_of_effective_bases_on_target:\t%.1f%%\n",100*$Fraction_of_effective_bases_on_target;
	printf STAT "Fraction_of_effective_bases_on_or_near_target:\t%.1f%%\n",100*$Fraction_of_effective_bases_on_or_near_target;
	printf STAT "Average_sequencing_depth_on_target:\t%.2f\n",$Average_sequencing_depth_on_target;
	printf STAT "Average_sequencing_depth_near_target:\t%.2f\n",$Average_sequencing_depth_near_target;
	printf STAT "Mismatch_rate_in_target_region:\t%.2f%%\n",100*$Mismatch_rate_in_target_region;
	printf STAT "Mismatch_rate_in_all_effective_sequence:\t%.2f%%\n",100*$Mismatch_rate_in_all_effective_sequence;
	print STAT "Base_covered_on_target:\t$Base_covered_on_target\n";
	printf STAT "Coverage_of_target_region:\t%.1f%%\n",100*$Coverage_of_target_region;
	print STAT "Base_covered_near_target:\t$Base_covered_near_target\n";
	printf STAT "Coverage_of_flanking_region:\t%.1f%%\n",100*$Coverage_of_flanking_region;
	printf STAT "Fraction_of_target_covered_with_at_least_20x:\t%.1f%%\n",100*$Fraction_of_target_covered_with_at_least_20x;
	printf STAT "Fraction_of_target_covered_with_at_least_10x:\t%.1f%%\n",100*$Fraction_of_target_covered_with_at_least_10x;
	printf STAT "Fraction_of_target_covered_with_at_least_4x:\t%.1f%%\n",100*$Fraction_of_target_covered_with_at_least_4x;
	printf STAT "Fraction_of_flanking_region_covered_with_at_least_20x:\t%.1f%%\n",100*$Fraction_of_flanking_region_covered_with_at_least_20x;
	printf STAT "Fraction_of_flanking_region_covered_with_at_least_10x:\t%.1f%%\n",100*$Fraction_of_flanking_region_covered_with_at_least_10x;
	printf STAT "Fraction_of_flanking_region_covered_with_at_least_4x:\t%.1f%%\n",100*$Fraction_of_flanking_region_covered_with_at_least_4x;
	close STAT;
	
	
	my %depth = ();
	my @depth = ();
	open DEPTH,"$outdir/target_region.depth" or die;
	while(<DEPTH>)
	{
		chomp;
		my @arr = split;
		next if ($arr[2] == 0);
		$depth{$arr[2]}+=1;		
	}
	@depth=sort {$a<=>$b} keys %depth;
	my %counts = ();
	
	open HIS,">$outdir/depth_frequency.xls" or die;
	open CUM,">$outdir/cumu.xls" or die;
	open CUN,">$outdir/$sample_name\.sample_cumulative_coverage_counts" or die;
	print CUM "Depth\tPercent\n0\t1\n";
	for my $i (0..500)
	{
		print CUN "\tget_$i";
	}
	print CUN "\n$sample_name";
	$counts{"0"} = $Initial_bases_on_target;


	foreach my $depth1 (@depth)
	{
		my $per=$depth{$depth1}/$Initial_bases_on_target;

		$maxCov=$per if($maxCov<$per);
		my $tmp=0;
		print HIS "$depth1\t$per\t$depth{$depth1}\n";
		foreach my $depth2(@depth)
		{
			$tmp+=$depth{$depth2} if($depth2 >= $depth1); 
		}
		$counts{$depth1} = $tmp;
		$tmp=$tmp/$Initial_bases_on_target;
		print CUM "$depth1\t$tmp\n";
	}
	
	for my $i (0..500)
	{
		print CUN "\t$counts{$i}";
	}

	close CUN;
	close HIS;
	close CUM;

	#remove
	`rm $outdir/flank_region.depth $outdir/target_region.depth $outdir/tmp.region.bed $outdir/flank_region.bed $outdir/flank_target_region.bed`;
}

open BYCHR,">$outdir/$sample_name".".coverage.bychr.txt" or die;
print BYCHR "chr\tlength\ttotal_depth\tmean_depth\tcovered_bases\tprop_covered_bases";
for my $i (1..22,'X','Y')
{
	my $len_ral = 0;
	if (defined($bedfile))
	{	
		$len_ral = $length_chr_WES{$i};
	}
	else
	{
		$len_ral = $length_chr{$i} - $length_nblock{$i};
	}
	my $aver_depth = sprintf("%.2f", $bychr_depth{$i}/$len_ral);
	my $aver_coverage = sprintf("%.3f", $bychr_coverage{$i}/$len_ral);
	print BYCHR "\n$i\t$len_ral\t$bychr_depth{$i}\t$aver_depth\t$bychr_coverage{$i}\t$aver_coverage";
}
close BYCHR;


if(1)
{
	my $ylim = 100*$maxCov;
	my ($xbin,$ybin);
	$ylim= int($ylim) + 1;
	if($ylim <= 3)
	{
		$ybin = 0.5;
	}else{
		$ybin=1;
	}
	my $xlim=0;
	if($Average_sequencing_depth<30)
	{
		$xlim=100;
		$xbin=20;
	}elsif($Average_sequencing_depth < 50)
	{
		$xlim=160;
		$xbin=20;
	}elsif($Average_sequencing_depth  < 120)
	{
		$xlim=250;
		$xbin=50;
	}else{
		$xlim=600;
		$xbin=100;
	}
	histPlot($outdir,"$outdir/depth_frequency.txt",$ylim,$ybin,$xlim,$xbin);
	cumuPlot($outdir,"$outdir/cumu.txt",$xlim,$xbin);
}

sub cumuPlot {
	my ($outdir, $dataFile, $xlim, $xbin) = @_;
	my $figFile = "$outdir/cumuPlot.pdf";
	my $Rline=<<Rline;
	pdf(file="$figFile",w=8,h=6)
	rt <- read.table("$dataFile")
	opar <- par()
	x <- rt\$V1[1:($xlim+1)]
	y <- 100*rt\$V2[1:($xlim+1)]
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col="red",type='l', lwd=2, bty="l",xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0, 100))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,100,by=20)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Cumulative sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1.4)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.4)
	par(opar)
	dev.off()
	png(filename="$outdir/cumuPlot.png",width = 480, height = 360)
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col="red",type='l', lwd=3, bty="l",xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0, 100))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,100,by=20)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Cumulative sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1.5)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.5)
	par(opar)
	dev.off()
	
Rline
	open (ROUT,">$figFile.R");
	print ROUT $Rline;
	close(ROUT);

	system("R CMD BATCH  $figFile.R");
}


sub histPlot {
	my ($outdir, $dataFile, $ylim, $ybin, $xlim, $xbin) = @_;
	my $figFile = "$outdir/histPlot.pdf";
	my $Rline=<<Rline;
	pdf(file="$figFile",w=8,h=6)
	rt <- read.table("$dataFile")
	opar <- par()
	t=sum(rt\$V2[($xlim+1):length(rt\$V2)])
	y=c(rt\$V2[1:$xlim],t)
	y <- y*100
	x <- rt\$V1[1:($xlim+1)]
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col="blue",type='h', lwd=1.5, xaxt="n",yaxt="n", xlab="", ylab="", bty="l",ylim=c(0,$ylim),xlim=c(0,$xlim))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,$ylim,by=$ybin)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	end <- length(xpos)-1
	mtext(c(xpos[1:end],"$xlim+"), side=1, las=1, at=xpos, line=0.3, cex=1.4)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.4)
	par(opar)
	dev.off()
	png(filename="$outdir/histPlot.png",width = 480, height = 360)
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col="blue",type='h', lwd=1.5, xaxt="n",yaxt="n", xlab="", ylab="", bty="l",ylim=c(0,$ylim),xlim=c(0,$xlim))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,$ylim,by=$ybin)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	end <- length(xpos)-1
	mtext(c(xpos[1:end],"$xlim+"), side=1, las=1, at=xpos, line=0.3, cex=1.5)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.5)
	par(opar)
	dev.off()
Rline
	open (ROUT,">$figFile.R");
	print ROUT $Rline;
	close(ROUT);

	system("R CMD BATCH  $figFile.R");
	system("rm  $figFile.R  $figFile.Rout");
}
