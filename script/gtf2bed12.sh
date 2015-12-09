# gtfToGenePred were download form UCSC  : http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
gtfToGenePred -genePredExt -geneNameAsName2 genes.gtf gene.tmp
awk '{print $2"\t"$4"\t"$5"\t"$1"\t0\t"$3"\t"$6"\t"$7"\t0\t"$8"\t"$9"\t"$10}' genes.tmp >  genes.bed12
rm gene.tmp
