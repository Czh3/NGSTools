zcat $1 | vcf-annotate -a ~/reference/butterfly/Melitaea_cinxia_all_v1_reform.sort.gtf.gz -c CHROM,FROM,TO,INFO/GN -d key=INFO,ID=GN,Number=1,Type=String,Description='Gene Name' | bgzip > $2
