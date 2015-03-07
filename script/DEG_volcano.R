library(ggplot2)

#usage:
#	Rscript DEG_volcano.R input<cuff_diff/gene_exp.dff> out.png

#read data
args = commandArgs()
input_file = args[6]
output_file = args[7]
df = as.data.frame(read.table(input_file, header=TRUE, row.names=1))

#high expression genes
H_genes = c("RpL31", "RpS9", "Nacalpha", "regucalcin", "ATPsyn-gamma", "CG1746")
df.H_genes = df[which(df$gene %in% H_genes),]

#middle expression genes
M_genes = c("RanBPM", "CG6707", "CG5991", "nct", "CG6340")
df.M_genes = df[which(df$gene %in% M_genes), ]

#low expression genes
L_genes = c("Rcd2", "Aplip1", "qua", "Nep2", "CG6847")
df.L_genes = df[which(df$gene %in% L_genes), ]

p = ggplot() + geom_point(data = df, aes(x = log2.fold_change., y = -log10(p_value))) +
geom_point(data = df.H_genes, aes(x = log2.fold_change., y = -log10(p_value)), color='red') +
geom_point(data = df.M_genes, aes(x = log2.fold_change., y = -log10(p_value)), color='green') +
geom_point(data = df.L_genes, aes(x = log2.fold_change., y = -log10(p_value)), color='blue') +
xlim(-3, 3)

ggsave(file=output_file)
