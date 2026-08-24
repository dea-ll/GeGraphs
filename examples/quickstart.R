library(GeGraphs)

sample_bam <- "path/to/sample.bam"
controls <- c(
  Control_1 = "path/to/control_1.bam",
  Control_2 = "path/to/control_2.bam"
)
sv_files <- list(
  Delly = "path/to/sample.delly.vcf",
  Dysgu = "path/to/sample.dysgu.vcf"
)

plot_gene(
  gene = "ERGIC1",
  bam_file = sample_bam,
  control_bams = controls,
  vcf_list = sv_files,
  genome = "hg19",
  flank = 5000,
  show_gene = TRUE,
  show_cov = TRUE,
  show_sv = TRUE,
  bin_size = 2000,
  out = "ERGIC1_review.pdf"
)

plot_region(
  chrom = "5",
  start = 172260000,
  end = 172360000,
  bam_file = sample_bam,
  control_bams = controls,
  vcf_list = sv_files,
  genome = "hg19",
  bin_size = 10000,
  out = "chr5_review.pdf"
)

genes_under <- get_genes_under_sv(
  bed_file = "path/to/sample.filtered_DEL.bed",
  genome = "hg19"
)

context <- get_gene_context_for_gene(
  gene = "ERGIC1",
  genome = "hg19",
  n_upstream = 2,
  n_downstream = 2
)
