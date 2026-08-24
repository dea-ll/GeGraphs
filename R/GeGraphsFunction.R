# GeGraphs: Gene-Centric Genomic Visualization Tools
#
# This file contains the core functions for the GeGraphs package,
# which provides gene-focused visualizations by integrating gene annotations,
# read coverage from BAM files, and structural variants from VCFs.


#' Plot CNV and Structural Variations for a Gene
#'
#' This function visualizes gene structure, coverage, and structural variations (SVs)
#' from BAM and VCF or BED files around a given gene, using MANE transcripts when desired.
#'
#' @param gene Gene symbol (HGNC) to plot (e.g., "NELL1"). Can be NULL if \code{hgnc_id} is provided.
#' @param hgnc_id HGNC identifier (e.g., "HGNC:7750"). If provided, it is used to determine the gene symbol.
#' @param bam_file BAM file with alignments.
#' @param control_bams A named list of control BAM file paths for normalization.
#' @param vcf_list A named list of VCF file paths for structural variants.
#' @param bed_list A named list of BED file paths for structural variants (alternative to \code{vcf_list}).
#' @param genome Genome assembly ("hg19" or "hg38").
#' @param out Output image filename (PNG).
#' @param chrom Chromosome (optional, inferred from the gene if NULL).
#' @param flank Number of base pairs to include upstream/downstream of the gene.
#' @param show_gene Logical, whether to show gene annotation.
#' @param show_cov Logical, whether to show read coverage.
#' @param show_sv Logical, whether to show structural variations.
#' @param sv_types Type(s) of SV to display for VCF input (e.g., "DEL", "DUP").
#' @param bin_size Bin size for normalized coverage.
#' @param title Custom plot title.
#' @param use_mane Logical, whether to restrict annotations to MANE transcripts.
#' @param include_clinical Logical, if \code{TRUE} also include MANE Plus Clinical transcripts.
#' @param control_colors Optional vector of colors for control coverage points.
#'   If NULL, controls are plotted in gray.
#'
#' @export
plot_gene <- function(
    gene                = NULL,
    hgnc_id             = NULL,
    bam_file,
    control_bams        = NULL,
    vcf_list            = NULL,
    bed_list            = NULL,
    genome              = "hg38",
    out                 = "gene_plot.png",
    chrom               = NULL,
    flank               = 5000,
    show_gene           = TRUE,
    show_cov            = TRUE,
    show_sv             = TRUE,
    sv_types            = "DEL",
    bin_size            = 2000,
    title               = NULL,
    use_mane            = TRUE,
    include_clinical    = TRUE,
    control_colors      = NULL
) {
  
  parts <- document_gene_parts(
    gene             = gene,
    hgnc_id          = hgnc_id,
    genome           = genome,
    use_mane         = use_mane,
    include_clinical = include_clinical
  )
  
  chrom <- if (!is.null(chrom)) {
    sub("^chr", "", chrom)
  } else {
    sub("^chr", "", as.character(GenomicRanges::seqnames(parts)[1]))
  }
  
  gene_label <- if (!is.null(gene)) gene else {
    sym <- get_symbol_from_hgnc(hgnc_id)
    if (length(sym) > 0 && !is.na(sym[1])) sym[1] else "Unknown_gene"
  }
  if (is.null(title)) title <- paste("Focus in", gene_label)
  
  roi <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges   = IRanges::IRanges(
      start = min(GenomicRanges::start(parts)) - flank,
      end   = max(GenomicRanges::end(parts))   + flank
    )
  )
  
  tracks <- list()
  track_sizes <- numeric()
  
  axis_track <- make_axis_track(genome)
  tracks <- c(tracks, list(axis_track))
  track_sizes <- c(track_sizes, 1.5)
  
  if (show_gene) {
    fill_colors <- assign_feature_colors(parts)
    track_gene <- Gviz::AnnotationTrack(
      range      = parts,
      genome     = genome,
      chromosome = chrom,
      name       = gene_label,
      feature    = parts$feature,
      stacking   = "dense",
      col        = "black",
      fill       = fill_colors
    )
    tracks <- c(tracks, list(track_gene))
    track_sizes <- c(track_sizes, 1.5)
  }
  
  if (show_cov) {
    cov_result <- make_patient_control_cov_tracks(
      patient_bam    = bam_file,
      control_bams   = control_bams,
      roi            = roi,
      bin_size       = bin_size,
      chrom          = chrom,
      genome         = genome,
      control_colors = control_colors
    )
    tracks <- c(tracks, cov_result$tracks)
    track_sizes <- c(track_sizes, cov_result$sizes)
  }
  
  if (show_sv) {
    sv_tracks <- list()
    
    if (!is.null(bed_list)) {
      sv_tracks <- make_sv_tracks_from_bed(
        bed_list = bed_list,
        genome   = genome,
        chrom    = chrom
      )
    } else if (!is.null(vcf_list)) {
      sv_tracks <- make_sv_tracks(
        vcf_list = vcf_list,
        sv_types = sv_types,
        genome   = genome,
        chrom    = chrom
      )
    }
    
    if (length(sv_tracks) > 0) {
      tracks <- c(tracks, sv_tracks)
      track_sizes <- c(track_sizes, rep(2.2, length(sv_tracks)))
    }
  }
  
  grDevices::png(out, width = 1600, height = 900, res = 150)
  Gviz::plotTracks(
    tracks,
    from                 = GenomicRanges::start(roi),
    to                   = GenomicRanges::end(roi),
    main                 = title,
    title.width          = 1.8,
    background.title     = "gray95",
    cex.title            = 0.9,
    fontcolor.title      = "black",
    col.axis             = "black",
    col.title            = "black",
    cex.axis             = 0.8,
    showTitle            = TRUE,
    margin               = 80,
    frame                = TRUE,
    sizes                = track_sizes,
    ucscChromosomeNames  = FALSE
  )
  grDevices::dev.off()
}

#' Plot CNV and Structural Variations for a Genomic Region
#'
#' This function visualizes coverage and structural variations for a specific genomic region
#' using BAM and VCF or BED inputs, without requiring gene annotation.
#'
#' @param chrom Chromosome (e.g., "1", "chr1").
#' @param start Start position in base pairs.
#' @param end End position in base pairs.
#' @param bam_file BAM file with alignments.
#' @param control_bams A named list of control BAM file paths for normalization.
#' @param vcf_list A named list of VCF file paths for structural variants.
#' @param bed_list A named list of BED file paths for structural variants (alternative to \code{vcf_list}).
#' @param genome Genome assembly ("hg19" or "hg38").
#' @param out Output image filename (PNG).
#' @param show_cov Logical, whether to show read coverage.
#' @param show_sv Logical, whether to show structural variations.
#' @param sv_types Type(s) of SV to display for VCF input.
#' @param bin_size Bin size for normalized coverage.
#' @param title Custom plot title.
#' @param control_colors Optional vector of colors for control coverage points.
#'   If NULL, controls are plotted in gray.
#'
#' @export
plot_region <- function(
    chrom,
    start,
    end,
    bam_file,
    control_bams        = NULL,
    vcf_list            = NULL,
    bed_list            = NULL,
    genome              = "hg38",
    out                 = "region_plot.png",
    show_cov            = TRUE,
    show_sv             = TRUE,
    sv_types            = "DEL",
    bin_size            = 2000,
    title               = NULL,
    control_colors      = NULL
) {
  
  chrom_clean <- sub("^chr", "", chrom)
  
  if (is.null(title)) {
    title <- paste0(
      "Region ", chrom, ":",
      format(start, big.mark = ","),
      "-",
      format(end,   big.mark = ",")
    )
  }
  
  roi <- GenomicRanges::GRanges(
    seqnames = chrom_clean,
    ranges   = IRanges::IRanges(start, end)
  )
  
  tracks <- list()
  track_sizes <- numeric()
  
  axis_track <- make_axis_track(genome)
  tracks <- c(tracks, list(axis_track))
  track_sizes <- c(track_sizes, 1.5)
  
  if (show_cov) {
    cov_result <- make_patient_control_cov_tracks(
      patient_bam    = bam_file,
      control_bams   = control_bams,
      roi            = roi,
      bin_size       = bin_size,
      chrom          = chrom_clean,
      genome         = genome,
      control_colors = control_colors
    )
    tracks <- c(tracks, cov_result$tracks)
    track_sizes <- c(track_sizes, cov_result$sizes)
  }
  
  if (show_sv) {
    sv_tracks <- list()
    
    if (!is.null(bed_list)) {
      sv_tracks <- make_sv_tracks_from_bed(
        bed_list = bed_list,
        genome   = genome,
        chrom    = chrom
      )
    } else if (!is.null(vcf_list)) {
      sv_tracks <- make_sv_tracks(
        vcf_list = vcf_list,
        sv_types = sv_types,
        genome   = genome,
        chrom    = chrom_clean
      )
    }
    
    if (length(sv_tracks) > 0) {
      tracks <- c(tracks, sv_tracks)
      track_sizes <- c(track_sizes, rep(2.2, length(sv_tracks)))
    }
  }
  
  grDevices::png(out, width = 1600, height = 900, res = 150)
  Gviz::plotTracks(
    tracks,
    from                 = start,
    to                   = end,
    main                 = title,
    title.width          = 1.8,
    background.title     = "gray95",
    cex.title            = 0.9,
    fontcolor.title      = "black",
    col.axis             = "black",
    col.title            = "black",
    cex.axis             = 0.8,
    showTitle            = TRUE,
    margin               = 80,
    frame                = TRUE,
    sizes                = track_sizes,
    ucscChromosomeNames  = FALSE
  )
  grDevices::dev.off()
}

#' Get genes overlapped by SVs
#'
#' Returns gene names that are overlapped by structural variants in a VCF file.
#'
#' @param vcf_file Path to VCF file.
#' @param sv_type Type of SV to consider (e.g., "DEL").
#' @param genome Genome assembly ("hg19" or "hg38").
#'
#' @export
#' @examples
#' # files <- GeGraphs:::load_example_files_HG02317()
#' # get_genes_SV(files$delly, sv_type = "DEL", genome = "hg38")

get_genes_SV <- function(vcf_file, sv_type = "DEL", genome = "hg38") {
  vcf  <- vcfR::read.vcfR(vcf_file)
  info <- vcfR::extract_info_tidy(vcf)
  
  keep <- which(info$SVTYPE == sv_type & !is.na(info$END))
  if (length(keep) == 0) return(character())
  
  gr <- GenomicRanges::GRanges(
    seqnames = paste0("chr", vcfR::getCHROM(vcf)[keep]),
    ranges   = IRanges::IRanges(
      start = vcfR::getPOS(vcf)[keep],
      end   = as.numeric(info$END[keep])
    )
  )
  
  if (genome == "hg38") {
    gff_file <- system.file(
      "extdata",
      "gencode.v47.basic.annotation.gff3.gz",
      package = "GeGraphs"
    )
  } else if (genome == "hg19") {
    gff_file <- system.file(
      "extdata",
      "gencode.v49lift37.basic.annotation.gff3.gz",
      package = "GeGraphs"
    )
  } else {
    stop("Genome must be 'hg38' or 'hg19'")
  }
  
  gff <- rtracklayer::import(gff_file)
  
  gff_genes <- gff[gff$type == "gene" & !is.na(gff$gene_name)]
  
  hits <- GenomicRanges::findOverlaps(gr, gff_genes)
  if (length(hits) == 0) return(character())
  overlapping_genes <- gff_genes[S4Vectors::subjectHits(hits)]
  
  unique(overlapping_genes$gene_name)
}


#' Get genes under SVs common across VCFs
#'
#' Returns genes impacted by SVs shared among multiple VCF files.
#'
#' @param vcf_files A character vector of VCF file paths.
#' @param sv_type Type of SV (e.g., "DEL").
#' @param genome Genome assembly ("hg19" or "hg38").
#'
#' @export
#' @examples
#' # files <- GeGraphs:::load_example_files_HG02317()
#' # get_common_genes_SV(
#' #   vcf_files = c(files$delly, files$manta, files$del),
#' #   sv_type   = "DEL",
#' #   genome    = "hg38"
#' # )

get_common_genes_SV <- function(vcf_files, sv_type = "DEL", genome = "hg38") {
  
  extract_gene_ranges <- function(vcf_path, sv_type, genome) {
    vcf  <- vcfR::read.vcfR(vcf_path)
    info <- vcfR::extract_info_tidy(vcf)
    
    keep <- which(info$SVTYPE == sv_type & !is.na(info$END))
    if (length(keep) == 0) {
      return(GenomicRanges::GRanges())
    }
    
    gr <- GenomicRanges::GRanges(
      seqnames = paste0("chr", vcfR::getCHROM(vcf)[keep]),
      ranges   = IRanges::IRanges(
        start = vcfR::getPOS(vcf)[keep],
        end   = as.numeric(info$END[keep])
      )
    )
    
    if (genome == "hg38") {
      gff_file <- system.file(
        "extdata",
        "gencode.v47.basic.annotation.gff3.gz",
        package = "GeGraphs"
      )
    } else if (genome == "hg19") {
      gff_file <- system.file(
        "extdata",
        "gencode.v49lift37.basic.annotation.gff3.gz",
        package = "GeGraphs"
      )
    } else {
      stop("Genome must be 'hg38' or 'hg19'")
    }
    
    gff <- rtracklayer::import(gff_file)
    
    genes <- gff[gff$type == "gene" & !is.na(gff$gene_name)]
    
    hits <- GenomicRanges::findOverlaps(gr, genes)
    if (length(hits) == 0) {
      return(GenomicRanges::GRanges())
    }
    genes[S4Vectors::subjectHits(hits)]
  }
  
  gene_hits_list  <- lapply(vcf_files, extract_gene_ranges, sv_type = sv_type, genome = genome)
  gene_names_list <- lapply(gene_hits_list, function(gr) gr$gene_name)
  
  if (length(gene_names_list) == 0) return(data.frame())
  
  common_gene_names <- Reduce(intersect, gene_names_list)
  if (length(common_gene_names) == 0) return(data.frame())
  
  common_gene_ranges <- gene_hits_list[[1]][
    gene_hits_list[[1]]$gene_name %in% common_gene_names
  ]
  
  if (length(common_gene_ranges) == 0) return(data.frame())
  
  result <- data.frame(
    gene  = common_gene_ranges$gene_name,
    chr   = as.character(GenomicRanges::seqnames(common_gene_ranges)),
    start = GenomicRanges::start(common_gene_ranges),
    end   = GenomicRanges::end(common_gene_ranges),
    stringsAsFactors = FALSE
  )
  
  unique(result)
}

#' Get genes upstream and downstream of SVs
#'
#' This function returns the closest upstream and downstream genes for structural variants (SVs)
#' provided either as a BED file or as a VCF file. Optionally, SVs can be restricted to those
#' overlapping a specific gene or genomic position.
#'
#' Exactly one of \code{bed_file} or \code{vcf_file} must be provided.
#' At most one of \code{gene} or \code{chrom}/\code{pos} can be specified.
#'
#' @param bed_file Path to BED file with SV coordinates (chr, start, end).
#' @param vcf_file Path to VCF file with SVs.
#' @param sv_type SV type to consider when \code{vcf_file} is used (e.g., "DEL").
#' @param genome Genome assembly ("hg19" or "hg38").
#' @param n_upstream Number of upstream genes to return per SV.
#' @param n_downstream Number of downstream genes to return per SV.
#' @param gene Optional gene symbol (HGNC), e.g. "ERGIC1". If provided, only SVs overlapping
#'   this gene are considered.
#' @param chrom Optional chromosome for position-based filtering (e.g., "5" or "chr5").
#' @param pos Optional 1-based genomic position for position-based filtering.
#'
#' @return
#' A data.frame with one row per SV–gene pair, including at least the columns:
#' \itemize{
#'   \item \code{gene}: gene symbol.
#'   \item \code{direction}: "upstream" or "downstream" relative to the SV midpoint.
#'   \item \code{hgnc_id}: HGNC identifier.
#'   \item \code{mane_type}: MANE type (e.g., "MANE Select", "MANE Plus Clinical").
#'   \item \code{distance}: signed distance (bp) from SV midpoint to gene midpoint.
#' }
#'
#' @export
get_gene_context <- function(
    bed_file     = NULL,
    vcf_file     = NULL,
    sv_type      = "DEL",
    genome       = "hg38",
    n_upstream   = 3,
    n_downstream = 3,
    gene         = NULL,
    chrom        = NULL,
    pos          = NULL
) {
  if (is.null(bed_file) && is.null(vcf_file)) {
    stop("You must provide either 'bed_file' or 'vcf_file'.")
  }
  if (!is.null(bed_file) && !is.null(vcf_file)) {
    stop("Provide only one of 'bed_file' or 'vcf_file', not both.")
  }
  if (!is.null(gene) && (!is.null(chrom) || !is.null(pos))) {
    stop("Specify either 'gene' or 'chrom'+'pos', not both.")
  }
  if (!is.null(chrom) && is.null(pos)) {
    stop("If 'chrom' is provided, 'pos' must also be provided.")
  }
  if (is.null(chrom) && !is.null(pos)) {
    stop("If 'pos' is provided, 'chrom' must also be provided.")
  }
  
  if (!is.null(bed_file)) {
    gr_sv <- bed_to_gr(bed_file)
  } else {
    gr_sv <- vcf_to_gr(vcf_file, type = sv_type)
  }
  if (length(gr_sv) == 0L) {
    return(data.frame())
  }
  
  if (!is.null(gene)) {
    gr_sv <- .filter_sv_by_gene(
      gr_sv  = gr_sv,
      gene   = gene,
      genome = genome
    )
    if (length(gr_sv) == 0L) {
      stop("No SV in input overlaps gene ", gene, " on genome ", genome)
    }
  }
  
  if (!is.null(chrom) && !is.null(pos)) {
    gr_sv <- .filter_sv_by_position(
      gr_sv = gr_sv,
      chrom = chrom,
      pos   = pos
    )
    if (length(gr_sv) == 0L) {
      stop("No SV in input overlaps position ", chrom, ":", pos, " on genome ", genome)
    }
  }
  
  get_gene_context_gr(
    gr_sv        = gr_sv,
    genome       = genome,
    n_upstream   = n_upstream,
    n_downstream = n_downstream
  )
}

#' Get candidate de novo SVs from trio
#'
#' This function identifies candidate de novo structural variants (SVs) in a trio
#' (proband, father, mother) and writes them to a BED file.
#'
#' Exactly one of the BED-based or VCF-based argument sets must be provided:
#' \itemize{
#'   \item BED mode: \code{bed_proband}, \code{bed_father}, \code{bed_mother}
#'   \item VCF mode: \code{vcf_proband}, \code{vcf_father}, \code{vcf_mother}, \code{sv_type}
#' }
#'
#' @param bed_proband BED file with proband SVs (chr, start, end).
#' @param bed_father BED file with father SVs.
#' @param bed_mother BED file with mother SVs.
#' @param vcf_proband VCF file with proband SVs.
#' @param vcf_father VCF file with father SVs.
#' @param vcf_mother VCF file with mother SVs.
#' @param sv_type SV type to consider for VCF input (e.g., "DEL").
#' @param out_bed_denovo Output BED file for candidate de novo SVs.
#' @param out_bed_proband Optional BED output for all proband SVs.
#' @param out_bed_parents Optional BED output for union of parental SVs.
#'
#' @return
#' A \code{GRanges} object containing candidate de novo SV intervals.
#'
#' @export
get_denovo_sv <- function(
    bed_proband     = NULL,
    bed_father      = NULL,
    bed_mother      = NULL,
    vcf_proband     = NULL,
    vcf_father      = NULL,
    vcf_mother      = NULL,
    sv_type         = "DEL",
    out_bed_denovo  = "proband_denovo.bed",
    out_bed_proband = NULL,
    out_bed_parents = NULL
) {
  bed_mode <- !is.null(bed_proband) || !is.null(bed_father) || !is.null(bed_mother)
  vcf_mode <- !is.null(vcf_proband) || !is.null(vcf_father) || !is.null(vcf_mother)
  
  if (!bed_mode && !vcf_mode) {
    stop("Provide either BED inputs (bed_proband/bed_father/bed_mother) ",
         "or VCF inputs (vcf_proband/vcf_father/vcf_mother).")
  }
  if (bed_mode && vcf_mode) {
    stop("Use only one mode: either BED or VCF, not both.")
  }
  
  if (bed_mode) {
    if (is.null(bed_proband) || is.null(bed_father) || is.null(bed_mother)) {
      stop("In BED mode, provide 'bed_proband', 'bed_father' and 'bed_mother'.")
    }
    
    gr_proband <- bed_to_gr(bed_proband)
    gr_father  <- bed_to_gr(bed_father)
    gr_mother  <- bed_to_gr(bed_mother)
  } else {
    if (is.null(vcf_proband) || is.null(vcf_father) || is.null(vcf_mother)) {
      stop("In VCF mode, provide 'vcf_proband', 'vcf_father' and 'vcf_mother'.")
    }
    
    gr_proband <- vcf_to_gr(vcf_proband, type = sv_type)
    gr_father  <- vcf_to_gr(vcf_father,  type = sv_type)
    gr_mother  <- vcf_to_gr(vcf_mother,  type = sv_type)
  }
  
  gr_parents <- c(gr_father, gr_mother)
  
  if (!is.null(out_bed_proband)) {
    gr_to_bed(gr_proband, out_bed_proband)
  }
  if (!is.null(out_bed_parents)) {
    gr_to_bed(gr_parents, out_bed_parents)
  }
  
  denovo_gr <- get_denovo_sv_gr(
    gr_proband = gr_proband,
    gr_father  = gr_father,
    gr_mother  = gr_mother
  )
  
  gr_to_bed(denovo_gr, out_bed_denovo)
  denovo_gr
}

#' Get genes overlapped by structural variants
#'
#' Returns genes whose gene bodies are overlapped by the provided structural
#' variants (SVs), given either as a BED file or a VCF file.
#'
#' Exactly one of \code{bed_file} or \code{vcf_file} must be provided.
#'
#' @param bed_file Path to BED file with SV coordinates (chr, start, end).
#' @param vcf_file Path to VCF file with SVs.
#' @param sv_type  SV type to consider when \code{vcf_file} is used (e.g., "DEL").
#' @param genome   Genome assembly ("hg19" or "hg38").
#'
#' @return
#' A data.frame with one row per gene, containing at least:
#' \itemize{
#'   \item \code{gene}: gene symbol.
#'   \item \code{chr}: chromosome (without "chr").
#'   \item \code{start}: gene start (bp).
#'   \item \code{end}: gene end (bp).
#' }
#'
#' @export
get_genes_under_sv <- function(
    bed_file = NULL,
    vcf_file = NULL,
    sv_type  = "DEL",
    genome   = "hg38"
) {
  if (is.null(bed_file) && is.null(vcf_file)) {
    stop("You must provide either 'bed_file' or 'vcf_file'.")
  }
  if (!is.null(bed_file) && !is.null(vcf_file)) {
    stop("Provide only one of 'bed_file' or 'vcf_file', not both.")
  }
  
  if (!is.null(bed_file)) {
    gr_sv <- bed_to_gr(bed_file)
    valid_chroms <- c(as.character(1:22), "X", "Y", "M", "MT")
    gr_sv <- gr_sv[as.character(GenomicRanges::seqnames(gr_sv)) %in% valid_chroms]
    gr_sv <- GenomeInfoDb::keepSeqlevels(
      gr_sv,
      intersect(GenomeInfoDb::seqlevels(gr_sv), valid_chroms),
      pruning.mode = "coarse"
    )
  } else {
    gr_sv <- vcf_to_gr(vcf_file, type = sv_type)
  }
  
  if (length(gr_sv) == 0L) {
    return(data.frame())
  }
  
  genes <- .get_clean_gene_annotations(genome = genome)
  
  hits <- GenomicRanges::findOverlaps(gr_sv, genes)
  if (length(hits) == 0L) {
    return(data.frame())
  }
  
  overlapping_genes <- genes[S4Vectors::subjectHits(hits)]
  
  res <- data.frame(
    gene  = overlapping_genes$gene,
    chr   = sub("^chr", "", as.character(GenomicRanges::seqnames(overlapping_genes))),
    start = GenomicRanges::start(overlapping_genes),
    end   = GenomicRanges::end(overlapping_genes),
    stringsAsFactors = FALSE
  )
  
  unique(res)
}

#' Find panel genes overlapped by structural variants from a data.frame
#'
#' This function takes a panel already loaded as a data.frame and a set of
#' structural variants (SVs) provided as a BED or VCF file, finds genes under
#' the SVs that are present in the selected panel column, and optionally
#' generates coverage plots for each overlapping gene using \code{plot_gene()}.
#'
#' Exactly one of \code{bedfile} or \code{vcffile} must be provided.
#'
#' @param paneldf A data.frame containing at least one column of gene symbols.
#' @param genecol Name of the column in \code{paneldf} containing gene symbols.
#' @param bedfile Path to BED file with SV coordinates (chr, start, end).
#' @param vcffile Path to VCF file with SVs.
#' @param svtype SV type to consider when \code{vcffile} is used
#'   (e.g., \code{"DEL"}).
#' @param genome Genome assembly (\code{"hg19"} or \code{"hg38"}).
#' @param makeplots Logical; if TRUE, plots are generated for each
#'   overlapping gene.
#' @param bamfile BAM file with alignments for the proband (required if
#'   \code{makeplots = TRUE}).
#' @param controlbams Optional named list of control BAM files.
#' @param plotdir Directory where PNG plots will be saved.
#' @param controlcolors Optional vector of colors for control coverage points.
#'   If NULL, controls are plotted in gray.
#' @param flank Number of base pairs to include upstream/downstream of each
#'   gene in plots.
#' @param binsize Bin size for normalized coverage in plots.
#'
#' @return
#' A data.frame with overlapping panel genes and SV-overlap coordinates, with
#' columns:
#' \itemize{
#'   \item \code{gene}: gene symbol.
#'   \item \code{chr}: chromosome.
#'   \item \code{start}: gene start (bp).
#'   \item \code{end}: gene end (bp).
#' }
#'
#' If \code{makeplots = TRUE}, PNG files are created in \code{plotdir}.
#'
#' @export
panels_overlap_df <- function(
    paneldf,
    genecol,
    bedfile = NULL,
    vcffile = NULL,
    svtype = "DEL",
    genome = "hg38",
    makeplots = FALSE,
    bamfile = NULL,
    controlbams = NULL,
    plotdir = "panelsvplots",
    controlcolors = NULL,
    flank = 5000,
    binsize = 2000
) {
  if (is.null(bedfile) && is.null(vcffile)) {
    stop("You must provide either bedfile or vcffile.")
  }
  if (!is.null(bedfile) && !is.null(vcffile)) {
    stop("Provide only one of bedfile or vcffile, not both.")
  }
  if (!is.data.frame(paneldf)) {
    stop("paneldf must be a data.frame.")
  }
  
  colnames(paneldf) <- trimws(colnames(paneldf))
  
  if (!genecol %in% colnames(paneldf)) {
    stop("Column ", genecol, " not found in paneldf.")
  }
  
  panelgenes <- unique(trimws(as.character(paneldf[[genecol]])))
  panelgenes <- toupper(panelgenes)
  panelgenes <- panelgenes[!is.na(panelgenes) & panelgenes != ""]
  
  if (length(panelgenes) == 0L) {
    stop("No non-empty gene symbols found in column ", genecol, ".")
  }
  
  genesunder <- get_genes_under_sv(
    bed_file = bedfile,
    vcf_file = vcffile,
    sv_type  = svtype,
    genome   = genome
  )
  
  if (nrow(genesunder) == 0L) {
    warning("No genes under SVs; returning empty data.frame.")
    return(data.frame())
  }
  
  genesunder$gene    <- as.character(genesunder$gene)
  genesunder$genenorm <- toupper(trimws(genesunder$gene))
  
  overlapidx <- which(genesunder$genenorm %in% panelgenes)
  
  if (length(overlapidx) == 0L) {
    warning("No overlapping genes between panel and genes under SVs.")
    return(data.frame())
  }
  
  overlapgenes <- unique(
    genesunder[overlapidx, c("gene", "chr", "start", "end"), drop = FALSE]
  )
  
  if (makeplots) {
    if (is.null(bamfile)) {
      stop("bamfile must be provided when makeplots = TRUE.")
    }
    if (!file.exists(bamfile)) {
      stop("bamfile not found: ", bamfile)
    }
    if (!dir.exists(plotdir)) {
      dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)
    }
    
    for (g in overlapgenes$gene) {
      outpng <- file.path(plotdir, paste0(g, "_panelsv.png"))
      
      tryCatch(
        plot_gene(
          gene           = g,
          bam_file       = bamfile,
          control_bams   = controlbams,
          genome         = genome,
          out            = outpng,
          flank          = flank,
          show_gene      = TRUE,
          show_cov       = TRUE,
          show_sv        = FALSE,
          bin_size       = binsize,
          control_colors = controlcolors
        ),
        error = function(e) {
          warning("Failed to plot gene ", g, ": ", conditionMessage(e))
        }
      )
    }
  }
  
  overlapgenes
}

#' Find panel genes overlapped by structural variants and optionally plot them
#'
#' This function reads a gene panel from a TSV, TXT or CSV file, then calls
#' \code{panels_overlap_df()} to find genes under the SVs that are present in
#' the selected panel column, and optionally generate coverage plots.
#'
#' Exactly one of \code{bedfile} or \code{vcffile} must be provided.
#'
#' @param panelfile Path to panel file (\code{.tsv}, \code{.txt} or
#'   \code{.csv}) containing at least one column of gene symbols.
#' @param panelgenecol Name of the column in \code{panelfile} containing gene
#'   symbols.
#' @param bedfile Path to BED file with SV coordinates (chr, start, end).
#' @param vcffile Path to VCF file with SVs.
#' @param svtype SV type to consider when \code{vcffile} is used
#'   (e.g., \code{"DEL"}).
#' @param genome Genome assembly (\code{"hg19"} or \code{"hg38"}).
#' @param makeplots Logical; if TRUE, plots are generated for each
#'   overlapping gene.
#' @param bamfile BAM file with alignments for the proband (required if
#'   \code{makeplots = TRUE}).
#' @param controlbams Optional named list of control BAM files.
#' @param plotdir Directory where PNG plots will be saved.
#' @param controlcolors Optional vector of colors for control coverage points.
#'   If NULL, controls are plotted in gray.
#' @param flank Number of base pairs to include upstream/downstream of each
#'   gene in plots.
#' @param binsize Bin size for normalized coverage in plots.
#'
#' @return
#' A data.frame with overlapping panel genes and SV-overlap coordinates, with
#' columns:
#' \itemize{
#'   \item \code{gene}: gene symbol.
#'   \item \code{chr}: chromosome.
#'   \item \code{start}: gene start (bp).
#'   \item \code{end}: gene end (bp).
#' }
#'
#' If \code{makeplots = TRUE}, PNG files are created in \code{plotdir}.
#'
#' @examples
#' files <- GeGraphs::loadexamplefiles("HG02317")
#' tf <- tempfile(fileext = ".tsv")
#' write.table(
#'   data.frame(gene = c("NELL1", "PCDH19"), stringsAsFactors = FALSE),
#'   file = tf,
#'   sep = "\t",
#'   quote = FALSE,
#'   row.names = FALSE
#' )
#' panels_overlap(
#'   panelfile = tf,
#'   panelgenecol = "gene",
#'   vcffile = files$delly,
#'   svtype = "DEL",
#'   genome = "hg38"
#' )
#'
#' @export
panels_overlap <- function(
    panelfile,
    panelgenecol = "gene",
    bedfile = NULL,
    vcffile = NULL,
    svtype = "DEL",
    genome = "hg38",
    makeplots = FALSE,
    bamfile = NULL,
    controlbams = NULL,
    plotdir = "panelsvplots",
    controlcolors = NULL,
    flank = 5000,
    binsize = 2000
) {
  if (!file.exists(panelfile)) {
    stop("panelfile not found: ", panelfile)
  }
  
  ext <- tolower(tools::file_ext(panelfile))
  
  if (ext %in% c("tsv", "txt")) {
    paneldf <- read.delim(
      panelfile,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else if (ext %in% c("csv")) {
    paneldf <- read.csv(
      panelfile,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    stop("panelfile must be .tsv, .txt or .csv")
  }
  
  panels_overlap_df(
    paneldf = paneldf,
    genecol = panelgenecol,
    bedfile = bedfile,
    vcffile = vcffile,
    svtype = svtype,
    genome = genome,
    makeplots = makeplots,
    bamfile = bamfile,
    controlbams = controlbams,
    plotdir = plotdir,
    controlcolors = controlcolors,
    flank = flank,
    binsize = binsize
  )
}

#' Get upstream/downstream gene context for a specific gene
#'
#' Convenience wrapper around get_gene_context() that restricts to SVs
#' overlapping a given gene symbol.
#'
#' @param gene Gene symbol (HGNC), e.g. "ERGIC1".
#' @param bed_file BED file with SVs (chr, start, end).
#' @param vcf_file VCF file with SVs (alternative to \code{bed_file}).
#' @param sv_type SV type for VCF input (e.g. "DEL").
#' @param genome Genome assembly ("hg19" or "hg38").
#' @param n_upstream Number of upstream genes to return.
#' @param n_downstream Number of downstream genes to return.
#'
#' @export
get_gene_context_for_gene <- function(
    gene,
    bed_file     = NULL,
    vcf_file     = NULL,
    sv_type      = "DEL",
    genome       = "hg38",
    n_upstream   = 3,
    n_downstream = 3
) {
  get_gene_context(
    bed_file     = bed_file,
    vcf_file     = vcf_file,
    sv_type      = sv_type,
    genome       = genome,
    n_upstream   = n_upstream,
    n_downstream = n_downstream,
    gene         = gene
  )
}


#' @keywords internal
assign_feature_colors <- function(parts) {
  feature_colors <- c(
    "five_prime_UTR" = "darkgreen",
    "CDS"            = "purple",
    "three_prime_UTR"= "orange"
  )
  features    <- as.character(parts$feature)
  fill_colors <- feature_colors[features]
  fill_colors[is.na(fill_colors)] <- "gray"
  fill_colors
}


#' @keywords internal
bed_to_gr <- function(bed_file) {
  bed <- read.table(
    bed_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (ncol(bed) < 3L) {
    stop("BED file must have at least 3 columns: chr, start, end")
  }
  
  GenomicRanges::GRanges(
    seqnames = sub("^chr", "", bed[[1]]),
    ranges   = IRanges::IRanges(
      start = as.integer(bed[[2]]) + 1L,
      end   = as.integer(bed[[3]])
    )
  )
}

#' @keywords internal
make_sv_tracks_from_bed <- function(bed_list, genome, chrom) {
  tracks <- list()
  
  if (is.null(names(bed_list))) {
    names(bed_list) <- paste0("BED", seq_along(bed_list))
  }
  
  chrom_clean <- sub("^chr", "", chrom)
  valid_chroms <- c(as.character(1:22), "X", "Y", "M", "MT")
  
  for (bed_name in names(bed_list)) {
    bed_path <- bed_list[[bed_name]]
    gr <- bed_to_gr(bed_path)
    gr <- gr[as.character(GenomicRanges::seqnames(gr)) %in% valid_chroms]
    gr <- GenomeInfoDb::keepSeqlevels(
      gr,
      intersect(GenomeInfoDb::seqlevels(gr), valid_chroms),
      pruning.mode = "coarse"
    )
    if (length(gr) == 0L) {
      next
    }
    
    gr_chr <- gr[as.character(GenomicRanges::seqnames(gr)) == chrom_clean]
    gr_chr <- GenomeInfoDb::keepSeqlevels(
      gr_chr,
      chrom_clean,
      pruning.mode = "coarse"
    )
    if (length(gr_chr) == 0L) {
      next
    }
    
    track <- Gviz::AnnotationTrack(
      range = gr_chr,
      genome = genome,
      chromosome = chrom_clean,
      name = bed_name,
      fill = "orange",
      col = "black"
    )
    
    tracks <- c(tracks, list(track))
  }
  
  tracks
}

#' @keywords internal
vcf_to_gr <- function(file, type) {
  vcf  <- vcfR::read.vcfR(file)
  info <- vcfR::extract_info_tidy(vcf)
  keep <- which(info$SVTYPE == type & !is.na(info$END))
  if (length(keep) == 0) {
    return(GenomicRanges::GRanges())
  }
  
  chroms <- sub("^chr", "", vcfR::getCHROM(vcf)[keep])
  valid_chroms <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                    "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                    "21", "22", "X", "Y", "M", "MT")
  valid_idx <- which(chroms %in% valid_chroms)
  
  if (length(valid_idx) == 0) {
    return(GenomicRanges::GRanges())
  }
  
  GenomicRanges::GRanges(
    seqnames = chroms[valid_idx],
    ranges   = IRanges::IRanges(
      start = vcfR::getPOS(vcf)[keep][valid_idx],
      end   = as.numeric(info$END[keep][valid_idx])
    )
  )
}


#' @keywords internal
make_sv_tracks <- function(vcf_list, sv_types = "DEL", genome, chrom) {
  tracks <- list()
  
  if (is.null(names(vcf_list))) {
    names(vcf_list) <- paste0("VCF", seq_along(vcf_list))
  }
  
  for (vcf_name in names(vcf_list)) {
    vcf_path <- vcf_list[[vcf_name]]
    
    for (sv_type in sv_types) {
      gr <- vcf_to_gr(vcf_path, type = sv_type)
      
      if (length(gr) > 0) {
        track_name <- paste(vcf_name, sv_type, sep = " ")
        
        track <- Gviz::AnnotationTrack(
          range      = gr,
          genome     = genome,
          chromosome = chrom,
          name       = track_name,
          fill       = "orange",
          col        = "black"
        )
        tracks <- c(tracks, list(track))
      }
    }
  }
  
  tracks
}

#' @keywords internal
make_axis_track <- function(genome) {
  Gviz::GenomeAxisTrack(
    genome           = genome,
    name             = "Genome Position",
    littleTicks      = TRUE,
    col              = "black",
    background.title = "white",
    fontcolor.title  = "black"
  )
}

#' @keywords internal
document_gene_parts <- function(gene            = NULL,
                                hgnc_id         = NULL,
                                genome          = "hg38",
                                use_mane        = TRUE,
                                include_clinical = TRUE) {
  if (!is.null(hgnc_id)) {
    gene_symbol <- get_symbol_from_hgnc(hgnc_id)
    if (length(gene_symbol) == 0 || is.na(gene_symbol[1])) {
      stop("No symbol found for HGNC_ID = ", hgnc_id)
    }
    gene_symbol <- gene_symbol[1]
  } else if (!is.null(gene)) {
    gene_symbol <- gene
  } else {
    stop("You must provide either 'gene' or 'hgnc_id'.")
  }
  
  if (genome == "hg38") {
    gff_file <- system.file("extdata", "gencode.v47.basic.annotation.gff3.gz",
                            package = "GeGraphs")
  } else if (genome == "hg19") {
    gff_file <- system.file("extdata", "gencode.v49lift37.basic.annotation.gff3.gz",
                            package = "GeGraphs")
  } else {
    stop("Genome must be 'hg38' or 'hg19'")
  }
  
  gff <- rtracklayer::import(gff_file)
  
  parts <- gff[gff$gene_name == gene_symbol &
                 gff$type %in% c("five_prime_UTR", "CDS", "three_prime_UTR")]
  
  if (length(parts) == 0) {
    stop("No annotation found for gene '", gene_symbol,
         "' in GFF ", basename(gff_file))
  }
  
  if (use_mane) {
    mane_tx <- get_mane_transcripts(gene_symbol, include_clinical = include_clinical)
    
    if (length(mane_tx) > 0) {
      parts <- parts[parts$transcript_id %in% mane_tx]
      if (length(parts) == 0) {
        warning("MANE transcripts found for ", gene_symbol,
                " but none match transcript_id in this GFF. ",
                "All GENCODE annotations for this gene are kept.")
        parts <- gff[gff$gene_name == gene_symbol &
                       gff$type %in% c("five_prime_UTR", "CDS", "three_prime_UTR")]
      }
    } else {
      warning("No MANE transcript found for ", gene_symbol,
              " -> all GENCODE transcripts are used.")
    }
  }
  
  seq_chr <- sub("^chr", "", as.character(GenomicRanges::seqnames(parts)))
  valid_chr <- c(as.character(1:22), "X", "Y", "M", "MT")
  
  parts <- parts[seq_chr %in% valid_chr]
  
  if (length(parts) == 0L) {
    stop("No standard-chromosome annotation left for gene '", gene_symbol, "'.")
  }
  
  parts$feature <- parts$type
  parts$group   <- parts$transcript_id
  parts$id      <- parts$exon_id
  parts
}

#' @keywords internal
load_mane_table <- function() {
  mane_file <- system.file("extdata", "mane_select_minimal_with_hgnc.tsv",
                           package = "GeGraphs")
  tab <- read.delim(mane_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(tab) <- c(
    "hgnc_id",
    "gene_name",
    "mane_refseq",
    "mane_ensembl_transcript",
    "mane_type"
  )
  tab
}

#' @keywords internal
get_mane_transcripts <- function(gene_symbol, include_clinical = TRUE) {
  tab <- load_mane_table()
  x <- tab[tab$gene_name == gene_symbol, , drop = FALSE]
  
  if (!include_clinical && "mane_type" %in% colnames(tab)) {
    x <- x[x$mane_type == "MANE Select", , drop = FALSE]
  }
  
  unique(x$mane_ensembl_transcript)
}

#' @keywords internal
get_symbol_from_hgnc <- function(hgnc_id) {
  tab <- load_mane_table()
  x <- tab[tab$hgnc_id == hgnc_id, , drop = FALSE]
  unique(x$gene_name)
}

#' @keywords internal
get_hgnc_id <- function(gene_symbol) {
  tab <- load_mane_table()
  x <- tab[tab$gene_name == gene_symbol, , drop = FALSE]
  unique(x$hgnc_id)
}

#' @keywords internal
get_gene_context_gr <- function(
    gr_sv,
    genome       = "hg38",
    n_upstream   = 3,
    n_downstream = 3
) {
  if (length(gr_sv) == 0L) {
    return(data.frame())
  }
  
  if (genome == "hg38") {
    gff_file <- system.file(
      "extdata",
      "gencode.v47.basic.annotation.gff3.gz",
      package = "GeGraphs"
    )
  } else if (genome == "hg19") {
    gff_file <- system.file(
      "extdata",
      "gencode.v49lift37.basic.annotation.gff3.gz",
      package = "GeGraphs"
    )
  } else {
    stop("Genome must be 'hg38' or 'hg19'")
  }
  
  gff <- rtracklayer::import(gff_file)
  genes <- gff[gff$type == "gene" & !is.na(gff$gene_name)]
  
  mane_tab <- load_mane_table()
  mane_tab <- mane_tab[, c("hgnc_id", "gene_name", "mane_type"), drop = FALSE]
  colnames(mane_tab) <- c("hgnc_id", "gene", "mane_type")
  
  genes_df <- data.frame(
    gene  = genes$gene_name,
    chr   = sub("^chr", "", as.character(GenomicRanges::seqnames(genes))),
    start = GenomicRanges::start(genes),
    end   = GenomicRanges::end(genes),
    stringsAsFactors = FALSE
  )
  
  genes_df <- merge(
    genes_df,
    mane_tab,
    by    = "gene",
    all.x = TRUE,
    sort  = FALSE
  )
  
  genes_df <- genes_df[!is.na(genes_df$hgnc_id), , drop = FALSE]
  
  sv_chr_all <- sub("^chr", "", as.character(GenomicRanges::seqnames(gr_sv)))
  
  res_list <- lapply(seq_along(gr_sv), function(i) {
    sv <- gr_sv[i]
    sv_chr <- sv_chr_all[i]
    
    chr_genes <- genes_df[genes_df$chr == sv_chr, , drop = FALSE]
    if (nrow(chr_genes) == 0L) {
      return(NULL)
    }
    
    sv_mid <- floor((GenomicRanges::start(sv) + GenomicRanges::end(sv)) / 2)
    gene_mid <- floor((chr_genes$start + chr_genes$end) / 2)
    
    dist <- gene_mid - sv_mid
    
    chr_genes$direction <- ifelse(dist < 0, "upstream", "downstream")
    chr_genes$distance  <- dist
    
    upstream <- chr_genes[chr_genes$direction == "upstream", , drop = FALSE]
    upstream <- upstream[order(abs(upstream$distance)), , drop = FALSE]
    if (nrow(upstream) > n_upstream) {
      upstream <- upstream[seq_len(n_upstream), , drop = FALSE]
    }
    
    downstream <- chr_genes[chr_genes$direction == "downstream", , drop = FALSE]
    downstream <- downstream[order(abs(downstream$distance)), , drop = FALSE]
    if (nrow(downstream) > n_downstream) {
      downstream <- downstream[seq_len(n_downstream), , drop = FALSE]
    }
    
    sel <- rbind(upstream, downstream)
    if (nrow(sel) == 0L) {
      return(NULL)
    }
    
    data.frame(
      gene      = sel$gene,
      direction = sel$direction,
      hgnc_id   = sel$hgnc_id,
      mane_type = sel$mane_type,
      stringsAsFactors = FALSE
    )
  })
  
  res_list <- res_list[!vapply(res_list, is.null, logical(1))]
  if (length(res_list) == 0L) {
    return(data.frame())
  }
  
  res <- unique(do.call(rbind, res_list))
  res
}

#' @keywords internal
get_denovo_sv_gr <- function(
    gr_proband,
    gr_father,
    gr_mother
) {
  gr_parents <- c(gr_father, gr_mother)
  
  if (length(gr_proband) == 0L) {
    return(GenomicRanges::GRanges())
  }
  
  if (length(gr_parents) == 0L) {
    return(gr_proband)
  }
  
  hits <- GenomicRanges::findOverlaps(gr_proband, gr_parents)
  keep_idx <- setdiff(seq_len(length(gr_proband)), S4Vectors::queryHits(hits))
  gr_proband[keep_idx]
}

#' @keywords internal
gr_to_bed <- function(gr, out_bed) {
  if (length(gr) == 0L) {
    write.table(
      data.frame(),
      file      = out_bed,
      sep       = "\t",
      quote     = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    return(invisible(NULL))
  }
  
  df <- data.frame(
    chr   = paste0("chr", as.character(GenomicRanges::seqnames(gr))),
    start = GenomicRanges::start(gr) - 1L,
    end   = GenomicRanges::end(gr)
  )
  
  write.table(
    df,
    file      = out_bed,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

#' @keywords internal
.filter_sv_by_gene <- function(gr_sv, gene, genome) {
  if (genome == "hg38") {
    gff_file <- system.file(
      "extdata",
      "gencode.v47.basic.annotation.gff3.gz",
      package = "GeGraphs"
    )
  } else if (genome == "hg19") {
    gff_file <- system.file(
      "extdata",
      "gencode.v49lift37.basic.annotation.gff3.gz",
      package = "GeGraphs"
    )
  } else {
    stop("Genome must be 'hg38' or 'hg19'")
  }
  
  gff <- rtracklayer::import(gff_file)
  genes <- gff[gff$type == "gene" & !is.na(gff$gene_name)]
  gene_gr <- genes[genes$gene_name == gene]
  if (length(gene_gr) == 0L) {
    stop("Gene '", gene, "' not found in GENCODE annotations for ", genome)
  }
  
  sv_chr_all   <- sub("^chr", "", as.character(GenomicRanges::seqnames(gr_sv)))
  gene_chr_all <- sub("^chr", "", as.character(GenomicRanges::seqnames(gene_gr)))
  
  same_chr <- which(sv_chr_all == gene_chr_all[1])
  if (length(same_chr) == 0L) {
    return(GenomicRanges::GRanges())
  }
  
  gr_sv_same_chr <- gr_sv[same_chr]
  sv_df <- data.frame(
    chr   = sv_chr_all[same_chr],
    start = GenomicRanges::start(gr_sv_same_chr),
    end   = GenomicRanges::end(gr_sv_same_chr),
    idx   = seq_along(gr_sv_same_chr),
    stringsAsFactors = FALSE
  )
  
  gene_chr   <- gene_chr_all[1]
  gene_start <- GenomicRanges::start(gene_gr)[1]
  gene_end   <- GenomicRanges::end(gene_gr)[1]
  
  overlaps_idx <- which(
    sv_df$chr == gene_chr &
      sv_df$end   >= gene_start &
      sv_df$start <= gene_end
  )
  
  if (length(overlaps_idx) == 0L) {
    return(GenomicRanges::GRanges())
  }
  
  gr_sv_same_chr[sv_df$idx[overlaps_idx]]
}

#' @keywords internal
.filter_sv_by_position <- function(gr_sv, chrom, pos) {
  chrom_clean <- sub("^chr", "", chrom)
  sv_chr_all  <- sub("^chr", "", as.character(GenomicRanges::seqnames(gr_sv)))
  
  sv_df <- data.frame(
    chr   = sv_chr_all,
    start = GenomicRanges::start(gr_sv),
    end   = GenomicRanges::end(gr_sv),
    idx   = seq_along(gr_sv),
    stringsAsFactors = FALSE
  )
  
  overlaps_idx <- which(
    sv_df$chr == chrom_clean &
      sv_df$end   >= pos &
      sv_df$start <= pos
  )
  
  if (length(overlaps_idx) == 0L) {
    return(GenomicRanges::GRanges())
  }
  
  gr_sv[sv_df$idx[overlaps_idx]]
}

#' @keywords internal
make_patient_control_cov_tracks <- function(
    patient_bam,
    control_bams,
    roi,
    bin_size       = 1500,
    chrom,
    genome,
    control_colors = NULL
) {
  
  tracks <- list()
  track_sizes <- numeric()
  
  region_start <- GenomicRanges::start(roi)
  region_end   <- GenomicRanges::end(roi)
  region_len   <- as.integer(region_end - region_start + 1L)
  
  n_bins <- ceiling(region_len / bin_size)
  bin_starts <- seq(region_start, region_end, by = bin_size)
  bin_pos <- as.integer(bin_starts + bin_size / 2)
  if (length(bin_pos) > n_bins) bin_pos <- bin_pos[1:n_bins]
  
  patient_raw <- get_norm_cov(patient_bam, roi, bin_size = bin_size, bin_pos = bin_pos)
  raw_pos <- patient_raw$raw_pos
  raw_cov <- patient_raw$vec
  
  raw_track <- Gviz::DataTrack(
    range = GenomicRanges::GRanges(
      seqnames = chrom,
      ranges   = IRanges::IRanges(raw_pos, width = 1)
    ),
    data  = raw_cov,
    genome     = genome,
    chromosome = chrom,
    name       = "Raw Coverage",
    type       = "histogram",
    col.histogram  = "gray30",
    fill.histogram = "gray70"
  )
  
  tracks <- c(tracks, list(raw_track))
  track_sizes <- c(track_sizes, 3)
  
  patient_cov <- get_norm_cov(patient_bam, roi, bin_size = bin_size, bin_pos = bin_pos)
  
  len_target <- length(patient_cov$norm)
  gr_patient <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges   = IRanges::IRanges(bin_pos, width = 1)
  )
  S4Vectors::mcols(gr_patient) <- data.frame(cov = patient_cov$norm)
  
  overlay_tracks <- list()
  
  if (!is.null(control_bams) && length(control_bams) > 0) {
    for (i in seq_along(control_bams)) {
      control_cov <- get_norm_cov(control_bams[[i]], roi, bin_size = bin_size, bin_pos = bin_pos)
      
      control_norm <- control_cov$norm
      if (length(control_norm) < len_target) {
        control_norm <- c(control_norm, rep(0, len_target - length(control_norm)))
      } else if (length(control_norm) > len_target) {
        control_norm <- control_norm[seq_len(len_target)]
      }
      
      gr_ctrl <- GenomicRanges::GRanges(
        seqnames = chrom,
        ranges   = IRanges::IRanges(bin_pos, width = 1)
      )
      S4Vectors::mcols(gr_ctrl) <- data.frame(cov = control_norm)
      
      ctrl_name <- names(control_bams)[i]
      if (is.null(ctrl_name) || ctrl_name == "") ctrl_name <- paste0("Control_", i)
      
      ctrl_col <- if (!is.null(control_colors) && length(control_colors) >= i) {
        control_colors[[i]]
      } else {
        rgb(0.5, 0.5, 0.5, 0.4)   ## gris plus clair / plus transparent
      }
      
      ctrl_track <- Gviz::DataTrack(
        range      = gr_ctrl,
        genome     = genome,
        chromosome = chrom,
        name       = ctrl_name,
        type       = "p",
        col        = ctrl_col,
        cex        = 1.5,
        pch        = 16,
        lwd        = 1
      )
      overlay_tracks <- c(overlay_tracks, list(ctrl_track))
    }
  }
  
  patient_track <- Gviz::DataTrack(
    range      = gr_patient,
    genome     = genome,
    chromosome = chrom,
    name       = "Norm. Coverage",
    type       = "p",
    col        = "red3",
    cex        = 2.2,
    pch        = 16,
    lwd        = 2
  )
  overlay_tracks <- c(overlay_tracks, list(patient_track))
  
  overlay_track <- Gviz::OverlayTrack(trackList = overlay_tracks)
  
  tracks <- c(tracks, list(overlay_track))
  track_sizes <- c(track_sizes, 4)
  
  list(tracks = tracks, sizes = track_sizes)
}

#' @keywords internal
get_norm_cov <- function(bam_file, roi, bin_size = 2000, bin_pos = NULL) {
  region_start <- GenomicRanges::start(roi)
  region_end   <- GenomicRanges::end(roi)
  region_len   <- as.integer(region_end - region_start + 1L)
  
  if (is.null(bin_pos)) {
    n_bins <- ceiling(region_len / bin_size)
    bin_starts <- seq(region_start, region_end, by = bin_size)
    bin_pos <- as.integer(bin_starts + bin_size / 2)
    if (length(bin_pos) > n_bins) bin_pos <- bin_pos[1:n_bins]
  }
  
  bam <- Rsamtools::BamFile(bam_file)
  pile <- Rsamtools::pileup(
    file         = bam,
    scanBamParam = Rsamtools::ScanBamParam(which = roi),
    pileupParam  = Rsamtools::PileupParam(include_deletions = TRUE, include_insertions = TRUE)
  )
  
  vec <- numeric(region_len)
  if (NROW(pile) > 0) {
    idx <- as.integer(pile$pos - region_start + 1L)
    valid <- which(idx >= 1L & idx <= region_len)
    if (length(valid) > 0) {
      vec[idx[valid]] <- pile$count[valid]
    }
  }
  
  total_cov <- sum(vec)
  if (!is.finite(total_cov) || total_cov <= 0) total_cov <- 1
  
  n_bins <- length(bin_pos)
  norm_summary <- numeric(n_bins)
  for (b in 1:n_bins) {
    idx_start <- (b - 1) * bin_size + 1
    idx_end   <- min(b * bin_size, region_len)
    norm_summary[b] <- sum(vec[idx_start:idx_end]) / total_cov * 1e6
  }
  
  max_norm <- max(norm_summary, na.rm = TRUE)
  if (max_norm > 0) {
    norm_summary_scaled <- norm_summary / max_norm
  } else {
    norm_summary_scaled <- norm_summary
  }
  
  list(
    norm    = norm_summary_scaled, 
    pos     = bin_pos,
    vec     = vec,
    raw_pos = seq(region_start, region_end)
  )
}

#' Internal: load and clean GENCODE gene annotations
#'
#' Returns a GRanges with one range per gene, including metadata:
#'   - gene      : gene symbol (gencode gene_name)
#'   - hgnc_id   : HGNC identifier (if available, otherwise NA)
#'   - mane_type : MANE type (e.g., "MANE Select", "MANE Plus Clinical", or NA)
#'
#' @param genome "hg19" or "hg38"
#'
#' @keywords internal
.get_clean_gene_annotations <- function(genome = "hg38") {
  if (genome == "hg38") {
    gff_file <- system.file(
      "extdata",
      "gencode.v47.basic.annotation.gff3.gz",
      package = "GeGraphs"
    )
  } else if (genome == "hg19") {
    gff_file <- system.file(
      "extdata",
      "gencode.v49lift37.basic.annotation.gff3.gz",
      package = "GeGraphs"
    )
  } else {
    stop("Genome must be 'hg38' or 'hg19'")
  }
  
  if (gff_file == "" || !file.exists(gff_file)) {
    stop("Could not find GENCODE GFF file for genome = ", genome)
  }
  
  gff <- rtracklayer::import(gff_file)
  
  gff_genes <- gff[gff$type == "gene" & !is.na(gff$gene_name)]
  if (length(gff_genes) == 0L) {
    stop("No gene entries found in GENCODE GFF for genome = ", genome)
  }
  
  old_levels <- GenomeInfoDb::seqlevels(gff_genes)
  new_levels <- sub("^chr", "", old_levels)
  GenomeInfoDb::seqlevels(gff_genes) <- new_levels
  
  keep_chr <- c(as.character(1:22), "X", "Y", "MT", "M")
  gff_genes <- gff_genes[
    as.character(GenomicRanges::seqnames(gff_genes)) %in% keep_chr
  ]
  if (length(gff_genes) == 0L) {
    stop("No standard chromosomes (1-22, X, Y, MT) found in GENCODE GFF for genome = ", genome)
  }
  
  hgnc_id <- if (!is.null(gff_genes$hgnc_id)) gff_genes$hgnc_id else NA_character_
  
  mane_type <- NA_character_
  if (!is.null(gff_genes$tag)) {
    tags_list <- gff_genes$tag
    mane_type <- vapply(
      X = as.list(tags_list),
      FUN = function(x) {
        if (is.null(x)) return(NA_character_)
        if ("MANE_Select" %in% x) return("MANE Select")
        if ("MANE_Plus_Clinical" %in% x) return("MANE Plus Clinical")
        NA_character_
      },
      FUN.VALUE = character(1)
    )
  }
  
  gff_genes$gene      <- as.character(gff_genes$gene_name)
  gff_genes$hgnc_id   <- as.character(hgnc_id)
  gff_genes$mane_type <- as.character(mane_type)
  
  gff_genes
}