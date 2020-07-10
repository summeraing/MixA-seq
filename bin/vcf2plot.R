#!/usr/bin/env Rscript

library(optparse)

option_list = list(
  make_option(c("-s", "--sample"), type = "character",
              help = "Required. Sample name.",
              metavar = "character"),
  make_option(c("-p", "--parent"), type = "character",
              help = "Required. Parent names, separated by ','.",
              metavar = "character"),
  make_option(c("-i", "--index"), type = "character",
              help = "Required. FASTA index file is made from a FASTA file using samtools.",
              metavar = "character"),
  make_option(c("-c", "--chr"),  type = "character",
              help = "Required. Sequence names to plot, separated by ','.",
              metavar = "character"),
  make_option(c("-g", "--genome"), type = "character", default = "hsa",
              help = "Genome abbreviate name [default %default].",
              metavar = "character"),
  make_option(c("-d", "--draw"), type = "integer", default = 1000,
              help = "Draw window size: 1 means 1KB [default %default].",
              metavar = "character"),
  make_option(c("-w", "--window"), type = "integer", default = 1000,
              help = "Genome window size: 1 means 1KB [default %default].",
              metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

#opt <- {}
#opt$sample <- "F2-9"
#opt$parent <- "AL8_78,T093"
#opt$index  <- ref.fasta.fai
#opt$genome <- "aps01"
#opt$window <- 10000
#opt$chr <- "Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7"

vcf2plot <- function(vcf_files, sample_names, genome_name, groups, ref_style = "UCSC", check_alleles = TRUE)
{
    # Check sample names
    if (length(vcf_files) != length(sample_names))
        stop("Please provide the same number of sample names as VCF files")

    library(VariantAnnotation)
    # Detect the number of available cores.
    num_cores = detectCores()
    if (!(.Platform$OS.type == "windows" || is.na(num_cores)))
        num_cores <- detectCores()
    else
        num_cores = 1

    # We handle errors and warnings separately for mclapply, because the error
    # reporting of mclapply is done through its return value(s).
    original_warn_state = getOption("warn")
    options(warn=-1)

    # Store warning messages in a vector.
    warnings <- NULL

    # Show the warning once for all VCF files
    if (!check_alleles)
    {
        warnings <- c(warnings,
                      paste("check_alleles is set to FALSE.  Make sure your",
                            "input VCF does not contain any positions with",
                            "insertions, deletions or multiple alternative",
                            "alleles, as these positions cannot be analysed",
                            "with MutationalPatterns and cause obscure",
                            "errors."))
    }

    vcf_list <- mclapply (seq_along(vcf_files), function (index)
    {
        file <- vcf_files[index]

        # Use VariantAnnotation's readVcf
        vcf <- rowRanges(readVcf (file, genome_name))

        groups <- unique(as.vector(t(groups)))

        # The provided VCF files may not contain all chromosomes that are
        # available in the reference genome.  Therefore, we only take the
        # chromosomes that are actually available in the VCF file,
        # belonging to the filter group.
        groups <- intersect(groups, seqlevels(vcf))

        # We use 'pruning.mode = "tidy"' to minimize the deleterious effect
        # on variants, yet, remove all variants that aren't in the filter
        # group.  By default, keepSeqlevels would produce an error.
        vcf <- keepSeqlevels(vcf, groups, pruning.mode = "tidy")

        if (check_alleles)
        {
            # Find and exclude positions with indels or multiple
            # alternative alleles.
            rem <- which(all(!( !is.na(match(vcf$ALT, DNA_BASES)) &
                                !is.na(match(vcf$REF, DNA_BASES)) &
                                (lengths(vcf$ALT) == 1) )))

            if (length(rem) > 0)
            {
                vcf = vcf[-rem]
                warnings <- c(warnings,
                              paste(length(rem),
                                    "position(s) with indels and/or multiple",
                                    "alternative alleles are excluded from",
                                    paste(sample_names[[index]], ".", sep = "")))
            }
        }

        # Pack GRanges object and the warnings to be able to display warnings
        # at a later time.
        return(list(vcf, warnings))
    }, mc.cores = num_cores)

    # Reset the option.
    options(warn=original_warn_state)

    vcf_list <- lapply (vcf_list, function (item) {
        # Handle errors.
        if (class (item) == "try-error") stop (item)
        # Handle warnings.
        if (!is.null(item[[2]]))
            for (i in item[[2]])
                warning (i)

        # Unpack the GRanges
        return(item[[1]])
    })

    vcf_list <- GRangesList(vcf_list)
    # Set the provided names for the samples.
    names(vcf_list) <- sample_names

    return(vcf_list)
}

percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

options(scipen = 999)
sample_names <- unlist(strsplit(opt$parent, ","))
str_out <- c()
df_out <- data.frame()
for (sp in sample_names) {
  type_o <- as.integer(system(paste0("grep -v '^#' ", opt$sample, ":", sp, ".diff.vcf | wc -l"), intern = T))
  str_out <- append(str_out, type_o)
  system(paste0("bedtools coverage -counts -a windows.bed -b ", opt$sample, ":", sp, ".diff.vcf > ", opt$sample, ":", sp, ".diff.bed"), intern = T)
  tmp <- read.table(paste0(opt$sample, ":", sp, ".diff.bed"),
    header = F, sep = "\t", dec = ".",
    row.names = NULL, check.names = FALSE, fill = TRUE,
    comment.char = "#", stringsAsFactors = FALSE, fileEncoding = "CSGB2312")
  tmp$name <- sp
  df_out <- rbind(df_out, tmp)
}
colnames(df_out)[1:4] <- c("chrom", "chromStart", "chromEnd", "count")
library(data.table)
df_out <- as.data.table(df_out)
df_out <- df_out[df_out[, .I[which.max(count)], by = list(chrom, chromStart, chromEnd)]$V1]
for (sp in sample_names) {
  tmp <- df_out[df_out$name == sp, ][, 1:2]
  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
  colnames(tmp) <- c("#CHROM", "POS")
  tmp$ID     <- "."
  tmp$REF    <- "G"
  tmp$ALT    <- "C"
  tmp$QUAL   <- "999"
  tmp$FILTER <- "."
  tmp$INFO   <- "."
  dfo <- data.frame()
  tpos <- as.numeric(tmp$POS)
  for (i in seq(opt$window / opt$draw)) {
    tmp$POS <- tpos + (i - 1) * opt$draw * 1000
    dfo <- rbind(dfo, tmp)
  }
  write.table("##fileformat=VCFv4.1", file = paste0(opt$sample, ":", sp, ".diff.post.vcf"),
    quote = F, sep = "\t", dec = ".", append = F,
    row.names = F, col.names = F, fileEncoding = "CSGB2312")
  suppressWarnings(
  write.table(dfo, file = paste0(opt$sample, ":", sp, ".diff.post.vcf"),
    quote = F, sep = "\t", dec = ".", append = T,
    row.names = F, col.names = T, fileEncoding = "CSGB2312")
  )
}

cat(paste(c(opt$sample, str_out), collapse = "\t"), file = "plot.tsv", sep = "\n", append = TRUE)
sample_perct <- percent(str_out / sum(str_out))
sample_label <- paste(sample_names, sample_perct, sep=", ")
vcf_files <- paste0(opt$sample, ":", sample_names, ".diff.post.vcf")
genome_name <- opt$genome
groups <- unlist(strsplit(opt$chr, ","))

obj <- vcf2plot(vcf_files, sample_names, genome_name, groups)
length(obj@unlistData)
table(obj@unlistData@elementMetadata@listData$REF@ranges@group)

library(ggbio)
group <- as.factor(obj@unlistData@elementMetadata@listData$REF@ranges@group)
levels(group) <- sample_label
#p <- autoplot(obj@unlistData, layout = "karyogram", aes(color = group, fill = group), alpha = 0.2)
p <- autoplot(obj@unlistData, layout = "karyogram", aes(color = group, fill = group), alpha = 1, na.value = "white") +
  scale_color_manual(values = c("#e4f96c", "#02d842", "#2160ff")) +
  scale_fill_manual(values = c("#e4f96c", "#02d842", "#2160ff"))
ggsave(paste0(opt$sample, ".diff.png"), p@ggplot, width = 10)
ggsave(paste0(opt$sample, ".diff.pdf"), p@ggplot, width = 10)
