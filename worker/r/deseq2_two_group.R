#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx == length(args)) {
    stop(sprintf("Missing required argument: %s", flag))
  }
  args[idx + 1]
}
parse_arg_default <- function(flag, default_value) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx == length(args)) {
    return(default_value)
  }
  args[idx + 1]
}

count_matrix <- parse_arg("--count-matrix")
group_a <- parse_arg("--group-a")
group_b <- parse_arg("--group-b")
out_dir <- parse_arg("--out-dir")
label_sig <- tolower(parse_arg_default("--label-significant-genes", "false")) %in% c("1", "true", "yes")
max_labels <- as.integer(parse_arg_default("--max-labels", "25"))
file_prefix <- parse_arg_default("--file-prefix", "")
run_gsea <- tolower(parse_arg_default("--run-gsea", "true")) %in% c("1", "true", "yes")
prepare_cibersortx <- tolower(parse_arg_default("--prepare-cibersortx", "true")) %in% c("1", "true", "yes")
if (is.na(max_labels) || max_labels < 0) max_labels <- 25

if (nchar(file_prefix) > 0) {
  file_prefix <- gsub("[^A-Za-z0-9._-]+", "_", file_prefix)
  file_prefix <- gsub("_+", "_", file_prefix)
  file_prefix <- gsub("^[_\\.-]+|[_\\.-]+$", "", file_prefix)
}
if (nchar(file_prefix) == 0) file_prefix <- "project"
outname <- function(base) paste0(file_prefix, "_", base)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

x <- read.csv(count_matrix, check.names = FALSE)
meta_cols <- c("sample_id", "file_name", "group")
missing_meta <- setdiff(meta_cols, colnames(x))
if (length(missing_meta) > 0) {
  stop(paste("Count matrix missing metadata columns:", paste(missing_meta, collapse = ", ")))
}

x <- x[x$group %in% c(group_a, group_b), ]
if (nrow(x) == 0) stop("No samples after group filter")

gene_cols <- setdiff(colnames(x), meta_cols)
counts <- as.matrix(round(t(x[, gene_cols, drop = FALSE])))
colnames(counts) <- x$sample_id

coldata <- data.frame(
  row.names = x$sample_id,
  group = factor(x$group, levels = c(group_a, group_b))
)

if (sum(coldata$group == group_a) < 2 || sum(coldata$group == group_b) < 2) {
  stop("Need at least two samples per group")
}

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ group)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", group_b, group_a))
res <- res[order(res$pvalue), ]
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]

sig <- !is.na(res_df$padj) & res_df$padj < 0.05
write.csv(res_df, file.path(out_dir, outname("deseq2_all_genes.csv")), row.names = FALSE)
write.csv(res_df[sig, ], file.path(out_dir, outname("deseq2_significant_genes.csv")), row.names = FALSE)

# QC table
qc <- data.frame(
  metric = c(
    "samples_total",
    paste0("samples_", group_a),
    paste0("samples_", group_b),
    "genes_before_filter",
    "genes_after_filter",
    "significant_padj_lt_0.05"
  ),
  value = c(
    nrow(coldata),
    sum(coldata$group == group_a),
    sum(coldata$group == group_b),
    nrow(counts),
    nrow(dds),
    sum(sig, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
write.csv(qc, file.path(out_dir, outname("qc_sample_counts.csv")), row.names = FALSE)

# Top genes tables
res_sig <- res_df[!is.na(res_df$padj), , drop = FALSE]
up <- res_sig[order(-res_sig$log2FoldChange, res_sig$padj), , drop = FALSE]
down <- res_sig[order(res_sig$log2FoldChange, res_sig$padj), , drop = FALSE]
if (nrow(up) > 50) up <- up[seq_len(50), , drop = FALSE]
if (nrow(down) > 50) down <- down[seq_len(50), , drop = FALSE]
write.csv(up, file.path(out_dir, outname("top_genes_up.csv")), row.names = FALSE)
write.csv(down, file.path(out_dir, outname("top_genes_down.csv")), row.names = FALSE)

# Run parameters JSON
json_path <- file.path(out_dir, outname("run_parameters.json"))
json_lines <- c(
  "{",
  paste0('  "file_prefix": "', file_prefix, '",'),
  paste0('  "group_a": "', group_a, '",'),
  paste0('  "group_b": "', group_b, '",'),
  paste0('  "label_significant_genes": ', ifelse(label_sig, "true", "false"), ","),
  paste0('  "max_labels": ', max_labels, ","),
  paste0('  "run_gsea": ', ifelse(run_gsea, "true", "false"), ","),
  paste0('  "prepare_cibersortx": ', ifelse(prepare_cibersortx, "true", "false"), ","),
  paste0('  "samples_total": ', nrow(coldata), ","),
  paste0('  "samples_', group_a, '": ', sum(coldata$group == group_a), ","),
  paste0('  "samples_', group_b, '": ', sum(coldata$group == group_b), ","),
  paste0('  "genes_before_filter": ', nrow(counts), ","),
  paste0('  "genes_after_filter": ', nrow(dds), ","),
  paste0('  "significant_padj_lt_0.05": ', sum(sig, na.rm = TRUE)),
  "}"
)
writeLines(json_lines, con = json_path)

# Volcano
plot_df <- res_df[!is.na(res_df$log2FoldChange) & !is.na(res_df$pvalue), ]
plot_df$neglog10p <- -log10(pmax(plot_df$pvalue, .Machine$double.xmin))
plot_df$sig <- ifelse(!is.na(plot_df$padj) & plot_df$padj < 0.05, "padj < 0.05", "ns")

p <- ggplot(plot_df, aes(x = log2FoldChange, y = neglog10p, color = sig)) +
  geom_point(alpha = 0.35, size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("padj < 0.05" = "#d62728", "ns" = "gray50")) +
  theme_bw(base_size = 11) +
  labs(
    title = sprintf("DESeq2: %s vs %s", group_b, group_a),
    subtitle = "positive log2FC = higher in group_b",
    x = sprintf("log2FC (%s / %s)", group_b, group_a),
    y = expression(-log[10](italic(p))),
    color = NULL
  )

if (label_sig && nrow(plot_df) > 0) {
  label_df <- plot_df[plot_df$sig == "padj < 0.05", , drop = FALSE]
  if (nrow(label_df) > 0) {
    label_df <- label_df[order(-label_df$neglog10p), , drop = FALSE]
    if (max_labels > 0 && nrow(label_df) > max_labels) {
      label_df <- label_df[seq_len(max_labels), , drop = FALSE]
    }
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = label_df,
        aes(label = gene_id),
        size = 2.7,
        max.overlaps = Inf,
        show.legend = FALSE
      )
    } else {
      p <- p + geom_text(
        data = label_df,
        aes(label = gene_id),
        size = 2.2,
        check_overlap = TRUE,
        show.legend = FALSE
      )
    }
  }
}

ggsave(file.path(out_dir, outname("volcano.png")), p, width = 8, height = 6, dpi = 200)
ggsave(file.path(out_dir, outname("volcano.pdf")), p, width = 8, height = 6)

# MA plot
png(file.path(out_dir, outname("ma_plot.png")), width = 1400, height = 1000, res = 160)
plotMA(res, main = sprintf("MA plot: %s vs %s", group_b, group_a), ylim = c(-5, 5))
dev.off()
pdf(file.path(out_dir, outname("ma_plot.pdf")), width = 8, height = 6)
plotMA(res, main = sprintf("MA plot: %s vs %s", group_b, group_a), ylim = c(-5, 5))
dev.off()

# CIBERSORTx prep
if (prepare_cibersortx) {
  norm_counts <- counts(dds, normalized = TRUE)
  mix <- as.data.frame(norm_counts)
  mix$gene_id <- rownames(mix)

  # Convert ENSG to symbols when possible
  mix$gene_symbol <- mix$gene_id
  if (requireNamespace("AnnotationDbi", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    ens_base <- sub("\\..*", "", mix$gene_id)
    # Only attempt human ENSG mapping when IDs look Ensembl-human-like.
    if (sum(grepl("^ENSG", ens_base)) >= max(50, floor(length(ens_base) * 0.1))) {
      syms <- tryCatch(
        AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db,
          keys = unique(ens_base),
          keytype = "ENSEMBL",
          column = "SYMBOL",
          multiVals = "first"
        ),
        error = function(e) NULL
      )
      if (!is.null(syms)) {
        mix$gene_symbol <- unname(syms[ens_base])
        mix$gene_symbol[is.na(mix$gene_symbol) | mix$gene_symbol == ""] <- mix$gene_id[is.na(mix$gene_symbol) | mix$gene_symbol == ""]
      }
    }
  }

  # Aggregate duplicate symbols
  sample_names <- colnames(norm_counts)
  agg <- aggregate(mix[, sample_names, drop = FALSE], by = list(gene_symbol = mix$gene_symbol), FUN = sum)
  agg <- agg[agg$gene_symbol != "" & !is.na(agg$gene_symbol), , drop = FALSE]
  write.csv(agg, file.path(out_dir, outname("cibersortx_mixture_normalized_counts.csv")), row.names = FALSE)

  pheno <- data.frame(sample_id = rownames(coldata), group = coldata$group, stringsAsFactors = FALSE)
  write.csv(pheno, file.path(out_dir, outname("cibersortx_sample_metadata.csv")), row.names = FALSE)

  cib_txt <- c(
    "CIBERSORTx input package",
    "",
    "Files:",
    paste0("- ", outname("cibersortx_mixture_normalized_counts.csv"), " (gene_symbol x sample matrix, linear scale)"),
    paste0("- ", outname("cibersortx_sample_metadata.csv")),
    "",
    "Suggested steps:",
    "1) Upload mixture matrix to CIBERSORTx.",
    "2) Use your selected signature matrix/mode.",
    "3) Join output fractions back to sample metadata for downstream comparison."
  )
  writeLines(cib_txt, con = file.path(out_dir, outname("cibersortx_instructions.txt")))
}

# GSEA
if (run_gsea) {
  gsea_status <- c()
  tryCatch({
    if (requireNamespace("fgsea", quietly = TRUE) && requireNamespace("msigdbr", quietly = TRUE)) {
    ranks <- res_df$stat
    names(ranks) <- sub("\\..*", "", res_df$gene_id)
    ranks <- ranks[is.finite(ranks)]
    ranks <- sort(ranks, decreasing = TRUE)

    collapse_ranks <- function(x) {
      x <- x[is.finite(x)]
      x <- x[!is.na(names(x)) & names(x) != ""]
      ord <- order(-abs(x))
      x <- x[ord]
      x <- x[!duplicated(names(x))]
      sort(x, decreasing = TRUE)
    }

    run_one <- function(species, category, subcat, out_stub) {
      msig <- msigdbr::msigdbr(species = species, category = category, subcategory = subcat)
      pathways <- split(msig$gene_symbol, msig$gs_name)
      fg <- fgsea::fgseaMultilevel(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500)
      fg <- as.data.frame(fg)
      if (nrow(fg) > 0) {
        # fgsea returns list columns (e.g. leadingEdge); serialize for CSV output.
        for (cn in names(fg)) {
          if (is.list(fg[[cn]])) {
            fg[[cn]] <- vapply(
              fg[[cn]],
              function(x) paste(as.character(x), collapse = ";"),
              character(1)
            )
          }
        }
      }
      if (nrow(fg) == 0) {
        write.csv(fg, file.path(out_dir, outname(paste0("gsea_", out_stub, ".csv"))), row.names = FALSE)
        return(fg)
      }
      fg <- fg[order(fg$padj, -abs(fg$NES)), ]
      write.csv(fg, file.path(out_dir, outname(paste0("gsea_", out_stub, ".csv"))), row.names = FALSE)

      top <- head(fg, 20)
      top$pathway <- factor(top$pathway, levels = rev(top$pathway))
      gp <- ggplot(top, aes(x = pathway, y = NES, fill = NES > 0)) +
        geom_col() +
        coord_flip() +
        theme_bw(base_size = 10) +
        scale_fill_manual(values = c("TRUE" = "#d62728", "FALSE" = "#1f77b4"), guide = "none") +
        labs(title = paste("GSEA", out_stub), x = NULL, y = "NES")
      ggsave(file.path(out_dir, outname(paste0("gsea_", out_stub, "_top.png"))), gp, width = 9, height = 6, dpi = 180)
      ggsave(file.path(out_dir, outname(paste0("gsea_", out_stub, "_top.pdf"))), gp, width = 9, height = 6)
      fg
    }

    # Need symbols for msigdbr pathways
    ens <- names(ranks)
    if (sum(grepl("^ENSG", ens)) >= max(50, floor(length(ens) * 0.1))) {
      # ENSEMBL-like input: map to HGNC symbols first.
      if (requireNamespace("AnnotationDbi", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        sym <- tryCatch(
          AnnotationDbi::mapIds(
            org.Hs.eg.db::org.Hs.eg.db,
            keys = unique(ens),
            keytype = "ENSEMBL",
            column = "SYMBOL",
            multiVals = "first"
          ),
          error = function(e) NULL
        )
        if (!is.null(sym)) {
          sym_names <- unname(sym[ens])
          keep <- !is.na(sym_names) & sym_names != ""
          ranks2 <- ranks[keep]
          names(ranks2) <- sym_names[keep]
          ranks <- collapse_ranks(ranks2)

          g_h <- run_one("Homo sapiens", "H", NULL, "hallmark")
          g_k <- run_one("Homo sapiens", "C2", "CP:KEGG_LEGACY", "kegg")
          g_r <- run_one("Homo sapiens", "C2", "CP:REACTOME", "reactome")

          summary_rows <- function(df, coll) {
            if (nrow(df) == 0) return(data.frame(collection=character(), pathway=character(), NES=numeric(), padj=numeric(), stringsAsFactors = FALSE))
            top <- head(df[order(df$padj), c("pathway", "NES", "padj")], 10)
            data.frame(collection = coll, pathway = top$pathway, NES = top$NES, padj = top$padj, stringsAsFactors = FALSE)
          }
          gsum <- rbind(summary_rows(g_h, "hallmark"), summary_rows(g_k, "kegg"), summary_rows(g_r, "reactome"))
          write.csv(gsum, file.path(out_dir, outname("gsea_summary_top10.csv")), row.names = FALSE)
          gsea_status <- c(gsea_status, "GSEA complete")
        } else {
          gsea_status <- c(gsea_status, "GSEA skipped: failed ENSEMBL->SYMBOL mapping")
        }
      } else {
        gsea_status <- c(gsea_status, "GSEA skipped: org.Hs.eg.db/AnnotationDbi not available for ENSEMBL mapping")
      }
    } else {
      # Symbol-like input: run directly with gene symbols.
      ranks <- collapse_ranks(ranks)
      g_h <- run_one("Homo sapiens", "H", NULL, "hallmark")
      g_k <- run_one("Homo sapiens", "C2", "CP:KEGG_LEGACY", "kegg")
      g_r <- run_one("Homo sapiens", "C2", "CP:REACTOME", "reactome")

      summary_rows <- function(df, coll) {
        if (nrow(df) == 0) return(data.frame(collection=character(), pathway=character(), NES=numeric(), padj=numeric(), stringsAsFactors = FALSE))
        top <- head(df[order(df$padj), c("pathway", "NES", "padj")], 10)
        data.frame(collection = coll, pathway = top$pathway, NES = top$NES, padj = top$padj, stringsAsFactors = FALSE)
      }
      gsum <- rbind(summary_rows(g_h, "hallmark"), summary_rows(g_k, "kegg"), summary_rows(g_r, "reactome"))
      write.csv(gsum, file.path(out_dir, outname("gsea_summary_top10.csv")), row.names = FALSE)
      gsea_status <- c(gsea_status, "GSEA complete (symbol input)")
    }
    } else {
      gsea_status <- c(gsea_status, "GSEA skipped: fgsea and/or msigdbr not installed")
    }
  }, error = function(e) {
    gsea_status <<- c(gsea_status, paste0("GSEA failed: ", conditionMessage(e)))
  })
  writeLines(gsea_status, con = file.path(out_dir, outname("gsea_status.txt")))
}

summary_path <- file.path(out_dir, outname("summary.txt"))
sink(summary_path)
cat("Two-group DESeq2 run\n")
cat(sprintf("group_a=%s, group_b=%s\n", group_a, group_b))
cat(sprintf("file_prefix=%s\n", file_prefix))
cat(sprintf("run_gsea=%s\n", run_gsea))
cat(sprintf("prepare_cibersortx=%s\n", prepare_cibersortx))
print(table(coldata$group))
cat(sprintf("genes after filter=%d\n", nrow(dds)))
cat(sprintf("significant padj<0.05=%d\n", sum(sig, na.rm = TRUE)))
summary(res)
sink()

# Bundle all artifacts into one zip if zip is available
zip_bin <- Sys.which("zip")
if (nzchar(zip_bin)) {
  old <- getwd()
  setwd(out_dir)
  files_now <- list.files(".", all.files = FALSE, no.. = TRUE)
  bundle_name <- outname("results_bundle.zip")
  files_to_zip <- files_now[files_now != bundle_name]
  if (length(files_to_zip) > 0) {
    utils::zip(zipfile = bundle_name, files = files_to_zip)
  }
  setwd(old)
}
