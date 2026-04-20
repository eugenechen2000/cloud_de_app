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

count_matrix <- parse_arg("--count-matrix")
group_a <- parse_arg("--group-a")
group_b <- parse_arg("--group-b")
out_dir <- parse_arg("--out-dir")

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
write.csv(res_df, file.path(out_dir, "deseq2_all_genes.csv"), row.names = FALSE)
write.csv(res_df[sig, ], file.path(out_dir, "deseq2_significant_genes.csv"), row.names = FALSE)

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

ggsave(file.path(out_dir, "volcano.png"), p, width = 8, height = 6, dpi = 200)
ggsave(file.path(out_dir, "volcano.pdf"), p, width = 8, height = 6)

sink(file.path(out_dir, "summary.txt"))
cat("Two-group DESeq2 run\n")
cat(sprintf("group_a=%s, group_b=%s\n", group_a, group_b))
print(table(coldata$group))
cat(sprintf("genes after filter=%d\n", nrow(dds)))
cat(sprintf("significant padj<0.05=%d\n", sum(sig, na.rm = TRUE)))
summary(res)
sink()
