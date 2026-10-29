# Does sequencing depth drive the V3-V4 vs V4 separation in the Aitchison PCoA?
#
# BSI's V3-V4 libraries are about 2.6x deeper than the CRT23 V4 samples, so V4 has far more zero
# counts on the shared genera, and zeros become large negative values in the CLR. This rebuilds the
# 03_compare.R PCoA on the same 39 cow pairs and 144 shared genera three ways: raw counts (as in
# 03), every sample rarefied to the lowest depth (several draws), and proportions. If depth were
# behind the split, equal-depth versions would weaken it.
#
#   micromamba run -p ~/envs/micro_prediction Rscript region_comparison_v3v4_v4/08_depth_check.R
#
# Output: outputs/region_comparison_v3v4_v4/depth_check.csv (and a summary on stdout).
# Result (2026-10-29): PC1 26.5% raw, 27.8-28.2% rarefied (18,237 reads), 32.1% proportions; the
# primer sets split completely on PC1 in every version, so the separation is not a depth artefact.
# Within V4 alone, deeper samples sit slightly towards V3-V4 on PC1 (Spearman about -0.4).

suppressPackageStartupMessages(library(vegan))
set.seed(42)
BASE <- "outputs/region_comparison_v3v4_v4"
MIN_REL <- 1e-4
PSEUDO <- 0.5
N_DRAWS <- 5

read_genus <- function(region) {
  g <- read.delim(file.path(BASE, sprintf("genus_%s_gtdb_r214.tsv", region)), check.names = FALSE)
  m <- as.matrix(g[, -1]); rownames(m) <- g$sample
  m[, setdiff(colnames(m), "Unassigned_genus"), drop = FALSE]
}
sheet <- read.csv(file.path(BASE, "samples.csv"), colClasses = c(cow = "character", sample = "character"))
v4 <- read_genus("V4"); v34 <- read_genus("V3V4")

# same genus space, shared set and pairing as 03_compare.R
genera <- union(colnames(v4), colnames(v34))
pad <- function(m) { x <- matrix(0, nrow(m), length(genera), dimnames = list(rownames(m), genera)); x[, colnames(m)] <- m; x }
allm <- rbind(pad(v4), pad(v34))
rel <- allm / rowSums(allm)
info <- sheet[match(rownames(allm), sheet$sample), ]
shared <- genera[colMeans(rel[info$region == "V4", ]) >= MIN_REL & colMeans(rel[info$region == "V3V4", ]) >= MIN_REL]
s34 <- sheet[sheet$region == "V3V4", ]
s34 <- s34[order(-s34$reads), ]
main34 <- s34[!duplicated(s34$cow), ]
t2 <- sheet[sheet$region == "V4" & sheet$timepoint == "T2", ]
cows <- intersect(main34$cow, t2$cow)
paired <- c(setNames(main34$sample, main34$cow)[cows], setNames(t2$sample, t2$cow)[cows])
region <- rep(c("V3V4", "V4"), each = length(cows))
X <- allm[paired, shared]

cat(length(cows), "cows,", length(shared), "shared genera\n")
cat("median reads on shared genera: V3-V4", median(rowSums(X[region == "V3V4", ])),
    "V4", median(rowSums(X[region == "V4", ])), "\n")
cat(sprintf("zero cells: V3-V4 %.1f%%, V4 %.1f%%\n",
            100 * mean(X[region == "V3V4", ] == 0), 100 * mean(X[region == "V4", ] == 0)))

clr <- function(m) { l <- log(m + PSEUDO); l - rowMeans(l) }
pcoa_row <- function(m, version) {
  fit <- cmdscale(dist(clr(m)), k = 2, eig = TRUE)
  pct <- 100 * fit$eig[1:2] / sum(fit$eig[fit$eig > 0])
  pc1 <- fit$points[, 1]
  split <- max(pc1[region == "V3V4"]) < min(pc1[region == "V4"]) || max(pc1[region == "V4"]) < min(pc1[region == "V3V4"])
  data.frame(version = version, depth = if (grepl("rarefied", version)) min(rowSums(m)) else NA,
             pc1_pct = pct[1], pc2_pct = pct[2], pc1_split_complete = split,
             spearman_depth_pc1_V3V4 = cor(rowSums(X[region == "V3V4", ]), pc1[region == "V3V4"], method = "spearman"),
             spearman_depth_pc1_V4 = cor(rowSums(X[region == "V4", ]), pc1[region == "V4"], method = "spearman"))
}
depth <- min(rowSums(X))
rows <- list(pcoa_row(X, "raw counts (as 03_compare.R)"))
for (i in seq_len(N_DRAWS)) rows[[length(rows) + 1]] <- pcoa_row(rrarefy(X, depth), sprintf("rarefied, draw %d", i))
rows[[length(rows) + 1]] <- pcoa_row(X / rowSums(X) * 1e4, "proportions x 1e4")
out <- do.call(rbind, rows)
write.csv(out, file.path(BASE, "depth_check.csv"), row.names = FALSE)
print(out, digits = 3, row.names = FALSE)
