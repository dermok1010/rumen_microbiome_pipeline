# Compare V3-V4 and V4 16S genus profiles of the same dairy cows (GTDB r214 for both).
#
# Uses R with vegan + ggplot2 from ~/envs/micro_prediction (this repo's own env has no vegan):
#   micromamba run -p ~/envs/micro_prediction Rscript region_comparison_v3v4_v4/03_compare.R
#
# Steps
#   1. Pairing: BSI's sampling date is unknown, so each cow's V3-V4 sample is compared with its
#      V4 T1 and T2 samples; the timepoint that is closer for most cows is used for all cows.
#      Where a cow has two V3-V4 libraries, the one with more reads is used (the other gives a
#      technical-replicate baseline).
#   2. Genus overlap: genera each region detects; shared common set = genera with mean relative
#      abundance >= 0.01% in both regions (as in the 16S vs RE-RRS work).
#   3. Same cow vs different cows: Aitchison and Bray-Curtis distance V3-V4 -> V4 for the same
#      cow vs other cows, against two baselines: V4 T1 vs T2 within a cow (time) and V3-V4
#      replicate libraries (technical).
#   4. Variation explained: PERMANOVA (cow + region) on the paired samples.
#   5. Per genus: mean abundance in each region, and Spearman correlation across cows.
#   6. Alpha: Shannon, inverse Simpson and observed genera, rarefied to one common depth.

suppressPackageStartupMessages({ library(vegan); library(ggplot2) })
set.seed(42)
options(warn = 1)
BASE <- "outputs/region_comparison_v3v4_v4"
# CORE=1: use the genus tables of both regions trimmed to the same V4 core (04/05), outputs in core/
CORE <- Sys.getenv("CORE", "0") == "1"
OUT <- if (CORE) file.path(BASE, "core") else BASE
dir.create(OUT, showWarnings = FALSE)
MIN_REL <- 1e-4
PSEUDO <- 0.5
N_RAREFY <- 20

read_genus <- function(region) {
  g <- read.delim(file.path(BASE, sprintf(if (CORE) "genus_%s_core_gtdb_r214.tsv" else "genus_%s_gtdb_r214.tsv", region)), check.names = FALSE)
  m <- as.matrix(g[, -1]); rownames(m) <- g$sample
  m[, setdiff(colnames(m), "Unassigned_genus"), drop = FALSE]
}
sheet <- read.csv(file.path(BASE, "samples.csv"), colClasses = c(cow = "character"))
v4 <- read_genus("V4"); v34 <- read_genus("V3V4")

# one V3-V4 library per cow (most reads); the rest are replicates
s34 <- sheet[sheet$region == "V3V4", ]
s34 <- s34[order(-s34$reads), ]
main34 <- s34[!duplicated(s34$cow), ]
reps34 <- s34[duplicated(s34$cow), ]
cows <- intersect(main34$cow, sheet$cow[sheet$region == "V4"])
cat(length(cows), "cows with both regions\n")

# common genus space: union of both regions' genera, zero where a region lacks one
genera <- union(colnames(v4), colnames(v34))
pad <- function(m) { x <- matrix(0, nrow(m), length(genera), dimnames = list(rownames(m), genera)); x[, colnames(m)] <- m; x }
allm <- rbind(pad(v4), pad(v34))
rel <- allm / rowSums(allm)
info <- sheet[match(rownames(allm), sheet$sample), ]

# shared common genera
mean_v4 <- colMeans(rel[info$region == "V4", ])
mean_34 <- colMeans(rel[info$region == "V3V4", ])
shared <- genera[mean_v4 >= MIN_REL & mean_34 >= MIN_REL]
overlap <- data.frame(
  set = c("V4 genera", "V3-V4 genera", "In both (any abundance)", "Shared, >= 0.01% in both"),
  n = c(sum(colSums(v4) > 0), sum(colSums(v34) > 0),
        length(intersect(colnames(v4)[colSums(v4) > 0], colnames(v34)[colSums(v34) > 0])), length(shared)),
  pct_reads_v4 = c(NA, NA, NA, 100 * mean(rowSums(rel[info$region == "V4", shared]))),
  pct_reads_v3v4 = c(NA, NA, NA, 100 * mean(rowSums(rel[info$region == "V3V4", shared]))))
write.csv(overlap, file.path(OUT, "genus_overlap.csv"), row.names = FALSE)

clr <- function(m) { l <- log(m + PSEUDO); l - rowMeans(l) }
dist_fns <- list(Aitchison = function(m) as.matrix(dist(clr(m))),
                 `Bray-Curtis` = function(m) as.matrix(vegdist(m / rowSums(m), "bray")))

# --- 1. which V4 timepoint matches BSI's sample
X <- allm[, shared]
D <- dist_fns$Aitchison(X)
pair_rows <- lapply(cows, function(cw) {
  b <- main34$sample[main34$cow == cw]
  t1 <- sheet$sample[sheet$region == "V4" & sheet$cow == cw & sheet$timepoint == "T1"]
  t2 <- sheet$sample[sheet$region == "V4" & sheet$cow == cw & sheet$timepoint == "T2"]
  data.frame(cow = cw, d_T1 = if (length(t1)) D[b, t1] else NA, d_T2 = if (length(t2)) D[b, t2] else NA)
})
pairing <- do.call(rbind, pair_rows)
pairing$closer <- ifelse(pairing$d_T1 < pairing$d_T2, "T1", "T2")
TP <- names(which.max(table(pairing$closer)))
write.csv(pairing, file.path(OUT, "timepoint_matching.csv"), row.names = FALSE)
cat("V3-V4 closer to T1 for", sum(pairing$closer == "T1", na.rm = TRUE), "cows, T2 for",
    sum(pairing$closer == "T2", na.rm = TRUE), "-> pairing with", TP, "\n")

v4_tp <- setNames(sheet$sample[sheet$region == "V4" & sheet$timepoint == TP],
                  sheet$cow[sheet$region == "V4" & sheet$timepoint == TP])
cows <- intersect(cows, names(v4_tp))
s_v34 <- setNames(main34$sample, main34$cow)[cows]
s_v4 <- v4_tp[cows]

# --- 3. same cow vs other cows, with time and technical baselines
dist_rows <- list()
for (dn in names(dist_fns)) {
  for (gs in c("shared", "all")) {
    cols <- if (gs == "shared") shared else genera
    Dm <- dist_fns[[dn]](allm[, cols])
    same <- Dm[cbind(s_v34, s_v4)]
    other <- unlist(lapply(seq_along(cows), function(i) Dm[s_v34[i], s_v4[-i]]))
    t1 <- sheet$sample[sheet$region == "V4" & sheet$timepoint == "T1"]; names(t1) <- sheet$cow[match(t1, sheet$sample)]
    t2 <- sheet$sample[sheet$region == "V4" & sheet$timepoint == "T2"]; names(t2) <- sheet$cow[match(t2, sheet$sample)]
    both <- intersect(names(t1), names(t2))
    time <- Dm[cbind(t1[both], t2[both])]
    rep <- if (nrow(reps34)) Dm[cbind(reps34$sample, setNames(main34$sample, main34$cow)[reps34$cow])] else NA
    for (k in list(list("V3-V4 vs V4, same cow", same), list("V3-V4 vs V4, different cows", other),
                   list("V4 T1 vs T2, same cow (time)", time), list("V3-V4 replicate libraries (technical)", rep)))
      dist_rows[[length(dist_rows) + 1]] <- data.frame(distance = dn, genus_set = gs, comparison = k[[1]],
        n = sum(!is.na(k[[2]])), median = median(k[[2]], na.rm = TRUE),
        q25 = quantile(k[[2]], 0.25, na.rm = TRUE), q75 = quantile(k[[2]], 0.75, na.rm = TRUE))
    # how often a cow's own V4 sample is its V3-V4 sample's nearest V4 neighbour
    nn <- sapply(seq_along(cows), function(i) names(which.min(Dm[s_v34[i], s_v4])) == s_v4[i])
    dist_rows[[length(dist_rows) + 1]] <- data.frame(distance = dn, genus_set = gs,
      comparison = "Own V4 sample is the nearest V4 sample (%)", n = length(nn),
      median = 100 * mean(nn), q25 = NA, q75 = NA)
  }
}
dists <- do.call(rbind, dist_rows)
write.csv(dists, file.path(OUT, "distances.csv"), row.names = FALSE)

# --- 4. PERMANOVA on the paired samples: cow vs region
pm_rows <- list()
paired <- c(s_v34, s_v4)
pdat <- data.frame(cow = c(cows, cows), region = rep(c("V3V4", "V4"), each = length(cows)))
for (dn in names(dist_fns)) {
  Dm <- as.dist(dist_fns[[dn]](allm[paired, shared]))
  a <- adonis2(Dm ~ cow + region, data = pdat, by = "terms", permutations = 999)
  pm_rows[[length(pm_rows) + 1]] <- data.frame(distance = dn, term = rownames(a)[1:3], R2 = a$R2[1:3], p = a$`Pr(>F)`[1:3])
}
write.csv(do.call(rbind, pm_rows), file.path(OUT, "permanova.csv"), row.names = FALSE)

# --- 5. per genus
rv34 <- rel[s_v34, shared]; rv4 <- rel[s_v4, shared]
per_genus <- data.frame(genus = shared, mean_pct_V4 = 100 * colMeans(rv4), mean_pct_V3V4 = 100 * colMeans(rv34),
  spearman_across_cows = sapply(shared, function(g) suppressWarnings(cor(rv34[, g], rv4[, g], method = "spearman"))))
per_genus$log2_ratio_V3V4_over_V4 <- log2(per_genus$mean_pct_V3V4 / per_genus$mean_pct_V4)
per_genus <- per_genus[order(-per_genus$mean_pct_V4), ]
write.csv(per_genus, file.path(OUT, "per_genus.csv"), row.names = FALSE)

# --- 6. alpha at one common depth
alpha_rows <- list()
for (gs in c("shared", "all")) {
  cols <- if (gs == "shared") shared else genera
  m <- allm[paired, cols]; depth <- min(rowSums(m))
  reps <- replicate(N_RAREFY, { r <- rrarefy(m, depth)
    cbind(observed = rowSums(r > 0), shannon = diversity(r), inv_simpson = diversity(r, "invsimpson")) },
    simplify = "array")
  a <- apply(reps, c(1, 2), mean)
  for (mt in colnames(a)) {
    x34 <- a[s_v34, mt]; x4 <- a[s_v4, mt]
    w <- suppressWarnings(wilcox.test(x34, x4, paired = TRUE))
    alpha_rows[[length(alpha_rows) + 1]] <- data.frame(genus_set = gs, rarefied_to = depth, metric = mt,
      mean_V4 = mean(x4), mean_V3V4 = mean(x34), mean_diff = mean(x34 - x4),
      paired_wilcoxon_p = w$p.value, spearman = cor(x34, x4, method = "spearman"))
  }
}
write.csv(do.call(rbind, alpha_rows), file.path(OUT, "alpha.csv"), row.names = FALSE)

# --- plots
theme_set(theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank()))
for (dn in names(dist_fns)) {
  Dm <- dist_fns[[dn]](allm[paired, shared])
  fit <- cmdscale(as.dist(Dm), k = 2, eig = TRUE); pct <- 100 * fit$eig[1:2] / sum(fit$eig[fit$eig > 0])
  pts <- data.frame(PC1 = fit$points[, 1], PC2 = fit$points[, 2], cow = pdat$cow,
                    region = ifelse(pdat$region == "V4", sprintf("V4 (%s)", TP), "V3-V4"))
  seg <- merge(pts[pts$region == "V3-V4", ], pts[pts$region != "V3-V4", ], by = "cow")
  g <- ggplot(pts, aes(PC1, PC2)) +
    geom_segment(data = seg, aes(x = PC1.x, y = PC2.x, xend = PC1.y, yend = PC2.y), colour = "grey70") +
    geom_point(aes(colour = region), size = 2.2) +
    scale_colour_manual(values = setNames(c("#E8954A", "#2F6690"), c("V3-V4", sprintf("V4 (%s)", TP)))) +
    labs(title = sprintf("Same cows, V3-V4%s vs V4: %s PCoA (%d shared genera)", if (CORE) " trimmed to V4" else "", dn, length(shared)),
         subtitle = "Grey lines join the two samples of each cow",
         x = sprintf("PCo1 (%.1f%%)", pct[1]), y = sprintf("PCo2 (%.1f%%)", pct[2]), colour = NULL)
  ggsave(file.path(OUT, sprintf("pcoa_%s.png", tolower(gsub("-", "", dn)))), g, width = 8, height = 6, dpi = 200)
}
pg <- per_genus
g <- ggplot(pg, aes(mean_pct_V4, mean_pct_V3V4)) + geom_abline(linetype = 2, colour = "grey50") +
  geom_point(aes(colour = spearman_across_cows), size = 2) + scale_x_log10() + scale_y_log10() +
  scale_colour_gradient2(low = "#B85A24", mid = "grey85", high = "#2F6690", midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Mean genus abundance, V3-V4 vs V4 (same cows)", x = "V4: mean % of reads (log)",
       y = "V3-V4: mean % of reads (log)", colour = "Spearman\nacross cows")
ggsave(file.path(OUT, "per_genus_abundance.png"), g, width = 7.5, height = 6, dpi = 200)
dd <- dists[dists$genus_set == "shared" & !grepl("nearest", dists$comparison), ]
g <- ggplot(dd, aes(median, comparison)) + geom_errorbarh(aes(xmin = q25, xmax = q75), height = 0.2) +
  geom_point(size = 2.5, colour = "#2F6690") + facet_wrap(~distance, scales = "free_x") +
  labs(title = "Distances between samples (shared genera): median and IQR", x = "Distance", y = NULL)
ggsave(file.path(OUT, "distance_baselines.png"), g, width = 10, height = 3.8, dpi = 200)
cat("done\n")
