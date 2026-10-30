# eSNP-Karyotyping MajorMinorCalc + PlotGenome on our SNP calls. Rscript scripts/bench_esnp_plot.R
suppressMessages({library(eSNPKaryotyping); library(zoo); library(gplots)})
in_dir <- "results/benchmark/esnp_input"; out_dir <- "results/benchmark/esnp_plots"; mm_dir <- "results/benchmark/esnp_mm"
for (d in c(out_dir, mm_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
for (f in sort(list.files(in_dir, pattern = "\\.tsv$", full.names = TRUE))) {
  id <- sub("\\.tsv$", "", basename(f))
  tab <- read.delim(f)
  tab$chr <- as.numeric(tab$chr)
  tab <- tab[order(tab$chr, tab$position), ]
  tab <- tab[tab$chr > 0, ]
  mm <- MajorMinorCalc(Table = tab, minDP = 20, maxDP = 1e9, minAF = 0.2)
  write.csv(mm, file.path(mm_dir, paste0(id, ".csv")), row.names = FALSE)
  png(file.path(out_dir, paste0(id, ".png")), width = 1800, height = 700, res = 120)
  PlotGenome(orderedTable = mm, Window = 151, Ylim = 3, PValue = TRUE, Organism = "Human")
  dev.off()
  cat(id, "done\n")
}
