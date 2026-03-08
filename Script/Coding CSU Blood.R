#Modul: Analisis Ekspresi Gen Pasien Chronic Spontaneous Urticaria (CSU) Melalui Sampel Darah
#Dataset: GSE72541 (CSU vs Healthy Control)
#Platform: Agilent SurePrint G3 Human GE v2 - GPL16699
#Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG)

#PART A. PENGANTAR KONSEP 
#Analisis ini membandingkan profil transkriptomik darah pasien CSU dengan individu sehat.
#Menggunakan model linear limma untuk stabilitas statistik pada data microarray.

#PART B. PERSIAPAN LINGKUNGAN KERJA (INSTALL & LOAD PACKAGE) 

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Paket Bioconductor & CRAN
bioc_pkgs <- c("GEOquery", "limma", "AnnotationDbi")
cran_pkgs <- c("pheatmap", "ggplot2", "dplyr", "umap")

for (pkg in bioc_pkgs) {
  if (!require(pkg, character.only = TRUE)) BiocManager::install(pkg, update = FALSE)
}
for (pkg in cran_pkgs) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
}

library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(umap)

#PART C. PENGAMBILAN DATA DARI GEO 
#AnnotGPL = TRUE untuk mendapatkan simbol gen langsung dari metadata GEO.
gset <- getGEO("GSE72541", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

#PART D. PRE-PROCESSING DATA EKSPRESI 
ex <- exprs(gset)

#Log2 Transformasi
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

#PART E. DEFINISI KELOMPOK SAMPEL 
group_info <- pData(gset)[["characteristics_ch1"]]

# Membersihkan teks: "diagnosis: CSU" menjadi "CSU"
groups <- make.names(gsub("diagnosis: ", "", group_info))
gset$group <- factor(groups)

levels(gset$group) <- c("CSU", "Control")

#PART F. DESIGN MATRIX & CONTRAST 
design <- model.matrix(~0 + gset$group)
colnames(design) <- levels(gset$group)

# Menentukan perbandingan: Kasus (CSU) - Kontrol (healthy.control)
grup_kasus <- levels(gset$group)[1] 
grup_kontrol <- levels(gset$group)[2]
contrast_formula <- paste(grup_kasus, "-", grup_kontrol, sep="")

contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

#PART G. ANALISIS DIFFERENTIAL EXPRESSION (LIMMA)
fit <- lmFit(ex, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Mengambil hasil DEG (p-value < 0.05)
topTableResults <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf, p.value = 1)

head(topTableResults)

# PART H. ANOTASI NAMA GEN (REVISI TIPE DATA & KONFLIK SELECT) 

# Menyiapkan data fitur dan ID Probe dari hasil limma
fdata <- fData(gset) 
topTableResults$PROBEID <- rownames(topTableResults)

# Deteksi kolom Simbol dan Deskripsi/Arti Gen secara dinamis
# Mencari kolom yang mengandung kata 'symbol' dan 'title/definition/description'
col_symbol <- colnames(fdata)[grep("symbol", colnames(fdata), ignore.case = TRUE)][1]
col_descr  <- colnames(fdata)[grep("title|definition|description", colnames(fdata), ignore.case = TRUE)][1]

# Fallback ke ID jika kolom tidak ditemukan
if (is.na(col_symbol)) col_symbol <- "ID"
if (is.na(col_descr))  col_descr  <- "ID"

# Membuat tabel pemetaan (Mapping Table)
# Menggunakan dplyr secara spesifik untuk menghindari konflik namespace
gene_mapping <- fdata %>%
  dplyr::select(
    ID, 
    SYMBOL = !!sym(col_symbol), 
    GENE_NAME = !!sym(col_descr)
  ) %>%
  # Konversi ID ke Character untuk sinkronisasi tipe data
  dplyr::mutate(ID = as.character(ID)) %>% 
  dplyr::filter(ID %in% topTableResults$PROBEID)

# Penggabungan (Merging) ke Tabel Hasil Utama
topTableResults <- topTableResults %>%
  dplyr::mutate(PROBEID = as.character(PROBEID)) %>%
  dplyr::inner_join(gene_mapping, by = c("PROBEID" = "ID"))

# Pembersihan Data (Filtering)
# Membuang transkrip XLOC, nilai kosong, atau NA pada simbol gen
topTableResults <- topTableResults %>%
  dplyr::filter(!grepl("XLOC", SYMBOL)) %>% 
  dplyr::filter(SYMBOL != "" & !is.na(SYMBOL))

# Verifikasi Akhir
# Mengecek apakah kolom SYMBOL dan GENE_NAME sudah muncul dengan benar
head(topTableResults[, c("PROBEID", "SYMBOL", "GENE_NAME")])

#PART I.1 BOXPLOT DISTRIBUSI 
boxplot(ex, col = as.numeric(gset$group), las = 2, outline = FALSE, 
        main = "Boxplot Distribusi Ekspresi GSE72541", 
        ylab = "Expression Value (log2)")

legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.8
)

#PART I.2 DENSITY PLOT 
expr_long <- data.frame(
  Expression = as.vector(ex), 
  Group = rep(gset$group, each = nrow(ex)))
ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) + theme_minimal() +
  labs(title = "Density Plot: CSU vs Healthy", 
       x = "Expression Value (log2)",
       y = "Density")

#PART I.3 UMAP 
umap_input <- t(na.omit(ex))
umap_result <- umap(umap_input)
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1], 
  UMAP2 = umap_result$layout[, 2], 
  Group = gset$group)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) + theme_minimal() +
  labs(title = "UMAP Plot: GSE72541 CSU Cluster")

#PART J.1 VISUALISASI VOLCANO PLOT 
topTableResults$status <- "NO"
topTableResults$status[topTableResults$logFC > 0.1 & topTableResults$adj.P.Val < 0.05] <- "UP"
topTableResults$status[topTableResults$logFC < -0.1 & topTableResults$adj.P.Val < 0.05] <- "DOWN"

ggplot(topTableResults, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() + ggtitle("Volcano Plot: CSU vs Healthy Control")

# PART J.2 VISUALISASI HEATMAP

topTableResults <- topTableResults[order(topTableResults$adj.P.Val), ]

topTable_Unique <- topTableResults %>%
  # Hapus baris yang tidak memiliki Symbol Gen (opsional, agar heatmap bersih)
  filter(!is.na(SYMBOL) & SYMBOL != "") %>% 
  # Pilih baris pertama untuk setiap Symbol (karena sudah diurutkan, ini adalah p-value terkecil)
  distinct(SYMBOL, .keep_all = TRUE)

top50 <- head(topTable_Unique, 50)

mat_heatmap <- ex[top50$PROBEID, ]

rownames(mat_heatmap) <- top50$SYMBOL

mat_heatmap <- mat_heatmap[complete.cases(mat_heatmap), ]
mat_heatmap <- mat_heatmap[apply(mat_heatmap, 1, var) > 0, ]

pheatmap(
  mat_heatmap,
  scale = "row",
  annotation_col = data.frame(Group = gset$group, row.names = colnames(ex)),
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,              
  annotation_legend = TRUE,      
  legend = TRUE,                 
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  main = "Top 50 Differentially Expressed Genes: CSU vs Control"
)

#PART K. MENYIMPAN HASIL 
write.csv(topTableResults, "Hasil_GSE72541_DEG.csv")
message("Analisis selesai. File hasil telah disimpan.")
