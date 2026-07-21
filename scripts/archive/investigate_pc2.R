# investigate_pc2.R
# Investigate what PC2 represents in chr21-only PCA

library(tidyverse)
library(DESeq2)

# Load data
count_data <- read_csv("data/processed/count_matrix.csv", show_col_types = FALSE)
metadata <- read_csv("data/processed/sample_metadata.csv", show_col_types = FALSE)
gene_annotations <- read_csv("data/processed/gene_annotations.csv", show_col_types = FALSE)

# Get count matrix
count_matrix <- count_data %>%
  select(-Gene_name, -Chr) %>%
  column_to_rownames("EnsemblID") %>%
  as.matrix()

# Round to integers
count_matrix <- round(count_matrix)

# Get chr21 genes only
chr21_genes <- gene_annotations %>%
  filter(Chr == "chr21") %>%
  pull(EnsemblID)

chr21_count_matrix <- count_matrix[chr21_genes, , drop = FALSE]

# Prepare colData
metadata <- metadata %>%
  filter(LabID %in% colnames(chr21_count_matrix)) %>%
  arrange(match(LabID, colnames(chr21_count_matrix)))

col_data <- DataFrame(
  sample_id = metadata$LabID,
  karyotype = factor(metadata$Karyotype, levels = c("Control", "T21")),
  sex = factor(metadata$Sex),
  age = metadata$Age_at_visit,
  bmi = metadata$BMI,
  visit = factor(metadata$Event_name)
)
rownames(col_data) <- metadata$LabID

# Create DESeqDataSet
dds_chr21 <- DESeqDataSetFromMatrix(
  countData = chr21_count_matrix,
  colData = col_data,
  design = ~ karyotype
)

dds_chr21 <- dds_chr21[rowSums(counts(dds_chr21)) > 10, ]

# Variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds_chr21, blind = FALSE)

# Get PCA data
pca_data <- plotPCA(vsd, intgroup = c("karyotype"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

# Add metadata to PCA data
pca_data$sex <- col_data$sex[match(rownames(pca_data), rownames(col_data))]
pca_data$age <- col_data$age[match(rownames(pca_data), rownames(col_data))]
pca_data$bmi <- col_data$bmi[match(rownames(pca_data), rownames(col_data))]
pca_data$visit <- col_data$visit[match(rownames(pca_data), rownames(col_data))]

cat("=== PC1 vs PC2 Analysis ===\n\n")

# Test PC1 associations
cat("PC1 associations:\n")
cat(sprintf("  Karyotype (T21 vs Control) p-value: %.2e\n",
            t.test(PC1 ~ karyotype, data = pca_data)$p.value))
cat(sprintf("  Sex (Male vs Female) p-value: %.2e\n",
            t.test(PC1 ~ sex, data = pca_data)$p.value))
cat(sprintf("  Age correlation: r=%.3f, p=%.2e\n",
            cor(pca_data$PC1, pca_data$age, use = "complete.obs"),
            cor.test(pca_data$PC1, pca_data$age)$p.value))
cat(sprintf("  BMI correlation: r=%.3f, p=%.2e\n",
            cor(pca_data$PC1, pca_data$bmi, use = "complete.obs"),
            cor.test(pca_data$PC1, pca_data$bmi)$p.value))

cat("\nPC2 associations:\n")
cat(sprintf("  Karyotype (T21 vs Control) p-value: %.2e\n",
            t.test(PC2 ~ karyotype, data = pca_data)$p.value))
cat(sprintf("  Sex (Male vs Female) p-value: %.2e\n",
            t.test(PC2 ~ sex, data = pca_data)$p.value))
cat(sprintf("  Age correlation: r=%.3f, p=%.2e\n",
            cor(pca_data$PC2, pca_data$age, use = "complete.obs"),
            cor.test(pca_data$PC2, pca_data$age)$p.value))
cat(sprintf("  BMI correlation: r=%.3f, p=%.2e\n",
            cor(pca_data$PC2, pca_data$bmi, use = "complete.obs"),
            cor.test(pca_data$PC2, pca_data$bmi)$p.value))

# Check visit number
if (length(unique(pca_data$visit)) > 1) {
  visit_anova_pc1 <- anova(lm(PC1 ~ visit, data = pca_data))
  visit_anova_pc2 <- anova(lm(PC2 ~ visit, data = pca_data))
  cat(sprintf("\n  Visit number (PC1) p-value: %.2e\n", visit_anova_pc1$`Pr(>F)`[1]))
  cat(sprintf("  Visit number (PC2) p-value: %.2e\n", visit_anova_pc2$`Pr(>F)`[1]))
}

# Create detailed PCA plots
pdf("results/figures/pca_chr21_detailed.pdf", width = 14, height = 10)
par(mfrow = c(2, 3))

# PC1 vs PC2 by karyotype
plot(pca_data$PC1, pca_data$PC2,
     col = ifelse(pca_data$karyotype == "T21", "red", "blue"),
     pch = 19, cex = 1.5, main = "Chr21 PCA: by Karyotype",
     xlab = paste0("PC1: ", percent_var[1], "% variance"),
     ylab = paste0("PC2: ", percent_var[2], "% variance"))
legend("topright", legend = c("Control", "T21"),
       col = c("blue", "red"), pch = 19)

# PC1 vs PC2 by sex
plot(pca_data$PC1, pca_data$PC2,
     col = ifelse(pca_data$sex == "Male", "green", "purple"),
     pch = 19, cex = 1.5, main = "Chr21 PCA: by Sex",
     xlab = paste0("PC1: ", percent_var[1], "% variance"),
     ylab = paste0("PC2: ", percent_var[2], "% variance"))
legend("topright", legend = c("Female", "Male"),
       col = c("purple", "green"), pch = 19)

# PC1 vs PC2 by age (color gradient)
age_colors <- colorRampPalette(c("lightblue", "darkblue"))(100)
age_bins <- cut(pca_data$age, breaks = 100)
plot(pca_data$PC1, pca_data$PC2,
     col = age_colors[as.numeric(age_bins)],
     pch = 19, cex = 1.5, main = "Chr21 PCA: by Age",
     xlab = paste0("PC1: ", percent_var[1], "% variance"),
     ylab = paste0("PC2: ", percent_var[2], "% variance"))

# PC1 vs PC2 by visit
visit_colors <- c("Visit 1" = "orange", "Visit 2" = "brown",
                  "Visit 3" = "pink", "Visit 4" = "yellow")
plot(pca_data$PC1, pca_data$PC2,
     col = visit_colors[as.character(pca_data$visit)],
     pch = 19, cex = 1.5, main = "Chr21 PCA: by Visit",
     xlab = paste0("PC1: ", percent_var[1], "% variance"),
     ylab = paste0("PC2: ", percent_var[2], "% variance"))
legend("topright", legend = names(visit_colors),
       col = visit_colors, pch = 19, cex = 0.8)

# Box plots
boxplot(PC1 ~ karyotype, data = pca_data,
        col = c("blue", "red"),
        main = "PC1 by Karyotype",
        ylab = "PC1")

boxplot(PC2 ~ karyotype, data = pca_data,
        col = c("blue", "red"),
        main = "PC2 by Karyotype",
        ylab = "PC2")

dev.off()

cat("\nSaved: results/figures/pca_chr21_detailed.pdf\n")

# Summary statistics
cat("\n=== Summary Statistics ===\n\n")
cat("PC1 by Karyotype:\n")
cat(sprintf("  Control: mean=%.2f, sd=%.2f, range=[%.2f, %.2f]\n",
            mean(pca_data$PC1[pca_data$karyotype == "Control"]),
            sd(pca_data$PC1[pca_data$karyotype == "Control"]),
            min(pca_data$PC1[pca_data$karyotype == "Control"]),
            max(pca_data$PC1[pca_data$karyotype == "Control"])))
cat(sprintf("  T21: mean=%.2f, sd=%.2f, range=[%.2f, %.2f]\n",
            mean(pca_data$PC1[pca_data$karyotype == "T21"]),
            sd(pca_data$PC1[pca_data$karyotype == "T21"]),
            min(pca_data$PC1[pca_data$karyotype == "T21"]),
            max(pca_data$PC1[pca_data$karyotype == "T21"])))

cat("\nPC2 by Karyotype:\n")
cat(sprintf("  Control: mean=%.2f, sd=%.2f, range=[%.2f, %.2f]\n",
            mean(pca_data$PC2[pca_data$karyotype == "Control"]),
            sd(pca_data$PC2[pca_data$karyotype == "Control"]),
            min(pca_data$PC2[pca_data$karyotype == "Control"]),
            max(pca_data$PC2[pca_data$karyotype == "Control"])))
cat(sprintf("  T21: mean=%.2f, sd=%.2f, range=[%.2f, %.2f]\n",
            mean(pca_data$PC2[pca_data$karyotype == "T21"]),
            sd(pca_data$PC2[pca_data$karyotype == "T21"]),
            min(pca_data$PC2[pca_data$karyotype == "T21"]),
            max(pca_data$PC2[pca_data$karyotype == "T21"])))

cat("\nPC2 by Sex:\n")
cat(sprintf("  Female: mean=%.2f, sd=%.2f\n",
            mean(pca_data$PC2[pca_data$sex == "Female"]),
            sd(pca_data$PC2[pca_data$sex == "Female"])))
cat(sprintf("  Male: mean=%.2f, sd=%.2f\n",
            mean(pca_data$PC2[pca_data$sex == "Male"]),
            sd(pca_data$PC2[pca_data$sex == "Male"])))
