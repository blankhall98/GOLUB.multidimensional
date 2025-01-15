library(ggplot2)
library(dplyr)
library(tidyr)
library(corrplot)
library(dplyr)
library(reshape2)

golub_data <- read.csv("data/golub.csv")

summary(golub_data)
dim(golub_data) # 72 7135
sum(is.na(golub_data))

# Identify columns with missing data
missing_data <- colSums(is.na(golub_data))
columns_with_missing <- missing_data[missing_data > 0]
cat("Columns with Missing Data and their Counts:\n")
print(columns_with_missing) # only missing values in GENDER


# Count columns by type
num_character <- sum(sapply(golub_data, is.character))  # Count character columns
num_numeric <- sum(sapply(golub_data, is.numeric))      # Count numeric columns

# Display the results
cat("Number of Character Columns:", num_character, "\n")
cat("Number of Numeric Columns:", num_numeric, "\n")

character_vars <- names(golub_data)[sapply(golub_data, is.character)]  # Identify character variables
cat("Character Variables:\n", character_vars, "\n\n")

# Check distribution of factors
table(golub_data$BM.PB)
table(golub_data$Gender)
table(golub_data$cancer)
table(golub_data$Source)
table(golub_data$tissue.mf)

# Combine categories 'allB' and 'allT' into 'all'
golub_data$cancer <- as.factor(ifelse(golub_data$cancer %in% c("allB", "allT"), "all", "aml"))


# Visualize the distribution of categorical variables
ggplot(golub_data, aes(x=BM.PB, fill=BM.PB)) +
  geom_bar() +
  labs(title = "Distribution of Samples by BM.PB", x="BM.PB", y="Count") +
  theme_minimal()

ggplot(golub_data, aes(x = cancer, fill = cancer)) +
  geom_bar() +
  labs(title = "Distribution of Samples by Cancer Type", x = "Cancer Type", y = "Count") +
  theme_minimal()

ggplot(golub_data, aes(x = Gender, fill = Gender)) +
  geom_bar() +
  labs(title = "Distribution of Samples by Gender", x = "Gender", y = "Count") +
  theme_minimal()

ggplot(golub_data, aes(x = Source, fill = Source)) +
  geom_bar() +
  labs(title = "Distribution of Samples by Source (Hospital)", x = "Hospital", y = "Count") +
  theme_minimal()

ggplot(golub_data, aes(x = tissue.mf, fill = tissue.mf)) +
  geom_bar() +
  labs(title = "Distribution of Samples by tissue", x = "Combination of source and gender", y = "Count") +
  theme_minimal()

# Correlation among gene probes
# Randomly select a subset of gene probes for correlation analysis (too many to visualize all at once)
set.seed(123)
sampled_genes <- sample(7:ncol(golub_data), 50)  # Randomly select 50 genes
gene_cor_matrix <- cor(golub_data[, sampled_genes])
corrplot(gene_cor_matrix, method = "circle", type = "upper", tl.cex = 0.7)

selected_genes <- sampled_genes[1:3]
for (gene in selected_genes) {
  gene_name <- colnames(golub_data)[gene]
  plot <- ggplot(golub_data, aes(x = cancer, y = golub_data[[gene]], fill = cancer)) +
    geom_boxplot() +
    labs(title = paste("Expression of", gene_name, "by Cancer Type"), x = "Cancer Type", y = "Expression Level") +
    theme_minimal()
  print(plot)
}

# QUANTITATIVE ANALYSIS
# Initialize a results data frame for t-test results
t_test_results <- data.frame(
  Gene = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# Perform t-test for each gene probe
for (gene in 7:ncol(golub_data)) {
  t_test <- t.test(
    golub_data[[gene]] ~ golub_data$cancer,
    data = golub_data
  )
  t_test_results <- rbind(t_test_results, data.frame(
    Gene = colnames(golub_data)[gene],
    P_Value = t_test$p.value
  ))
}

# Adjust p-values for multiple comparisons using Bonferroni correction
t_test_results$Adjusted_P_Value <- p.adjust(t_test_results$P_Value, method = "bonferroni")

# Filter significant genes (e.g., Adjusted P-Value < 0.05)
significant_genes <- t_test_results[t_test_results$Adjusted_P_Value < 0.05, ]
print(significant_genes)

nrow(significant_genes)


# HEATMAP
# Extract significant gene data
significant_gene_data <- golub_data[, colnames(golub_data) %in% significant_genes$Gene, drop = FALSE]
rownames(significant_gene_data) <- golub_data$SampleID  # Ensure row names correspond to samples

# Add cancer type as a column for annotations
significant_gene_data <- cbind(Cancer = golub_data$cancer, significant_gene_data)

# Melt the data for ggplot2
melted_data <- melt(significant_gene_data, id.vars = "Cancer")
colnames(melted_data) <- c("Cancer", "Gene", "Expression")

# Order genes by significance (optional)
gene_order <- significant_genes$Gene[order(significant_genes$Adjusted_P_Value)]
melted_data$Gene <- factor(melted_data$Gene, levels = gene_order)

# Split genes into smaller groups
gene_groups <- split(levels(melted_data$Gene), ceiling(seq_along(levels(melted_data$Gene)) / 10))  # Groups of ~30 genes

# 1. Using all observations
# Generate heatmaps for each group
for (i in seq_along(gene_groups)) {
  # Select genes in the current group
  genes_to_plot <- gene_groups[[i]]
  
  # Subset data for the selected genes
  subset_data <- melted_data[melted_data$Gene %in% genes_to_plot, ]
  
  # Add a placeholder for samples (e.g., row number as unique identifier)
  subset_data$Sample <- factor(1:nrow(subset_data))
  
  # Create the heatmap
  p <- ggplot(subset_data, aes(x = Gene, y = Sample, fill = Expression)) +
    geom_tile() +
    facet_wrap(~Cancer, scales = "free_y") +  # Facet by cancer type
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(
      title = paste("Heatmap of Significant Genes (Group", i, ")"),
      x = "Gene",
      y = "Sample"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 6),
      strip.text = element_text(size = 10, face = "bold")
    )
  
  # Display the heatmap for the current group
  print(p)
}


# 2. SHOWING 15 OBSERVATIONS OF EACH CANCER TYPE
# Generate heatmaps for each group
for (i in seq_along(gene_groups)) {
  # Select genes in the current group
  genes_to_plot <- gene_groups[[i]]
  
  # Subset data for the selected genes
  subset_data <- melted_data[melted_data$Gene %in% genes_to_plot, ]
  
  # Select 15 random samples for each cancer type
  subset_data_filtered <- subset_data %>%
    group_by(Cancer) %>%
    sample_n(15)  # Select 15 random samples per cancer type
  
  # Add a placeholder for samples (e.g., row number as unique identifier)
  subset_data_filtered$Sample <- factor(1:nrow(subset_data_filtered))
  
  # Create the heatmap
  p <- ggplot(subset_data_filtered, aes(x = Gene, y = Sample, fill = Expression)) +
    geom_tile() +
    facet_wrap(~Cancer, scales = "free_y") +  # Facet by cancer type
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(
      title = paste("Heatmap of Significant Genes (Group", i, ")"),
      x = "Gene",
      y = "Sample"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 6),
      strip.text = element_text(size = 10, face = "bold")
    )
  
  # Display the heatmap
  print(p)
}
