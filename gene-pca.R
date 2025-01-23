library(ggbiplot)
library(ggplot2)
library(dplyr)

gene_data <- read.csv("C:\\Users\\gatit\\OneDrive\\Desktop\\BIOAID_tpm_PC0.001_log2_genesymbol_dedup.csv")

transposed_data <- as.data.frame( t(gene_data) ) # Transposes gene data

new_headers <- transposed_data[1,] # The column headings (gene names) are stored
transposed_data <- transposed_data[-1,] # Removes the first row (gene names) from this dataframe
colnames(transposed_data) <- new_headers # The column headings are replaced with new_headers
transposed_data[] <- lapply(transposed_data, as.numeric) # Treats str characters as numeric

PCA_ready_data <- transposed_data %>%
  mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

zero_varience_cols <- apply(PCA_ready_data,2, var) == 0

if (any(zero_varience_cols)) {
  cat("Zero-variance columns detected and removed:\n")
  print(colnames(PCA_ready_data)[zero_varience_cols])
}



# Perform PCA
pca_result <- prcomp(normalized_data, center = TRUE, scale. = TRUE)

# Summary of PCA
print(summary(pca_result))  # Variance explained by each principal component

# Variance explained by each principal component
explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)
print("Explained variance by each principal component:")
print(explained_variance)


# Visualize PCA using ggbiplot
pca_plot <- ggbiplot(pca_result, obs.scale = 1, var.scale = 1,
                     groups = gene_data$X,  # Grouping by gene IDs
                     ellipse = TRUE, circle = TRUE) +
  labs(title = "PCA of Gene Expression Data") +
  theme_minimal()

print(pca_plot)

# Scree plot to visualize variance explained
scree_plot <- ggplot(data.frame(PC = 1:length(explained_variance), Variance = explained_variance),
                     aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(title = "Scree Plot", x = "Principal Components", y = "Variance Explained") +
  theme_minimal()

print(scree_plot)
  
  
  
  
  