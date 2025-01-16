#STEP 1: LOAD REQUIRED LIBRARIES
# Install necessary packages if not already installed
install.packages(c("tidyverse", "factoextra", "caret", "e1071", "glmnet"))

# Load the libraries
library(tidyverse)    # Data manipulation and visualization
library(factoextra)   # PCA visualization
library(caret)        # Machine learning framework
library(e1071)        # SVM modeling
library(glmnet)       # Lasso and Ridge regression

#STEP 2: LOAD THE DATA
data <- read.csv("data/golub.csv", row.names = 1)

#obtain only gene data
gene_data <- data %>% select(-BM.PB, -Gender, -Source, -tissue.mf, -cancer)
#obtain target vector
target <- data$cancer

#STEP 3: PERFORM PCA FOR DIMENSIONAL REDUCTION

#Principal Component Analysis (PCA) reduces high-dimensional data into fewer dimensions 
#(principal components, PCs) while retaining the maximum variance. 
#This simplifies the data for classification and visualization.

# Perform PCA
pca_result <- prcomp(gene_data, center = TRUE, scale. = TRUE)

# Summary of PCA
summary(pca_result)

# Scree plot to visualize explained variance
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 100))

# The scree plot shows how much variance each principal component explains.
# Choose enough PCs to explain ~90% of the variance.

# Retain the top 10 PCs
pca_data <- as.data.frame(pca_result$x[, 1:10])

# Add the cancer classification to the PCA-transformed data
pca_data$cancer <- as.factor(data$cancer)

# STEP 4: CLASSIFICATION USING SVM

#Why SVM?
#Support Vector Machines (SVM) are powerful for high-dimensional data, particularly with a linear kernel
#SVM aims to find the hyperplane that best separates classes

# Fit SVM model
svm_model <- svm(cancer ~ ., data = pca_data, kernel = "linear")

# Display model details
summary(svm_model)

#view coefficients
# Extract SVM weights for linear kernel
svm_weights <- t(svm_model$coefs) %*% svm_model$SV

# Display weights
print(svm_weights)

#The weights indicate the contribution of each principal component to the classification.

#Evaluate SVM model
# Predict on the same data
predictions <- predict(svm_model, newdata = pca_data)

# Confusion matrix and accuracy
confusion <- table(Predicted = predictions, Actual = pca_data$cancer)
print(confusion)

accuracy <- sum(diag(confusion)) / sum(confusion)
print(paste("Accuracy:", round(accuracy * 100, 2), "%"))

# Extract coefficients for all support vectors
svm_weights <- t(svm_model$coefs) %*% svm_model$SV

# Convert to a data frame
svm_weights_df <- as.data.frame(svm_weights)

# Assign PCs as feature names
colnames(svm_weights_df) <- paste0("PC", 1:ncol(svm_weights))

# Add a column for class or decision boundary
svm_weights_df$Class <- rownames(svm_weights_df)

# Reshape for visualization
svm_weights_long <- pivot_longer(svm_weights_df, cols = starts_with("PC"), 
                                 names_to = "Principal Component", 
                                 values_to = "Weight")

# Plot weights by PC
library(ggplot2)
ggplot(svm_weights_long, aes(x = `Principal Component`, y = Weight, fill = Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "SVM Weights Across Principal Components",
       x = "Principal Component",
       y = "Weight") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#report section
# Plot decision boundaries on the first two PCs
library(ggplot2)

# Add predictions to the PCA data
pca_data$Prediction <- predictions

# Plot PCA with predictions
ggplot(pca_data, aes(x = PC1, y = PC2, color = Prediction, shape = cancer)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "SVM Decision Boundaries on First Two Principal Components",
       x = "Principal Component 1",
       y = "Principal Component 2",
       color = "Predicted Class",
       shape = "Actual Class")
###
# Summarize SVM weights for each PC
svm_weights_summary <- svm_weights_long %>%
  group_by(`Principal Component`) %>%
  summarize(Mean_Weight = mean(Weight), SD_Weight = sd(Weight)) %>%
  arrange(desc(abs(Mean_Weight)))

# Display the summary table
print(svm_weights_summary)

# Plot mean weights for principal components
ggplot(svm_weights_summary, aes(x = reorder(`Principal Component`, -Mean_Weight), y = Mean_Weight)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Mean SVM Weights Across Principal Components",
       x = "Principal Component",
       y = "Mean Weight") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
###
# Extract support vector indices
support_vector_indices <- svm_model$index

# Display the support vectors
support_vectors <- pca_data[support_vector_indices, ]
print(support_vectors)

# Visualize support vectors on the first two PCs
ggplot(pca_data, aes(x = PC1, y = PC2, color = cancer)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_point(data = support_vectors, aes(x = PC1, y = PC2), color = "red", size = 3, shape = 17) +
  theme_minimal() +
  labs(title = "Support Vectors on First Two Principal Components",
       x = "Principal Component 1",
       y = "Principal Component 2",
       color = "Cancer Type")
###
# Create a confusion matrix as a table
confusion_df <- as.data.frame(as.table(confusion))

# Display the confusion matrix
print(confusion_df)

# Add percentages to the confusion matrix
confusion_df <- confusion_df %>%
  mutate(Percentage = Freq / sum(Freq) * 100)

# Create a heatmap for the confusion matrix
ggplot(confusion_df, aes(x = Predicted, y = Actual, fill = Percentage)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), color = "white") +
  theme_minimal() +
  labs(title = "Confusion Matrix Heatmap",
       x = "Predicted Class",
       y = "Actual Class",
       fill = "Percentage")
###
# Extract PCA loadings (genes by PCs)
loadings <- as.data.frame(pca_result$rotation[, 1:10])  # Loadings for the top 10 PCs
colnames(loadings) <- paste0("PC", 1:10)

# Extract SVM weights for the PCs
svm_weights_vector <- as.numeric(svm_weights[1, ])  # SVM weights for each PC

# Multiply loadings by SVM weights to compute gene contributions
gene_contributions <- loadings %>% 
  mutate(CombinedWeight = rowSums(loadings * svm_weights_vector))

# Add gene names
gene_contributions <- gene_contributions %>%
  rownames_to_column(var = "Gene")

# Sort by the absolute value of the combined weight
gene_contributions <- gene_contributions %>%
  arrange(desc(abs(CombinedWeight)))

# View top contributing genes
head(gene_contributions, 10)

#STEP 5: CLASSIFICATION USING LASSO AND RIDGE REGRESSION

#Lasso Regression adds an L1 penalty, shrinking less important coefficients to zero,
#enabling feature selection.

#Ridge Regression adds an L2 penalty, shrinking all coefficients but retaining them,
#which reduces multicollinearity.

#prepare data for glmnet
# Prepare matrix input for glmnet
x <- as.matrix(gene_data)
y <- as.factor(data$cancer)

#fit lasso model
# Fit multinomial lasso regression
lasso_model <- cv.glmnet(x, y, alpha = 1, family = "multinomial")

# Plot cross-validation curve
plot(lasso_model)

# Coefficients at the best lambda
lasso_coefficients <- coef(lasso_model, s = lasso_model$lambda.min)
lasso_coefficients  # Displays coefficients for each class (allB, allT, aml)

#fit ridge model
# Fit multinomial ridge regression
ridge_model <- cv.glmnet(x, y, alpha = 0, family = "multinomial")

# Plot cross-validation curve
plot(ridge_model)

# Coefficients at the best lambda
ridge_coefficients <- coef(ridge_model, s = ridge_model$lambda.min)
ridge_coefficients  # Displays coefficients for each class (allB, allT, aml)

#STEP 6: EVALUATE LASSO AND RIDGE MODELS

#predict and evaluate lasso
# Predict using lasso
lasso_predictions <- predict(lasso_model, newx = x, s = lasso_model$lambda.min, type = "class")

# Confusion matrix and accuracy
lasso_confusion <- table(Predicted = lasso_predictions, Actual = y)
print(lasso_confusion)

lasso_accuracy <- sum(diag(lasso_confusion)) / sum(lasso_confusion)
print(paste("Lasso Accuracy:", round(lasso_accuracy * 100, 2), "%"))

#predict and evaluate ridge
# Predict using ridge
ridge_predictions <- predict(ridge_model, newx = x, s = ridge_model$lambda.min, type = "class")

# Confusion matrix and accuracy
ridge_confusion <- table(Predicted = ridge_predictions, Actual = y)
print(ridge_confusion)

ridge_accuracy <- sum(diag(ridge_confusion)) / sum(ridge_confusion)
print(paste("Ridge Accuracy:", round(ridge_accuracy * 100, 2), "%"))

