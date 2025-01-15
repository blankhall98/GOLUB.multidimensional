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
# lasso_coefficients  # Displays coefficients for each class (allB, allT, aml)

# Combine coefficients for all classes
lasso_results <- lapply(names(lasso_coefficients), function(class) {
  coef_matrix <- as.matrix(lasso_coefficients[[class]])
  coef_df <- data.frame(
    Gene_Probe = rownames(coef_matrix),
    Coefficient = coef_matrix[, 1],
    Class = class
  )
  return(coef_df)
})

# Combine results for all classes
lasso_results_df <- bind_rows(lasso_results)

# Filter for non-zero coefficients
non_zero_coefficients <- lasso_results_df %>%
  filter(Coefficient != 0)

# Count unique gene probes
unique_genes <- non_zero_coefficients %>%
  distinct(Gene_Probe) %>%
  nrow()

# Print unique gene probes and their counts
print(paste("Number of unique gene probes:", unique_genes))

# View non-zero coefficients with classes
print(non_zero_coefficients)


#fit ridge model
# Fit multinomial ridge regression
ridge_model <- cv.glmnet(x, y, alpha = 0, family = "multinomial")

# Plot cross-validation curve
plot(ridge_model)

# Coefficients at the best lambda
ridge_coefficients <- coef(ridge_model, s = ridge_model$lambda.min)

# Combine coefficients for all classes into a single data frame
coefficients_df <- do.call(rbind, lapply(1:length(ridge_coefficients), function(i) {
  coef_data <- as.data.frame(as.matrix(ridge_coefficients[[i]]))
  coef_data$Gene <- rownames(coef_data)
  coef_data$Class <- names(ridge_coefficients)[i]
  colnames(coef_data)[1] <- "Coefficient"
  coef_data
}))

# Filter out the intercept and rank by absolute coefficient value
coefficients_df <- coefficients_df[coefficients_df$Gene != "(Intercept)", ]
coefficients_df$AbsCoefficient <- abs(coefficients_df$Coefficient)
top_coefficients <- coefficients_df %>%
  arrange(desc(AbsCoefficient)) %>%
  group_by(Gene) %>%
  slice_max(order_by = AbsCoefficient, n = 1) %>%
  ungroup() %>%
  arrange(desc(AbsCoefficient))

# Display the top 20 most influential coefficients
top_coefficients <- head(top_coefficients, 20)

# Print results
print(top_coefficients)

# Visualize the top coefficients
ggplot(top_coefficients, aes(x = reorder(Gene, -AbsCoefficient), y = Coefficient, fill = Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(
    title = "Top 20 Most Influential Coefficients in Ridge Regression",
    x = "Gene",
    y = "Coefficient"
  ) +
  theme_minimal()

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

