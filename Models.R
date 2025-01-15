#STEP 1: Load Necessary Libraries

# Install if not already installed
install.packages(c("tidyverse", "caret", "e1071", "factoextra"))

# Load libraries
library(tidyverse)  # Data manipulation and visualization
library(caret)      # Classification and modeling
library(e1071)      # SVM model
library(factoextra) # Dimensionality reduction visualization

#STEP 2: Load the Data
data <- read.csv("data/golub.csv", row.names = 1)

#STEP 3: Preprocessing
#Convert categorical variables into numerical form using encoding
data$Gender <- as.numeric(factor(data$Gender))
data$Source <- as.numeric(factor(data$Source))
# Separate gene expressions and additional variables
# Assuming "cancer" is the target classification column
gene_data <- data %>% select(-BM.PB, -Gender, -Source, -tissue.mf, -cancer)
#additional_features <- data %>% select(BM.PB, Gender, Source, tissue.mf)
target <- data$cancer

#STEP 4: Perfrom PCA on Gene expression Data
# Perform PCA on gene expression data
pca_result <- prcomp(gene_data, center = FALSE, scale. = FALSE)

# Summary of PCA to check explained variance
summary(pca_result)

# Scree plot to visualize explained variance
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 100))

# Retain enough PCs to explain ~90% variance
pca_data <- as.data.frame(pca_result$x[, 1:10])  # Adjust the number of PCs as needed

# Add additional features and target variable
final_data <- cbind(pca_data, cancer = as.factor(target))
final_data <- cbind(pca_data, additional_features, cancer = as.factor(target))

#STEP 5: Fit classification Model
# Fit an SVM model using all data
svm_model <- svm(cancer ~ ., data = final_data, kernel = "linear")

# Model summary
summary(svm_model)

#STEP 6: Evaluate the Model
# Predict on the same data
predictions <- predict(svm_model, newdata = final_data)

# Confusion matrix and accuracy
# Create confusion matrix
confusion <- table(Predicted = predictions, Actual = final_data$cancer)
print(confusion)

# Calculate accuracy
accuracy <- sum(diag(confusion)) / sum(confusion)
print(paste("Accuracy:", round(accuracy * 100, 2), "%"))

# PCA biplot with class labels
fviz_pca_ind(pca_result,
             geom.ind = "point",
             col.ind = target,  # Color by class
             palette = c("#00AFBB", "#FC4E07"),
             addEllipses = TRUE,
             legend.title = "Cancer Type")
