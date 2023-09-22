file_path <- file.choose()
data2 <- read.csv(file_path, header = TRUE)

# Remove null data
data2$X <- NULL

# Reshape the dataset
data2 <- data2[, -1]
data2$diagnosis <- factor(ifelse(data2$diagnosis == "B", "Benign", "Malignant"))

# Display the structure, head, and summary of the dataset
cat("Structure of the dataset:\n")
str(data2)

cat("\nHead of the dataset:\n")
head(data2)

cat("\nSummary of the dataset:\n")
summary(data2)

# Load libraries
library(PerformanceAnalytics)
library(psych)
library(ggplot2)
library(GGally)

# Function to create correlation chart
createCorrelationChart <- function(data, title) {
  chart.Correlation(
    data[, c(2:11)],
    histogram = TRUE,
    col = "grey10",
    pch = 1,
    main = title
  )
}

# Function to create scatterplot matrix
createScatterplotMatrix <- function(data, title) {
  pairs.panels(
    data[, c(12:21)],
    method = "pearson",
    hist.col = "#1fbbfa",
    density = TRUE,
    ellipses = TRUE,
    show.points = TRUE,
    pch = 1,
    lm = TRUE,
    cex.cor = 1,
    smoother = FALSE,
    stars = TRUE,
    main = title
  )
}

# Function to create ggpairs plot
createGGPairsPlot <- function(data, title) {
  ggpairs(data[, c(22:31)]) +
    theme_bw() +
    labs(title = title) +
    theme(
      plot.title = element_text(
        face = 'bold',
        color = 'black',
        hjust = 0.5,
        size = 13
      )
    )
}

# Data
data2 <- data2

# Create and display correlation chart
createCorrelationChart(data2, "Cancer Mean")

# Create and display scatterplot matrix
createScatterplotMatrix(data2, "Cancer SE")

# Create and display ggpairs plot
createGGPairsPlot(data2, "Cancer Worst")
##Variables
# Load libraries
library(ggplot2)
library(GGally)

# Function for creating correlation plots
create_corr_plot <- function(data2, vars, title) {
  chart.Correlation(data[, vars], histogram = TRUE, col = "grey10", pch = 1, main = title)
}

# Function for creating ggpairs plots
create_ggpairs_plot <- function(data, vars, title) {
  ggpairs(data2[, vars], aes(color = diagnosis, alpha = 0.75), lower = list(continuous = "smooth")) +
    theme_bw() +
    labs(title = title) +
    theme(plot.title = element_text(face = 'bold', color = 'black', hjust = 0.5, size = 12))
}

# Function for creating ggcorr plots
create_ggcorr_plot <- function(data, vars, title) {
  ggcorr(data2[, vars], name = "corr", label = TRUE) +
    theme(legend.position = "none") +
    labs(title = title) +
    theme(plot.title = element_text(face = 'bold', color = 'black', hjust = 0.5, size = 12))
}

# Plots
data2 <- data2
create_corr_plot(data2, c(2:11), "Cancer Mean")
create_ggpairs_plot(data2, c(2:11), "Cancer Mean")
create_ggcorr_plot(data2, c(2:11), "Cancer Mean")
create_corr_plot(data2, c(12:21), "Cancer SE")
create_ggpairs_plot(data2, c(12:21), "Cancer SE")
create_ggcorr_plot(data2, c(12:21), "Cancer SE")
create_corr_plot(data2, c(22:31), "Cancer Worst")
create_ggpairs_plot(data2, c(22:31), "Cancer Worst")
create_ggcorr_plot(data2, c(22:31), "Cancer Worst")
##PCA
# Load necessary libraries
library(factoextra)
library(corrplot)
library(gridExtra)

data2_pca <- transform(data2)
# Function to perform PCA and return summary
all_pca <- prcomp(data2_pca[,-1], cor=TRUE, scale = TRUE)
summary(all_pca)


  # Screeplot
fviz_eig(all_pca, addlabels=TRUE, ylim=c(0,60), geom = c("bar", "line"), barfill = "pink", barcolor="grey",linecolor = "red", ncp=10)+
  labs(title = "Cancer All Variances - PCA",
       x = "Principal Components", y = "% of variances")

  # Get PCA Variables
  all_var <- get_pca_var(all_pca)
  all_var

  # Correlation between variables and PCA
  library("corrplot")
  corrplot(all_var$cos2, is.corr=FALSE)

  # Contributions of variables to PC1 & PC2
  library(gridExtra)
  p1 <- fviz_contrib(all_pca, choice="var", axes=1, fill="pink", color="grey", top=10)
  p2 <- fviz_contrib(all_pca, choice="var", axes=2, fill="skyblue", color="grey", top=10)
  grid.arrange(p1,p2,ncol=2)

# Function to perform k-means clustering and visualize PCA variables by groups
  set.seed(218)
  res.all <- kmeans(all_var$coord, centers = 6, nstart = 25)
  grp <- as.factor(res.all$cluster)

  fviz_pca_var(all_pca, col.var = grp,
               palette = "jco",
               legend.title = "Cluster")

# Split data into train and test sets
nrows <- NROW(data2)
set.seed(218)                           # Fix random seed for reproducibility
index <- sample(1:nrows, 0.7 * nrows)   # Shuffle and divide data

train <- data2[index,]                   # 70% training data
test <- data2[-index,]                  # 30% test data

# Check class distribution in train and test sets
prop.table(table(train$diagnosis))
prop.table(table(test$diagnosis))

##1.KNN
library(caret)
library(class)
library(highcharter)

# Define the range of k values to test
k_values <- 1:30

# Create a dataframe to store cross-validation results
cv_results <- data.frame(K = integer(), Accuracy = numeric())

# Define the training control parameters for cross-validation
ctrl <- trainControl(method = "cv",  # Use k-fold cross-validation
                     number = 5,     # Number of folds (adjust as needed)
                     verboseIter = TRUE)

# Perform cross-validated training and hyperparameter tuning
set.seed(123)  # Set a random seed for reproducibility
for (k_value in k_values) {
  # Create a container to store accuracy for each fold
  fold_accuracies <- numeric()

  # Perform k-fold cross-validation
  folds <- createFolds(train$diagnosis, k = 5)  # Adjust the number of folds as needed
  for (fold in 1:5) {  # Loop through the folds
    train_indices <- unlist(folds[-fold])
    validation_indices <- unlist(folds[fold])

    # Train the KNN model with the current k value
    predict_knn <- knn(train = train[train_indices, -1],
                       test = train[validation_indices, -1],
                       cl = train[train_indices, 1],
                       k = k_value,
                       prob = TRUE)

    # Calculate accuracy for the current fold
    fold_accuracies[fold] <- mean(predict_knn == train[validation_indices, 1])
  }

  # Compute the mean accuracy across all folds for the current k value
  mean_accuracy <- mean(fold_accuracies)

  # Store the k value and its corresponding mean accuracy
  cv_results <- rbind(cv_results, data.frame(K = k_value, Accuracy = mean_accuracy))
}

# Find the k value with the highest mean accuracy
optimal_k <- cv_results$K[which.max(cv_results$Accuracy)]
optimal_accuracy <- max(cv_results$Accuracy)

# Visualize the cross-validation results
hchart(cv_results, 'line', hcaes(K, Accuracy)) %>%
  hc_title(text = "Cross-Validation Accuracy With Varying K (KNN)") %>%
  hc_subtitle(paste("Optimal K value is", optimal_k, "(Accuracy:", optimal_accuracy, ")")) %>%
  hc_add_theme(hc_theme_google()) %>%
  hc_xAxis(title = list(text = "Number of Neighbors (K)")) %>%
  hc_yAxis(title = list(text = "Accuracy"))

# Apply the optimal K to show the best predictive performance in KNN
predict_knn <- knn(train = train[,-1], test = test[,-1], cl = train[,1], k = optimal_k, prob = TRUE)
# Evaluate the model using confusion matrix
cm_knn <- confusionMatrix(predict_knn, test$diagnosis)

# Extract relevant metrics (accuracy, F1-score, recall)
accuracy <- cm_knn$overall["Accuracy"]
f1_score <- cm_knn$byClass["F1"]
recall <- cm_knn$byClass["Recall"]
precision <- cm_knn$byClass["Precision"]

# Print the confusion matrix and metrics
print(cm_knn)
cat("Accuracy:", accuracy, "\n")
cat("F1-score:", f1_score, "\n")
cat("Recall:", recall, "\n")
cat("Precision:", precision, "\n")
##2.SVM
library(e1071)
library(caret)

# Define the range of gamma and cost values to test
gamma <- seq(0, 0.1, 0.005)
cost <- 2^(0:5)
parms <- expand.grid(cost = cost, gamma = gamma)  ## 231

# Create a dataframe to store cross-validation results
cv_results <- data.frame(Gamma = numeric(), Cost = numeric(), Accuracy = numeric())

# Define the training control parameters for cross-validation
ctrl <- trainControl(method = "cv",  # Use k-fold cross-validation
                     number = 5,     # Number of folds (adjust as needed)
                     verboseIter = TRUE)

# Perform cross-validated training and hyperparameter tuning
for (i in 1:NROW(parms)) {
  # Create a container to store accuracy for each fold
  fold_accuracies <- numeric()

  # Perform k-fold cross-validation
  set.seed(123)  # Set a random seed for reproducibility
  folds <- createFolds(train$diagnosis, k = 5)  # Adjust the number of folds as needed
  for (fold in 1:5) {  # Loop through the folds
    train_indices <- unlist(folds[-fold])
    validation_indices <- unlist(folds[fold])

    # Train the SVM model with the current gamma and cost
    learn_svm <- svm(diagnosis ~ ., data = train[train_indices, ], gamma = parms$gamma[i], cost = parms$cost[i])

    # Make predictions on the validation fold
    pre_svm <- predict(learn_svm, train[validation_indices, -1])

    # Calculate accuracy for the current fold
    fold_accuracies[fold] <- sum(pre_svm == train[validation_indices, "diagnosis"]) / length(pre_svm)
  }

  # Compute the mean accuracy across all folds for the current gamma and cost
  mean_accuracy <- mean(fold_accuracies)

  # Store the gamma, cost, and corresponding mean accuracy in cv_results
  cv_results <- rbind(cv_results, data.frame(Gamma = parms$gamma[i], Cost = parms$cost[i], Accuracy = mean_accuracy))
}

# Find the combination of gamma and cost with the highest mean accuracy
optimal_params <- cv_results[which.max(cv_results$Accuracy), ]
# Visualize the cross-validation results
library(ggplot2)

# Define the color scale for better readability
color_scale <- scale_fill_gradient(
  low = "white",
  high = "blue",
  limits = c(0.5, 1)  # Adjust the limits based on your data range
)

# Create the plot
ggplot(cv_results, aes(x = Gamma, y = Cost, fill = Accuracy)) +
  geom_tile() +
  color_scale +
  labs(
    title = "Grid Search for SVM Hyperparameters with Cross-Validation",
    x = "Gamma",
    y = "Cost",
    fill = "Accuracy"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Apply the optimal parameters (gamma, cost) to train the final SVM model
learn_imp_svm <- svm(diagnosis ~ ., data = train, gamma = optimal_params$Gamma, cost = optimal_params$Cost)
pre_imp_svm <- predict(learn_imp_svm, test[, -1])
# Create the confusion matrix
cm_imp_svm <- confusionMatrix(pre_imp_svm, reference = test$diagnosis)

# Extract relevant metrics (accuracy, recall, F1-score)
accuracy <- cm_imp_svm$overall["Accuracy"]
recall <- cm_imp_svm$byClass["Recall"]
f1_score <- cm_imp_svm$byClass["F1"]
precision <- cm_imp_svm$byClass["Precision"]

# Print the confusion matrix and metrics
print(cm_imp_svm)
cat("Accuracy:", accuracy, "\n")
cat("Recall:", recall, "\n")
cat("Precision:", precision, "\n")
cat("F1-score:", f1_score, "\n")
##3.NB
library(e1071)
library(caret)

# Define the range of Laplace values to test
laplace_values <- seq(0.1, 2.0, by = 0.1)  # Adjust the range as needed

# Create a dataframe to store the cross-validation results
cv_results <- data.frame(Laplace = numeric(), Accuracy = numeric())

# Perform k-fold cross-validation for each Laplace value
num_folds <- 5  # Adjust the number of folds as needed
set.seed(123)   # Set a random seed for reproducibility
for (laplace_value in laplace_values) {
  # Create a container to store accuracy for each fold
  fold_accuracies <- numeric()

  # Perform k-fold cross-validation (stratified)
  folds <- createFolds(train$diagnosis, k = num_folds, list = TRUE, returnTrain = FALSE)
  for (fold in 1:num_folds) {
    validation_indices <- unlist(folds[[fold]])
    train_indices <- setdiff(1:nrow(train), validation_indices)

    # Train the naive Bayes model with the current Laplace value
    learn_nb <- naiveBayes(train[train_indices, -1], train[train_indices, "diagnosis"], laplace = laplace_value)

    # Make predictions on the validation fold
    pre_nb <- predict(learn_nb, train[validation_indices, -1])

    # Calculate accuracy for the current fold
    fold_accuracies[fold] <- sum(pre_nb == train[validation_indices, "diagnosis"]) / length(pre_nb)
  }

  # Compute the mean accuracy across all folds for the current Laplace value
  mean_accuracy <- mean(fold_accuracies)

  # Store the Laplace value and its corresponding mean accuracy
  cv_results <- rbind(cv_results, data.frame(Laplace = laplace_value, Accuracy = mean_accuracy))
}

# Find the Laplace value with the highest mean accuracy
optimal_laplace <- cv_results$Laplace[which.max(cv_results$Accuracy)]
optimal_accuracy <- max(cv_results$Accuracy)

# Visualize the cross-validation results
library(ggplot2)
ggplot(cv_results, aes(x = Laplace, y = Accuracy)) +
  geom_line() +
  labs(title = "Cross-Validation Accuracy With Varying Laplace (naiveBayes)",
       x = "Laplace", y = "Accuracy") +
  theme_minimal()

# Train the final naive Bayes model using the optimal Laplace value
learn_imp_nb <- naiveBayes(train[, -1], train$diagnosis, laplace = optimal_laplace)
pre_imp_nb <- predict(learn_imp_nb, test[, -1])

# Create the confusion matrix
cm_imp_nb <- confusionMatrix(pre_imp_nb, reference = test$diagnosis)

# Extract relevant metrics (F1-score, precision, recall, accuracy)
f1_score <- cm_imp_nb$byClass["F1"]
precision <- cm_imp_nb$byClass["Precision"]
recall <- cm_imp_nb$byClass["Recall"]
accuracy <- cm_imp_nb$overall["Accuracy"]

# Print the confusion matrix and metrics
print(cm_imp_nb)
cat("F1-score:", f1_score, "\n")
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("Accuracy:", accuracy, "\n")
##4.Rainforest
library(randomForest)
library(caret)

# Define the training control parameters for cross-validation
ctrl <- trainControl(method = "cv",  # Use k-fold cross-validation
                     number = 5,     # Number of folds (adjust as needed)
                     verboseIter = TRUE)

# Create a data frame to store cross-validation results
cv_results <- data.frame(Accuracy = numeric(), Kappa = numeric())

# Define the hyperparameters to tune
tune_grid <- expand.grid(mtry = c(1, 2, 3))  # Adjust the values of 'mtry' as needed

# Perform cross-validated training and hyperparameter tuning
set.seed(123)  # Set a random seed for reproducibility
learn_rf_cv <- train(diagnosis ~ ., data = train, method = "rf", ntree = 500,
                     tuneGrid = tune_grid, trControl = ctrl, proximity = TRUE, importance = TRUE)

# Get the best model from cross-validation
best_rf <- learn_rf_cv$finalModel

# Make predictions on the test set using the best model
pre_rf_cv <- predict(best_rf, test[, -1])

# Create the confusion matrix
cm_rf_cv <- confusionMatrix(pre_rf_cv, reference = test$diagnosis)

# Extract relevant metrics (recall, precision, accuracy, F1-score)
recall <- cm_rf_cv$byClass["Recall"]
precision <- cm_rf_cv$byClass["Precision"]
accuracy <- cm_rf_cv$overall["Accuracy"]
f1_score <- cm_rf_cv$byClass["F1"]

# Print the confusion matrix and metrics
print(cm_rf_cv)
cat("Recall:", recall, "\n")
cat("Precision:", precision, "\n")
cat("Accuracy:", accuracy, "\n")
cat("F1-score:", f1_score, "\n")

# Plot the error rate vs. the number of trees for the best model
plot(learn_rf_cv, main = "Random Forest (Error Rate vs. Number of Trees)")

# Plot the margin plot for the best model
plot(margin(best_rf, test$diagnosis))

# Plot the variable importance plot for the best model
varImpPlot(best_rf)

# Function to generate and visualize a fourfold plot

generate_fourfold_plot <- function(cm, method_name) {
  col <- c("#ed3b3b", "#0099ff")
  fourfoldplot(cm$table, color = col, conf.level = 0, margin = 1, main = paste(method_name, "(", round(cm$overall[1] * 100), "%)", sep = ""))
}

# Function to select the best prediction model based on accuracy
select_best_prediction_model <- function(confusion_matrices) {
  accuracy <- sapply(confusion_matrices, function(cm) cm$overall[1])
  best_model <- names(accuracy[accuracy == max(accuracy)])
  return(best_model)
}

# Generate and visualize fourfold plots
par(mfrow = c(3, 5))
generate_fourfold_plot(cm_imp_nb, "NB")
generate_fourfold_plot(cm_knn, "KNN")
generate_fourfold_plot(cm_imp_svm, "SVM")
generate_fourfold_plot(cm_rf_cv, "RF")

# Create a list of confusion matrices for different methods
confusion_matrices <- list(cm_imp_nb, cm_knn, cm_imp_svm, cm_rf)

#Select a best prediction model according to high accuracy
opt_predict <- c(cm_imp_nb$overall[1],cm_rf$overall[1], cm_knn$overall[1], cm_rf_cv$overall[1], cm_imp_svm$overall[1])
names(opt_predict) <- c("nb","knn","rf","imp_svm")
best_predict_model <- subset(opt_predict, opt_predict==max(opt_predict))
best_predict_model
# 1. Import patient data
# Replace 'your_data.csv' with the actual file path or use file.choose() to interactively select the file
patient <- read.csv(file_path, header = TRUE)
#Malignant
M <- patient[30,]               ## 30th patient
M[,c(1,2)]                      ## Malignant
#belign
B <- patient[25,]               ## 25th patient
B[,c(1,2)]                      ## Benign
M$diagnosis <- NULL
B$diagnosis <- NULL
# for print output
cancer_diagnosis_predict_p <- function(new, method=learn_imp_svm) {
  new_pre <- predict(method, new[,-1])
  new_res <- as.character(new_pre)
  return(paste("Patient ID: ",new[,1],"  =>  Result: ", new_res, sep=""))
}

# for submission output
cancer_diagnosis_predict_s <- function(new, method=learn_imp_svm) {
  new_pre <- predict(method, new[,-1])
  new_res <- as.character(new_pre)
  return(new_res)
}
#Belign
cancer_diagnosis_predict_p(B)
#Malignant
cancer_diagnosis_predict_p(M)
# Exclude rows with specific indices (index) from the patient data
filtered_patient <- patient[-index, ]

# Extract the actual diagnosis values
origin_diagnosis <- filtered_patient$diagnosis

# Remove the diagnosis column to prepare data for prediction
filtered_patient$diagnosis <- NULL

# Make predictions using cancer_diagnosis_predict_p
predicted_results <- cancer_diagnosis_predict_p(filtered_patient)

# Create a data frame to store the results
results_df <- data.frame(
  id = filtered_patient$id,
  predict_diagnosis = ifelse(predicted_results == 'Malignant', 'M', 'B'),
  orgin_diagnosis = origin_diagnosis
)

# Calculate whether the predictions are correct
results_df$correct <- ifelse(results_df$predict_diagnosis == results_df$orgin_diagnosis, "True", "False")
library(knitr)
# Display the first 20 rows using kable
kable(head(results_df, 20))


##For Web Page
# Load necessary libraries and define your prediction function and model here
library(shiny)
library(knitr)
library(DT)  # For interactive data tables

# Define the UI for the Shiny app
ui <- fluidPage(
  titlePanel("Cancer Diagnosis Prediction"),

  sidebarLayout(
    sidebarPanel(
      fileInput("datafile", "Upload Data File (CSV):", accept = c(".csv")),
      actionButton("predictBtn", "Predict"),
      downloadButton("downloadBtn", "Download Predictions")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Prediction Results", DTOutput("predictionTable"))
      )
    )
  )
)

# Define the server logic for the Shiny app
server <- function(input, output) {
  # Define a reactive dataset to store the loaded data
  data <- reactive({
    req(input$datafile)
    read.csv(input$datafile$datapath, header = TRUE, stringsAsFactors = FALSE)
  })

  # Define a reactive dataset to store the prediction results
  predictionResults <- eventReactive(input$predictBtn, {
    req(data())
    tryCatch({
      # Make predictions using your prediction function and model here
      # Replace this section with your actual prediction code
      # Example: r <- my_predict_function(data(), method = my_model)

      # Placeholder code for testing
      t <- data()
      orgin <- t$diagnosis
      t$diagnosis <- NULL
      r <- sample(c('Malignant', 'Benign'), nrow(t), replace = TRUE)

      # Create the submission dataframe
      sub <- data.frame(
        id = t$id,
        predict_diagnosis = ifelse(r == 'Malignant', 'M', 'B'),
        orgin_diagnosis = orgin
      )
      sub$correct <- ifelse(sub$predict_diagnosis == sub$orgin_diagnosis, "True", "False")
      return(sub)
    }, error = function(e) {
      # Handle errors gracefully and display an error message
      errorMessage <- paste("Error:", e$message)
      showNotification(errorMessage, type = "error")
      return(data.frame(id = NA, predict_diagnosis = NA, orgin_diagnosis = NA, correct = NA))
    })
  })

  # Render the prediction results as an interactive table
  output$predictionTable <- renderDT({
    predictionResults()
  })

  # Define a download handler to download the prediction results
  output$downloadBtn <- downloadHandler(
    filename = function() {
      paste("prediction_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(predictionResults(), file)
    }
  )
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
