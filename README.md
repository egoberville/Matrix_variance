---
title: 'Pairwise common variance analysis using modified RV coefficient'
author: "Eric GOBERVILLE - eric.goberville@sorbonne-universite.fr"
output:
  html_document:
    includes:
      in_header: header.html
    theme: cosmo
    highlight: tango
    toc: true
    toc_float: 
      collapsed: false
      smooth_scroll: false
    toc_depth: 5
    code_folding: hide
    css: custom.css
date: "`r format(Sys.time(), '%d %B, %Y')`"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
```

#
##### <i class="fas fa-book"></i> Introduction
<div style="border-bottom: 1px solid #ddd; margin-bottom: 10px;"></div>

This document provides an analysis of pairwise common variance between multiple matrices using a modified RV coefficient. The analysis includes cleaning matrices, calculating covariances, and visualizing the results with a heatmap and Venn diagram.

#
##### <i class="fas fa-book"></i> RV coefficient: a short history
<div style="border-bottom: 1px solid #ddd; margin-bottom: 10px;"></div>

The RV coefficient, introduced by Escoufier (1973) and further elaborated by Robert and Escoufier (1976), measures the similarity between squared symmetric matrices, specifically positive semi-definite matrices. It serves as a theoretical tool for analyzing multivariate techniques. Escoufier (1973) defined the RV coefficient as a similarity measure for positive semi-definite matrices, and Robert and Escoufier (1976) highlighted its significant mathematical properties, noting that many multivariate analysis techniques effectively maximize this coefficient under suitable constraints. To compare rectangular matrices using the RV coefficient, the initial step is to transform them into square matrices. The RV coefficient ranges between 0 and +1, as it is used with positive semi-definite matrices.


#
##### <i class="fas fa-cogs"></i> Load Libraries
<div style="border-bottom: 1px solid #ddd; margin-bottom: 10px;"></div>

```{r warning = FALSE, message=FALSE}
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggforce)
library(kableExtra)
library(ggthemes)
library(viridis)
theme_set(theme_minimal())
```


#
##### <i class="fas fa-database"></i> Generate Example Data
<div style="border-bottom: 1px solid #ddd; margin-bottom: 10px;"></div>

To make this analysis reproducible, we will generate example matrices with random data and introduce some NA values.


```{r}
set.seed(123)

# Function to generate a matrix with random NAs
generate_matrix_with_na <- function(nrow, ncol, na_prop) {
  mat <- matrix(rnorm(nrow * ncol), nrow, ncol)
  na_count <- floor(nrow * ncol * na_prop)
  na_indices <- sample(nrow * ncol, na_count)
  mat[na_indices] <- NA
  return(mat)
}

# Example matrices with different proportions of NAs
X <- generate_matrix_with_na(400, 100, 0.5)  # 50% NAs
Y <- generate_matrix_with_na(400, 100, 0.35) # 35% NAs
Z <- generate_matrix_with_na(400, 100, 0.2)  # 20% NAs
W <- generate_matrix_with_na(400, 100, 0.05) # 5% NAs
```


#
##### <i class="fas fa-cogs"></i> Function Definitions
<div style="border-bottom: 1px solid #ddd; margin-bottom: 10px;"></div>

Here, we define functions to clean matrices, calculate the modified RV coefficient, and compare multiple matrices.

*   ***function calc_modified_rv***

The `calc_modified_rv` function is designed to calculate the modified RV coefficient. The function cleans input matrices, calculates their covariances, and computes the modified RV coefficient. This coefficient measures the similarity between two matrices while handling missing data. The function includes options for debugging to provide detailed information about the steps and intermediate results.


```{r}
# Function to calculate the modified RV coefficient
calc_modified_rv <- function(X, Y, debug = FALSE, na_threshold = 0.005) {
  trace <- function(A) sum(diag(A))
  
  clean_matrix <- function(M, threshold) {
    row_keep <- rowSums(!is.na(M)) / ncol(M) > threshold
    col_keep <- colSums(!is.na(M)) / nrow(M) > threshold
    M[row_keep, col_keep, drop = FALSE]
  }
  
  X <- clean_matrix(X, na_threshold)
  Y <- clean_matrix(Y, na_threshold)
  
  if (debug) {
    cat("Dimensions after cleaning: X", dim(X), ", Y", dim(Y), "\n")
    cat("NA count: X", sum(is.na(X)), ", Y", sum(is.na(Y)), "\n")
  }
  
  if (nrow(X) != nrow(Y) || nrow(X) == 0 || ncol(X) == 0 || ncol(Y) == 0) {
    warning("Matrices have incompatible dimensions after cleaning.")
    return(NA)
  }
  
  cov_X <- cov(X, use = "pairwise.complete.obs")
  cov_Y <- cov(Y, use = "pairwise.complete.obs")
  cov_XY <- cov(X, Y, use = "pairwise.complete.obs")
  
  if (debug) {
    cat("NA count in covariance matrices: cov_X", sum(is.na(cov_X)), 
        ", cov_Y", sum(is.na(cov_Y)), ", cov_XY", sum(is.na(cov_XY)), "\n")
  }
  
  complete_indices <- complete.cases(cov_X) & complete.cases(cov_Y) & complete.cases(cov_XY)
  cov_X <- cov_X[complete_indices, complete_indices]
  cov_Y <- cov_Y[complete_indices, complete_indices]
  cov_XY <- cov_XY[complete_indices, complete_indices]
  
  numerator <- trace(cov_XY %*% t(cov_XY))
  denominator <- sqrt(trace(cov_X %*% cov_X) * trace(cov_Y %*% cov_Y))
  
  if (debug) {
    cat("Numerator:", numerator, "Denominator:", denominator, "\n")
  }
  
  if (is.na(denominator) || denominator == 0) {
    warning("Unable to calculate RV coefficient: denominator is NA or zero after data cleaning.")
    return(NA)
  }
  
  (numerator / denominator) * 100
}
```

* ***Parameters***

    + `X`       &nbsp; The first input matrix
    + `Y`       &nbsp; The second input matrix
    + `debug`   &nbsp; A boolean flag for printing debug information. Default is FALSE. If debug is TRUE, prints the dimensions of the cleaned matrices and the count of NA values in each. Prints also the count of NA values in the covariance matrices and the values of the numerator and denominator
    + `na_threshold`   &nbsp; A threshold for retaining rows and columns based on the proportion of non-NA values. Default is 0.005
    + `trace`   &nbsp;Computes the trace of a matrix, which is the sum of the diagonal elements
    + `clean_matrix`   &nbsp; Cleans a matrix by removing rows and columns where the proportion of non-NA values is below the given threshold. The cleaning is applied to both matrices X and Y
    + `dimension Check`   &nbsp; Checks if the cleaned matrices have the same number of rows and non-zero dimensions. If not, a warning is issued and NA is returned
    + `covariance Calculation`   &nbsp; Calculates the covariance matrices for X, Y, and the cross-covariance between X and Y using pairwise complete observations
    + `covariance Matrix Cleaning`   &nbsp; Filters the covariance matrices to keep only the rows and columns with complete cases across all three matrices
    + `RV Coefficient Calculation`   &nbsp; Computes the numerator of the RV coefficient as the trace of the product of cov_XY and its transpose. Computes the denominator as the square root of the product of the traces of cov_X squared and cov_Y squared
    + `denominator Check`   &nbsp; Checks if the denominator is NA or zero, which would make the RV coefficient calculation invalid. If so, a warning is issued and NA is returned
    + `result`   &nbsp; Returns the modified RV coefficient as a percentage, which represents the similarity between the matrices X and Y
    
#

*   ***function compare_multiple_matrices***

The `compare_multiple_matrices` function compares several matrices by calculating the modified RV coefficient for each pair. It uses a nested loop to iterate over all possible pairs of matrices, computes the coefficient using the `calc_modified_rv` function, and stores the results in a symmetric matrix. The final matrix provides a comprehensive view of the common variance between all pairs of input matrices.


```{r}
# Function to compare multiple matrices
compare_multiple_matrices <- function(...) {
  matrices <- list(...)
  n <- length(matrices)
  
  result <- matrix(NA, nrow = n, ncol = n)
  rownames(result) <- colnames(result) <- paste("Matrix", 1:n)
  
  for (i in 1:n) {
    for (j in i:n) {
      common_var <- calc_modified_rv(matrices[[i]], matrices[[j]])
      result[i, j] <- result[j, i] <- common_var
    }
  }
  
  return(result)
}
```


#
##### <i class="fas fa-cogs"></i> Calculate common variance
<div style="border-bottom: 1px solid #ddd; margin-bottom: 10px;"></div>

We calculate the pairwise common variance between the matrices to visualize the results using a heatmap and Venn diagram.

```{r}
# Calculate pairwise common variance
result <- compare_multiple_matrices(X, Y, Z, W)
kable(round(result, 2)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))
```


#
##### <i class="fas fa-project-diagram"></i> Heatmap visualization
<div style="border-bottom: 1px solid #ddd; margin-bottom: 10px;"></div>

The `create_heatmap` function is designed to create a heatmap based on a matrix of pairwise common variance values. This visualization helps to quickly identify the degree of common variance between different pairs of matrices.


```{r}
# Function to create a heatmap
create_heatmap <- function(result_matrix) {
  melted_matrix <- melt(result_matrix)

  ggplot(melted_matrix, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +  # Remove the trace color to make it invisible
    scale_fill_viridis_c(option = "magma", 
                         name = "Common variance (%)", 
                         na.value = "white",  # Set NA values to white
                         limits = c(0, 100)) +  # Cap the color scale at 100%
    theme_minimal() +
    labs(title = "Pairwise common variance between matrices",
         x = "", y = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Create heatmap
heatmap <- create_heatmap(result)
print(heatmap)

```


#
##### <i class="fas fa-project-diagram"></i> Venn diagram visualization
<div style="border-bottom: 1px solid #ddd; margin-bottom: 10px;"></div>

We then visualize the common variance using a Venn diagram. This code generates a Venn diagram based on the common variance among matrices, taking W as the reference matrix to represent the common variance with the other matrices considered (X, Y, and Z). It involves extracting common variance values, normalizing them, creating a data frame for circle properties, and using ggplot2 to plot the circles and labels. The final plot visually represents the pairwise common variance.

```{r venn-diagram}
# Calculate sizes and positions based on common variance with W
w_common_var <- result[4, 1:3]
max_size <- max(w_common_var)
sizes <- w_common_var / max_size

# Create a data frame for the circles
circle_data <- data.frame(
  x0 = c(0, cos(0) * (1 - sizes[1]/2), cos(2*pi/3) * (1 - sizes[2]/2), cos(4*pi/3) * (1 - sizes[3]/2)),
  y0 = c(0, sin(0) * (1 - sizes[1]/2), sin(2*pi/3) * (1 - sizes[2]/2), sin(4*pi/3) * (1 - sizes[3]/2)),
  r = c(0.5, sizes[1]/2, sizes[2]/2, sizes[3]/2),
  name = c("W", "X", "Y", "Z")
)

# Create the plot with improved readability
venn_plot <- ggplot() +
  geom_circle(data = circle_data, aes(x0 = x0, y0 = y0, r = r, fill = name), alpha = 0.4) +
  geom_text(data = circle_data, aes(x = x0, y = y0, label = name), fontface = "bold", size = 6) +
  scale_fill_brewer(palette = "Set2") +
  theme_void() +
  theme(legend.position = "none") +
  coord_fixed(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2))

# Add common variance labels with slight offset
for(i in 2:4) {
  common_var <- round(result[4,i-1], 1)
  angle <- (i-2) * 2*pi/3
  offset <- 0.1
  x <- cos(angle) * (0.5 + circle_data$r[i]/2 + offset)
  y <- sin(angle) * (0.5 + circle_data$r[i]/2 + offset)
  venn_plot <- venn_plot +
    annotate("text", x = x, y = y, 
             label = paste0(common_var, "%"),
             size = 5)
}

# Print the plot
print(venn_plot)

```

* ***Parameters***

    + `result`       &nbsp; is assumed to be a matrix containing the pairwise common variance (modified RV coefficient) between multiple matrices
    + `max_size`       &nbsp; calculates the maximum common variance value among the extracted values
    + `sizes`   &nbsp; normalizes the common variance values by dividing each value by the maximum value, resulting in values between 0 and 1
#    
    A data frame `circle_data` is created to store the properties of the circles to be plotted
#    
    + `x0` and `y0` represent the center coordinates of the circles. The positions are calculated based on the angles and the normalized sizes
    + The `cos` and `sin` functions are used to place the circles at equal angles around the origin
    + `r` represents the radii of the circles, set as half of the normalized sizes, with the first circle (representing W) having a fixed radius of 0.5
    + `name` contains the names of the matrices
    

#
##### <i class="fas fa-lightbulb"></i> Conclusion
<div style="border-bottom: 1px solid #ddd; margin-bottom: 10px;"></div>

In this analysis, we have calculated the pairwise common variance between multiple matrices and visualized the results using a heatmap and Venn diagram. This approach helps in understanding the shared variance and similarities between different data sets.




