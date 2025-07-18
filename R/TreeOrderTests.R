#' Likelihood Ratio Test for Tree Ordered Alternatives
#'
#' Performs a likelihood ratio test for testing the equality of means against tree ordered alternatives.
#'
#' This test compares the null hypothesis of equal means for all groups to the alternative that the control group mean is less than or equal to the treatment group means under the tree order restriction.
#'
#' @param sample_data A list of numeric vectors. The first element represents the control group, the others represent treatment groups.
#' @param significance_level A numeric value between 0 and 1 specifying the significance level for the test (e.g., 0.05).
#'
#' @return A character string with the critical value, the LRT test statistic, and the test decision.
#'
#' @details The likelihood ratio statistic is computed using constrained maximum likelihood estimates under the null and tree ordered alternative hypotheses. The critical value is estimated by a bootstrap procedure.
#'
#' @importFrom stats quantile rnorm var
#' @export
#'
#' @author Subha Halder

TreeLRT <- function(sample_data, significance_level){
  set.seed(456)
  R_MLE <- function(X, n) {
    X1 <- X[-1]
    n1 <- n[-1]
    sorted_indices <- order(X1)
    X1_sorted <- X1[sorted_indices]
    n1_sorted <- n1[sorted_indices]
    A <- numeric(length(X1_sorted))
    for (j in 2:length(X)) {
      A[j-1] <- (n[1] * X[1] + sum(n1_sorted[1:(j - 1)] * X1_sorted[1:(j - 1)])) /
        (n[1] + sum(n1_sorted[1:(j - 1)]))
    }
    if (all(X1 >= X[1])) {
      new_X <- X
    } else if (A[length(X)-2] >= X1_sorted[length(X)-1]) {
      X <- rep(A[length(X)-1], length(X))
      new_X <- X
    } else {
      comparisons <- logical(length(X1_sorted) - 1)
      comparisons1 <- logical(length(X1_sorted) - 1)
      stored_values <- numeric(0)
      for (k in 1:(length(X1_sorted) - 1)) {
        comparisons[k] <- A[k] < X1_sorted[k + 1]
        if(comparisons1[k] <- A[k] < X1_sorted[k + 1]) {
          for (s in 1:k) {
            stored_values[s] <- X1_sorted[s]
          }
          break
        }
      }
      selected_A_values <- A[comparisons]
      X[1] <- selected_A_values[1]
      for (l in 2:length(X)) {
        if (X[l] %in% stored_values) {
          X[l] <- selected_A_values[1]
        }
      }
      new_X <- X
    }
    return(new_X)
  }

  LRT_H0_new <- function(sample_data_list) {
    means <- sapply(sample_data_list, mean)
    sample_sizes <- sapply(sample_data_list, length)
    S <- unlist(sample_data_list)
    mu1 <- mean(S)
    var1 <- sapply(1:length(sample_data_list), function(i) (sum((sample_data_list[[i]] - means[i])^2)) / sample_sizes[i])
    u1 <- sample_sizes / var1

    repeat {
      new_mu1 <- (sum(u1 * means)) / sum(u1)
      new_var1 <- sapply(1:length(sample_data_list), function(i) (sum((sample_data_list[[i]] - new_mu1)^2)) / sample_sizes[i])
      new_u1 <- sample_sizes / new_var1

      if (max(abs(new_mu1 - mu1)) <= 0.000000000000001) {
        break  # Exit the loop if the difference is less than epsilon
      }

      u1 <- new_u1
      mu1 <- new_mu1
      var1 <- new_var1
    }

    return(var1)
  }

  LRT_H1_new <- function(sample_data_list) {
    n <- sapply(sample_data_list, length)
    mu0 <- sapply(sample_data_list, mean)
    var0 <- sapply(1:length(sample_data_list), function(i) (sum((sample_data_list[[i]] - mu0[i])^2)) / n[i])
    w0 <- n / var0
    repeat {
      new_mu0 <- R_MLE(sapply(sample_data_list, mean), w0)
      new_var0 <- sapply(1:length(sample_data_list), function(i) (sum((sample_data_list[[i]] - new_mu0[i])^2)) / n[i])
      new_w0 <- n / new_var0

      if (max(abs(new_mu0 - mu0)) <= 0.000000000000001) {
        break  # Exit the loop if the difference is less than epsilon
      }

      w0 <- new_w0
      mu0 <- new_mu0
      var0 <- new_var0
    }

    return(var0)
  }

  sample_data <- lapply(sample_data, function(x) x[!is.na(x)])
  num_samples = 100000
  num_datasets <- length(sample_data)
  n <- sapply(sample_data, length)
  lambda_values_star <- numeric(num_samples)
  for (i in 1:num_samples) {
    bootstrap_samples <- lapply(sample_data, function(x) rnorm(n = length(x), mean = 0, sd = sqrt(var(x))))
    V_R_star <- LRT_H1_new(bootstrap_samples) / LRT_H0_new(bootstrap_samples)
    weights <- sapply(1:num_datasets, function(i) V_R_star[i]^(length(bootstrap_samples[[i]]) / 2))
    lambda_values_star[i] <- prod(weights)
  }
  sort_lambda_star <- sort(lambda_values_star)
  quantile_value <- quantile(sort_lambda_star, probs = significance_level)
  V_R <- LRT_H1_new(sample_data) / LRT_H0_new(sample_data)
  weights <- sapply(1:num_datasets, function(i) V_R[i]^(length(sample_data[[i]]) / 2))
  lambda <- prod(weights)
  if (lambda < quantile_value) {
    result <- "Reject null hypothesis"
  } else {
    result <- "Do not reject null hypothesis"
  }
  return(paste("Critical value:", quantile_value, "; LRT Test statistic:", lambda, "; Result:", result))
}


#' Maximum Difference Test for Tree Ordered Alternatives
#'
#' Computes a test statistic based on the maximum standardized difference between the treatment means and the control mean under the tree order restriction.
#'
#' @param sample_data A list of numeric vectors. The first element represents the control group, the others represent treatment groups.
#' @param significance_level A numeric value between 0 and 1 specifying the significance level for the test (e.g., 0.05).
#'
#' @return A character string with the critical value, the Max-D test statistic, and the test decision.
#'
#' @details The test statistic is the maximum of standardized differences between each treatment mean and the control mean. The critical value is estimated by a bootstrap procedure.
#'
#' @importFrom stats quantile rnorm var
#' @export
#'
#' @author Subha Halder

TreeMaxD <- function(sample_data, significance_level){
  set.seed(456)
  sample_data <- lapply(sample_data, function(x) x[!is.na(x)])
  num_samples = 100000
  num_datasets <- length(sample_data)
  n <- sapply(sample_data, length)
  D_star_max <- numeric(num_samples)
  for (i in 1:num_samples) {
    bootstrap_samples <- lapply(sample_data, function(x) rnorm(n = length(x), mean = 0, sd = sqrt(var(x))))
    D_star_max[i] <- max(sapply(2:num_datasets, function(j) {
      (mean(bootstrap_samples[[j]]) - mean(bootstrap_samples[[1]])) /
        sqrt((var(bootstrap_samples[[j]]) / length(bootstrap_samples[[j]])) +
               (var(bootstrap_samples[[1]]) / length(bootstrap_samples[[1]])))
    }))
  }
  sort_D_star_max <- sort(D_star_max)
  quantile_value <- quantile(sort_D_star_max, probs = 1 - significance_level)
  D_Max <- max(sapply(2:num_datasets, function(i) {
    (mean(sample_data[[i]]) - mean(sample_data[[1]])) /
      sqrt(
        (var(sample_data[[i]]) / length(sample_data[[i]])) +
          (var(sample_data[[1]]) / length(sample_data[[1]]))
      )
  }))
  if (D_Max > quantile_value) {
    result <- "Reject null hypothesis"
  } else {
    result <- "Do not reject null hypothesis"
  }
  return(paste("Critical value:", quantile_value, "; D_Max Test statistic:", D_Max, "; Result:", result))
}


#' Minimum Difference Test for Tree Ordered Alternatives
#'
#' Computes a test statistic based on the minimum standardized difference between the treatment means and the control mean under the tree order restriction.
#'
#' @param sample_data A list of numeric vectors. The first element represents the control group, the others represent treatment groups.
#' @param significance_level A numeric value between 0 and 1 specifying the significance level for the test (e.g., 0.05).
#'
#' @return A character string with the critical value, the Min-D test statistic, and the test decision.
#'
#' @details The test statistic is the minimum of standardized differences between each treatment mean and the control mean. The critical value is estimated by a bootstrap procedure.
#'
#' @importFrom stats quantile rnorm var
#' @export
#'
#' @author Subha Halder

TreeMinD <- function(sample_data, significance_level){
  set.seed(456)
  sample_data <- lapply(sample_data, function(x) x[!is.na(x)])
  num_samples = 100000
  num_datasets <- length(sample_data)
  n <- sapply(sample_data, length)
  D_star_min <- numeric(num_samples)
  for (i in 1:num_samples) {
    bootstrap_samples <- lapply(sample_data, function(x) rnorm(n = length(x), mean = 0, sd = sqrt(var(x))))
    D_star_min[i] <- min(sapply(2:num_datasets, function(j) {
      (mean(bootstrap_samples[[j]]) - mean(bootstrap_samples[[1]])) /
        sqrt((var(bootstrap_samples[[j]]) / length(bootstrap_samples[[j]])) +
               (var(bootstrap_samples[[1]]) / length(bootstrap_samples[[1]])))
    }))
  }
  sort_D_star_min <- sort(D_star_min)
  quantile_value <- quantile(sort_D_star_min, probs = 1 - significance_level)
  D_Min <- min(sapply(2:num_datasets, function(i) {
    (mean(sample_data[[i]]) - mean(sample_data[[1]])) /
      sqrt(
        (var(sample_data[[i]]) / length(sample_data[[i]])) +
          (var(sample_data[[1]]) / length(sample_data[[1]]))
      )
  }))
  if (D_Min > quantile_value) {
    result <- "Reject null hypothesis"
  } else {
    result <- "Do not reject null hypothesis"
  }
  return(paste("Critical value:", quantile_value, "; D_Min Test statistic:", D_Min, "; Result:", result))
}


