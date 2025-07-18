# TreeOrderTests

**TreeOrderTests** is an R package for testing the equality of means against **tree ordered alternatives** in a one-way ANOVA model.

---

## 📊 Overview

The package provides three tests:
- **TreeLRT**: Likelihood Ratio Test  
- **TreeMaxD**: Maximum Difference Based Test  
- **TreeMinD**: Minimum Difference Based Test

These tests address the null hypothesis:

\[
H_0: \mu_0 = \mu_1 = \cdots = \mu_k
\quad \text{vs} \quad
H_1: \mu_0 \leq \mu_i, \quad i = 1,\ldots,k, \quad \text{with at least one strict inequality.}
\]

Here, \(\mu_0\) is the mean of the control group and \(\mu_i\) are the means of the treatment groups.

---

## 📦 Installation

install.packages("devtools")
library(devtools)
install_github("SubhaHalder-spec/TreeOrderTests", force = TRUE)
library(TreeOrderTests)

# Example data: control + 2 treatments
control <- rnorm(30, mean = 10, sd = 2)
treatment1 <- rnorm(30, mean = 11, sd = 2)
treatment2 <- rnorm(30, mean = 12, sd = 2)

data_list <- list(control, treatment1, treatment2)

# Run TreeLRT
TreeLRT(data_list, significance_level = 0.05)

# Run TreeMaxD
TreeMaxD(data_list, significance_level = 0.05)

# Run TreeMinD
TreeMinD(data_list, significance_level = 0.05)
