# **intdag**: An R Package For Interventional DAG Estimation and Inference

This repository contains an implementation of the paper [*Inference for a Large Directed Acyclic Graph with Unspecified Interventions*](https://arxiv.org/abs/2110.03805).

## Installation

### R Packages

First, install the package [**glmtlp**](https://github.com/chunlinli/glmtlp) which implements the DC projection algorithm in the paper. 
```r
devtools::install_github("chunlinli/glmtlp")
```
Then install the package [**intdag**](https://github.com/chunlinli/intdag) which offers the peeling causal discovery method and data-perturbation/asymptotic inference.
```r
devtools::install_github("chunlinli/intdag/intdag")
```

### R Scripts of utilities

Before proceeding, source the following R scripts for utility functions required in the following illustration.

Assume the working directory is the cloned repository https://github.com/chunlinli/intdag.
```r
source("utility.R")
set.seed(1110)
```

## Examples

First, we generate a random graph with `p=10` and `q=20`.
```r
p <- 10
q <- 20
graph <- graph_generation(p = p, q = q, graph_type = "random", iv_sufficient = FALSE)
```
Note that we use the option `iv_sufficient = FALSE`, which means a crucial assumption in the paper -- Assumption 1C is violated. 
Let's print the $U$ matrix.
```r
graph$u
```
This is a `p` by `p` adjacency matrix of the DAG that we want to recover and/or make inference. Note that $U_{jk}\neq 0$ means an edge from $Y_j$ to $Y_k$.

Then print the $W$ matrix.
```r
graph$w
```
This is a `q` by `p` matrix, indicating the interventional relations of an intervention varibale $X_l$ and a primary variable $Y_j$, where $W_{lj}\neq 0$ means $X_l$ directly intervenes on $Y_j$.
In above, `w` corresponds to the simulation Setup B in the paper. 

Next, generate a random sample of size `n=200` based on the graph. According to the analysis in the paper, the distribution of intervention variables $X$ does not matter too much. Here we generate $X$ so that they have an AR(1) correlation structure. 
```r
n <- 200
x <- matrix(rnorm(n * q), nrow = n, ncol = q)
if (rho != 0) {
    for (j in 2:q) {
        x[, j] <- sqrt(1 - rho^2) * x[, j] + rho * x[, j - 1]
    }
}
```
Note that `x` is an `n` by `q` matrix.

Then we generate the $Y$ variables, the ones of primary interest. 
```r
y <- (x %*% graph$w + rmvnorm(n, sigma = diag(seq(from = 0.5, to = 1, length.out = p), p, p))) %*% solve(diag(p) - graph$u)
```
Here `y` is an `n` by `p` matrix.

Now fit the model using `intdag`. Note that `intdag` calls `glmtlp` for solving multi-response regression in its first step.
```r
v_out <- v_estimation(y = y, x = x, model_selection = "bic")
```
Second, use peeling algorithm to recover the topological layers.
```r
top_out <- topological_order(v_out$v)
an_mat <- top_out$an_mat
in_mat <- top_out$in_mat
```
Third, refit to estimate $U$ and $W$.
```r
discovery_out <- causal_discovery(y = y, x = x, an_mat = an_mat, in_mat = in_mat)
```
The final estimate for $U$ is:
```r
discovery_out$u
```
We compare the final estimate with the true graph in the structural Hamming distance.
```r
sum(abs((abs(discovery_out$u) > 0.05) - (graph$u != 0)))
```
Here we use a truncation threshold `0.05` to screen the small (noisy) coefficients.

## Contents

The R package is placed in directory `./intdag/`.

The extensive simulation code is placed in directory `./simulation/`.

## Citing information

If you find the code useful, please consider citing 
```
@article{
    author = {Chunlin Li, Xiaotong Shen, Wei Pan},
    title = {Inference for a large directed acyclic graph with unspecified interventions},
    year = {2020}
}
```
Implementing these algorithms is error-prone and this code project is still in development. 
Please file an issue if you encounter any error when using the code. I will be grateful to be informed.