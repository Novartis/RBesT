#'
#' # SBC Graphical Results
#'
#'

#'
#' ## Dense Scenario, $\mu$
#'
print(pl_dense$mu$hist)
print(pl_dense$mu$ecdf_diff)


#'
#' ## Dense Scenario, $\tau$
#'

print(pl_dense$tau$hist)
print(pl_dense$tau$ecdf_diff)

#'
#' ## Sparse Scenario, $\mu$
#'

print(pl_sparse$mu$hist)
print(pl_sparse$mu$ecdf_diff)

#'
#' ## Sparse Scenario, $\tau$
#'

print(pl_sparse$tau$hist)
print(pl_sparse$tau$ecdf_diff)

