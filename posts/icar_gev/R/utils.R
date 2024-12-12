#' @export
update_names <- function(table, variable, new_names) {
  box::use(
    dplyr[inner_join, join_by, mutate, select]
  )
  table |> 
    inner_join(
      new_names,
      by = join_by({{ variable }} == station)
    ) |> 
    mutate(
      "{{variable}}" := new_name
    ) |> 
    select(-new_name)
}

#' @export
get_scaling_factor <- function(edges, N_stations) {
  box::use(
    dplyr[filter, rename],
    Matrix[sparseMatrix, Diagonal, rowSums, diag],
    INLA[inla.qinv]
  )
  
  nbs <- edges |>
    filter(neighbor > station) |>
    rename(node1 = station, node2 = neighbor)
  
  N <- nrow(model_stations)
  
  adj.matrix <- sparseMatrix(i = nbs$node1, j = nbs$node2, x = 1, symmetric = TRUE)
  # The ICAR precision matrix (note! This is singular)
  Q <- Diagonal(N, rowSums(adj.matrix)) - adj.matrix
  # Add a small jitter to the diagonal for numerical stability (optional but recommended)
  Q_pert <- Q + Diagonal(N) * max(diag(Q)) * sqrt(.Machine$double.eps)
  
  # Compute the diagonal elements of the covariance matrix subject to the
  # constraint that the entries of the ICAR sum to zero.
  # See the inla.qinv function help for further details.
  Q_inv <- inla.qinv(Q_pert, constr = list(A = matrix(1, 1, N), e = 0))
  
  # Compute the geometric mean of the variances, which are on the diagonal of Q.inv
  scaling_factor <- exp(mean(log(diag(Q_inv))))
  
  scaling_factor
}