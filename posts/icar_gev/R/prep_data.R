#' @export
prep_data <- function(x_range = c(100, 150), y_range = c(100, 150)) {
  if (interactive()) {
    box::use(
      posts/icar_gev/r/utils
    )
  } else {
    box::use(
      ./utils
    )
  }
  
  box::use(
    dplyr[filter, between, mutate, select, row_number, distinct, inner_join, semi_join, join_by],
    tidyr[pivot_wider],
    tibble[column_to_rownames],
    readr[write_csv],
    here[here],
    bggjphd[stations, precip, twelve_neighbors]
  )
  model_stations <- stations |>
    filter(
      between(proj_x, x_range[1], x_range[2]),
      between(proj_y, y_range[1], y_range[2])
    )
  
  model_stations |> 
    write_csv(
      here("posts", "icar_gev", "data", "model_stations.csv")
    )
  
  new_names <- model_stations |> 
    mutate(new_name = row_number()) |> 
    distinct(station, new_name)
  
  new_names |> 
    write_csv(
      here("posts", "icar_gev", "data", "new_names.csv")
    )
  
  model_precip <- precip |>
    semi_join(
      model_stations,
      by = join_by(station)
    )
  
  
  precip_matrix <- model_precip |>
    pivot_wider(names_from = station, values_from = precip) |>
    column_to_rownames("year") |>
    as.matrix()
  
  N_stations <- ncol(precip_matrix)
  N_years <- nrow(precip_matrix)
  
  
  edges <- twelve_neighbors |>
    filter(
      type %in% c("e", "n", "w", "s")
    ) |>
    inner_join(
      model_stations,
      by = join_by(station)
    ) |>
    semi_join(
      model_stations,
      by = join_by(neighbor == station)
    ) |>
    select(station, neighbor) |> 
    utils$update_names(table = _, variable = station, new_names = new_names) |> 
    utils$update_names(table = _, variable = neighbor, new_names = new_names)
  
  N_neighbors = nrow(edges)
  node1 <- edges$station
  node2 <- edges$neighbor
  scaling_factor <- utils$get_scaling_factor(edges, model_stations)
  
  stan_data <- list(
    N_stations = N_stations,
    N_years = N_years,
    precip = precip_matrix,
    N_neighbors = N_neighbors,
    node1 = node1,
    node2 = node2,
    scaling_factor = scaling_factor
  )
  
  stan_data
}