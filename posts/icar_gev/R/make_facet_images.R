library(sf)
library(arrow)
library(dplyr)
library(tidyr)
library(bggjphd)
library(here)
library(readr)
library(tmap)
library(scales)
library(ggplot2)
library(forcats)
library(patchwork)

theme_set(theme_bggj())

box::use(
  posts/icar_gev/R/utils[update_names]
)

gev_params <- open_dataset(
  here("posts", "icar_gev", "results", "gev_results.parquet")
)

model_stations <- read_csv(
  here("posts", "icar_gev", "data", "model_stations.csv")
)

new_names <- read_csv(
  here("posts", "icar_gev", "data", "new_names.csv")
)

gev_params <- gev_params |> 
  collect() |>
  pivot_wider(names_from = variable, values_from = mean) |> 
  inner_join(
    model_stations |> 
      update_names(station, new_names),
    by = join_by(station)
  ) 

n_x <- length(unique(gev_params$proj_x))
n_y <- length(unique(gev_params$proj_y))


gev_params <- gev_params |> 
  bggjphd::stations_to_sf() |> 
  sf::st_transform("WGS84") |> 
  bggjphd::points_to_grid(n_x = n_x, n_y = n_y)


uk_map <- get_uk_spatial(scale = "large")

bbox <- st_geometry(gev_params) |> st_bbox()


p <- gev_params |> 
  pivot_longer(c(psi:gamma), names_to = "variable") |> 
  mutate(
    variable = fct_relevel(
        variable,
        "psi", "tau", "phi", "gamma"
      )
  ) |> 
  group_by(variable2 = variable) |> 
  group_map(
    \(x, ...) {
      variable <- unique(x$variable |> as.character())
      x |> 
        ggplot() +
        geom_sf(data = uk_map) +
        geom_sf(aes(fill = value), alpha = 0.7, col = NA) +
        scale_fill_viridis_c(
          labels = label_percent()
        ) +
        coord_sf(
          xlim = bbox[c(1, 3)],
          ylim = bbox[c(2, 4)], 
          expand = FALSE
        ) +
        theme(
          legend.position = "none"
        ) +
        labs(
          subtitle = label_parse()(variable)
        )
    }
  ) |> 
  wrap_plots() +
  plot_annotation(
    title = "Spatial distribution of posterior means",
    subtitle = "GEV parameters on unconstrained scales"
  )

p

ggsave(
  plot = p,
  filename = here("posts", "icar_gev", "Figures", "facet_unconstrained.png"),
  width = 8, height = 8, scale = 1.3
)



p <- gev_params |> 
  pivot_longer(c(mu0:delta), names_to = "variable") |> 
  mutate(
    variable = fct_recode(
      variable,
      "mu[0]" = "mu0",
      "Delta" = "delta"
    ) |> 
      fct_relevel(
        "mu[0]", "sigma", "xi", "Delta"
      )
  ) |> 
  group_by(variable2 = variable) |> 
  group_map(
    \(x, ...) {
      variable <- unique(x$variable |> as.character())
      x |> 
        ggplot() +
        geom_sf(data = uk_map) +
        geom_sf(aes(fill = value), alpha = 0.7, col = NA) +
        scale_fill_viridis_c(
          labels = label_percent()
        ) +
        coord_sf(
          xlim = bbox[c(1, 3)],
          ylim = bbox[c(2, 4)], 
          expand = FALSE
        ) +
        theme(
          legend.position = "none"
        ) +
        labs(
          subtitle = label_parse()(variable)
        )
    }
  ) |> 
  wrap_plots() +
  plot_annotation(
    title = "Spatial distribution of posterior means",
    subtitle = "GEV parameters on constrained scales"
  )



ggsave(
  plot = p,
  filename = here("posts", "icar_gev", "Figures", "facet_constrained.png"),
  width = 8, height = 8, scale = 1.3
)
