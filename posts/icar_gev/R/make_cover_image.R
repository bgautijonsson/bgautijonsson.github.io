library(sf)
library(arrow)
library(dplyr)
library(tidyr)
library(bggjphd)
library(here)
library(readr)
library(tmap)
library(scales)

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


plot_dat <- gev_params |> 
  mutate(
    delta = exp(10 * delta) - 1
  ) 
  
p <- plot_dat |> 
  ggplot() +
  geom_sf(data = uk_map) +
  geom_sf(aes(fill = delta), alpha = 0.7, linewidth = 0) +
  scale_fill_viridis_c(
    labels = label_percent(),
    breaks = c(range(plot_dat$delta), 0.045, 0.07)
  ) +
  coord_sf(
    xlim = bbox[c(1, 3)],
    ylim = bbox[c(2, 4)], 
    expand = FALSE
  ) +
  labs(
    fill = "Increase per decade",
    title = "Spatial distribution of increase in maximum hourly rainfall in the UK",
    subtitle = "Shown as % increase in the location parameter of a GEV distribution"
  )

p

ggsave(
  plot = p,
  filename = here("posts", "icar_gev", "Figures", "cover.png"),
  width = 8, height = 1 * 8, scale = 1.3
)
