if (interactive()) {
  box::use(
    posts/icar_gev/r/prep_data
  )
} else {
  box::use(
    ./prep_data
  )
}



library(cmdstanr)
library(here)

output_dir <- here("posts", "icar_gev", "Stan", "Draws", Sys.time())
dir.create(output_dir)

stan_data <- prep_data$prep_data(
  x_range = c(50, 150),
  y_range = c(100, 200)
)

model <- cmdstan_model(
  here("posts", "icar_gev", "Stan", "BYM_GEV_LOGIT.stan")
)

fit <- model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  save_warmup = FALSE,
  refresh = 10,
  output_dir = output_dir
)
