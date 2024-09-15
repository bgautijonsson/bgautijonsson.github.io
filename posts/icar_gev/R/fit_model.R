if (interactive()) {
  box::use(
    posts / icar_gev / r / prep_data
  )
} else {
  box::use(
    . / prep_data
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




here("posts", "icar_gev", "results", "bym_results.csv") |>
  read_csv() |>
  mutate(
    variable = str_c("<b>&", variable, ";</sub></b>") |>
      str_replace("_", ";<sub>&")
  ) |>
  arrange(rep(1:3, times = 4)) |>
  gt() |>
  cols_hide(
    columns = c(mean, sd, mad, rhat, ess_bulk, ess_tail)
  ) |>
  cols_label(
    variable = "Parameter",
    median = "Median",
    q5 = "5th Percentile",
    q95 = "95th Percentile"
  ) |>
  tab_header(
    title = "BYM2 hyperparameters"
  ) |>
  fmt_markdown(
    columns = variable
  ) |>
  fmt_number(
    columns = mean:rhat,
    decimals = 3
  ) |>
  fmt_number(
    columns = ess_bulk:ess_tail,
    decimals = 0
  ) |>
  tab_options(
    table.font.size = 20,
    table.width = pct(100)
  )
