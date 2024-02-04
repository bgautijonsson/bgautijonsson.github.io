library(tidyverse)
library(arrow)
library(here)
library(bggjphd)
library(posterior)

d <- open_dataset(
  here("posts", "icar_gev", "posterior")
)

#### BYM Variables ####
bym_vars <- paste(
  rep(
    c(
      "logit_rho",
      "sigma",
      "mu"
    ),
    each = 4
  ),
  rep(
    c(
      "psi",
      "tau",
      "phi",
      "gamma"
    ),
    times = 3
  ),
  sep = "_"
)

inv_logit <- function(x) 1 / (1 + exp(-x))

post <- d |>
  filter(
    name %in% bym_vars
  ) |>
  collect() |>
  rename(.chain = chain) |>
  pivot_wider() |>
  mutate_at(
    vars(starts_with("logit")),
    inv_logit
  ) |>
  rename_with(
    \(x) str_replace(x, "logit_", ""),
  ) |>
  as_draws_df()


bym <- post |>
  summarise_draws()

bym |> 
  write_csv(
    here("posts", "icar_gev", "results", "bym_results.csv")
  )

#### GEV Parameters ####

gev_params <- c("psi", "tau", "phi", "gamma", "mu0", "sigma", "xi", "delta")

gev_regex <- str_c("^", gev_params, "\\[") |> 
  str_c(collapse = "|")

param_names <- here("posts", "icar_gev", "data", "new_names.csv") |> 
  read_csv() |> 
  select(new_name) |> 
  crossing(
    param = gev_params
  ) |> 
  mutate(
    string = str_c(param, "[", new_name, "]")
  )

read_gev_param <- function(gev_param) {
  
  search_string <- param_names |> 
    filter(param == gev_param) |>
    pull(string)
  
  gev_params <- d |>
    # filter(
    #   str_detect(name, "^psi\\[|^tau\\[|^phi\\[|^gamma\\[")
    # ) |>
    filter(
      name %in% search_string
    ) |>
    group_by(name) |>
    summarise(
      mean = mean(value),
      median = median(value),
      sd = sd(value)
    ) |>
    collect() |>
    mutate(
      station = parse_number(str_replace(name, "mu0", "mu")),
      variable = str_replace(name, "\\[.*\\]$", "")
    ) |>
    select(variable, station, mean)
  
  gc()
  Sys.sleep(1)
  
  gev_params
}

gev_param_results <- map_dfr(gev_params, read_gev_param)


gev_param_results |> 
  write_parquet(
    here("posts", "icar_gev", "results", "gev_results.parquet")
  )

