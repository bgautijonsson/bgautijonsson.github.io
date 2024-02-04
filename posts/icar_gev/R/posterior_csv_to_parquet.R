library(tidyverse)
library(cmdstanr)
library(posterior)
library(arrow)
library(here)

dir.create(here("posts", "icar_gev", "posterior"))
for (chain in 1:4) {
  dir.create(here("posts", "icar_gev", "posterior", paste0("chain=", chain)))
}


runs <- here("posts", "icar_gev", "Stan", "Draws") |>
  list.files(
    full.names = TRUE
  )

files <- runs |> 
  tail(1) |> 
  list.files(
    full.names = TRUE
  )

process_output <- function(file, i) {
  d <- read_cmdstan_csv(file)$post_warmup_draws |>
    as_draws_df() |>
    as_tibble() |>
    select(-.chain) |>
    pivot_longer(c(-.iteration, -.draw)) |>
    write_parquet(
      here(
        "posts",
        "icar_gev",
        "posterior",
        paste0("chain=", i),
        "part-0.parquet"
      )
    )

  rm(d)
  gc()
  Sys.sleep(1)
}

files |>
  iwalk(process_output)

