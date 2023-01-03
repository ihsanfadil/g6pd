
# Preamble ----------------------------------------------------------------

# Author       : Ihsan Fadilah
# Email        : ifadilah@oucru.org
# Project      : Prevalence of G6PD deficiency across Indonesia
# Last Updated : 3 January 2023

# Load packages
library(tidyverse)
library(broom)
library(raster)
library(here)
library(janitor)
library(extrafont); loadfonts()

# General plot formatting
theme_set(theme_minimal())
theme_update(
  text = element_text(size = 9.5, family = "Fira Code"), # Font
  plot.title = element_text(hjust = 0),      # Centre-align title
  plot.subtitle = element_text(hjust = 0),   # Centre-align subtitle
  legend.title = element_blank(),            # Remove legend title
  # legend.position = c(0.80, 0.15),           # Move legend to bottom right
  legend.background = element_blank(),       # Remove legend background
  legend.box.background = element_blank(),   # Remove lengend-box background
  legend.spacing.y = unit(0.01, 'mm'),       # Make legend closer
  legend.key.height = unit(0.4, "cm"),       # Make legend closer
  # panel.grid.minor = element_blank(),      # Remove minor lines
  panel.grid.minor.x = element_blank(),      # Remove minor lines on the x axis
  axis.title.x = element_text(hjust = 1),    # Move title for x-axis
  axis.title.y = element_text(hjust = 0.5)   # Move title for y-axis
)
custom_colours <- c("#b37486", "#e7a29c", "#7c98a6", 
                    "#d67d53", "#5b9877", "#534d6b", "#e6bd57")

# Code --------------------------------------------------------------------
## Preprocessing ----------------------------------------------------------

# Prevalence
# Load in the data
prevalence <- readxl::read_xlsx(
  path = here('data', 'Lab 3 Data - Mapping Abstraction (1).xlsx'),
  sheet = 'PREV_DATA'
)

# Analysis
prevalence <- clean_names(prevalence) |>  # Columns in lower-case
  mutate(prev_male = 100 * male_def / n_male,
         prev_female = 100 * fem_def / n_female,
         site_name = if_else(site_name == '2 cities in South Kalimantan',
                             'Banjarmasin & Banjarbaru',
                             site_name),
         site_name = if_else(site_name == 'Padira Tana' & year_start == 2012,
                             'Padira Tana (2012)', site_name),
         site_name = if_else(site_name == 'Padira Tana' & year_start == 2016,
                             'Padira Tana (2016)', site_name),
         island = case_when(
           admin1 == 'Bangka Belitung Islands' ~ 'Bangka Belitung Islands',
           admin1 == 'Bengkulu' ~ 'Sumatra',
           admin1 == 'Central Kalimantan' ~ 'Kalimantan',
           admin1 == 'East Nusa Tenggara' ~ 'East Nusa Tenggara',
           admin1 == 'Jambi' ~ 'Sumatra',
           admin1 == 'Maluku' ~ 'Maluku',
           admin1 == 'North Maluku' ~ 'North Maluku',
           admin1 == 'Papua' ~ 'Papua',
           admin1 == 'South Kalimantan' ~ 'Kalimantan',
           admin1 == 'West Sumatra' ~ 'Sumatra',
           TRUE ~ NA_character_
           ),
         pm = prev_male / 100,
         pf = prev_female / 100,
         lci_male = 100 * (pm - (qnorm(0.975) * sqrt((pm * (1 - pm))/n_male))),
         lci_male = if_else(lci_male < 0, 0, lci_male),
         uci_male = 100 * (pm + (qnorm(0.975) * sqrt((pm * (1 - pm))/n_male))),       
         )

# CI for males
ci_male <- prevalence |> 
  drop_na(prev_male) |> 
  group_by(site_name) |> 
  nest() |> 
  mutate(out = map(
    .x = data,
    .f = ~summarise(.data = .x,
                    out = list(prop.test(male_def, n_male) |> tidy()))
  )) |> 
  ungroup() |> 
  unnest(out) |> 
  unnest(cols = c(data, out)) |> 
  mutate(estimate = 100 * estimate,
         conf.low = 100 * conf.low,
         conf.high = 100 * conf.high)

# CI for females
ci_female <- prevalence |> 
  drop_na(prev_female) |> 
  group_by(site_name) |> 
  nest() |> 
  mutate(out = map(
    .x = data,
    .f = ~summarise(.data = .x,
                    out = list(prop.test(fem_def, n_female) |> tidy()))
  )) |> 
  ungroup() |> 
  unnest(out) |> 
  unnest(cols = c(data, out)) |> 
  mutate(estimate = 100 * estimate,
         conf.low = 100 * conf.low,
         conf.high = 100 * conf.high)

# Map
# Get the level-2 (district-level) Indonesia-map data
ina_data <- getData('GADM', country = 'IDN', level = 1)

# Convert the data into ggplot-friendly data
ina_data_fortified <- fortify(ina_data)


## Plot: Points -----------------------------------------------------------
point_male <- ci_male |> 
  drop_na(prev_male) |> # 3 with NA
  filter(n_male >= 5) |> # Pund 0/2
  mutate(site_name = factor(site_name),
         site_name = fct_reorder(site_name, estimate)) |> 
  ggplot(aes(x = estimate, y = site_name)) +
    geom_point(colour = "#5b9877", alpha = 0.7, size = 1.8, shape = 19) +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), colour = "#5b9877",
                  alpha = 0.15, size = 1.5, linetype = 1, width = 0.5) +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)) +
    labs(x = "\nEstimated prevalence (%) in males",
         y = "")

point_male

point_female <- ci_female |> 
  drop_na(prev_female) |> # 7 with NA
  filter(n_female >= 5) |> # None excluded
  mutate(site_name = factor(site_name),
         site_name = fct_reorder(site_name, estimate)) |> 
  ggplot(aes(x = estimate, y = site_name)) +
  geom_point(colour = "#b37486", alpha = 0.7, size = 1.8, shape = 19) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), colour = "#b37486",
                alpha = 0.15, size = 1.5, linetype = 1, width = 0.5) +
  scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)) +
  labs(x = "\nEstimated prevalence (%) in females",
       y = "")

point_female

## Plot: Trend ------------------------------------------------------------
trend_summary <- prevalence |> 
  mutate(year = factor(year_start)) |> 
  group_by(year, island) |> 
  count()

trend <- trend_summary |> 
  ggplot(aes(x = year, y = n)) +
    geom_col(aes(fill = island)) +
    scale_y_continuous(limits = c(0, 16), breaks = seq(0, 16, by = 2)) +
    scale_fill_manual(values = custom_colours) +
    theme(legend.position = c(0.75, 0.80)) +
    labs(x = "\nYear",
         y = "Number of published prevalence studies\n")

trend

## Plot: Map --------------------------------------------------------------
ina_map <- ina_data_fortified |> 
  ggplot(aes(x = long, y = lat, group = group)) +
  geom_polygon(aes(fill = 'lightgray', colour = 'white'), size = 1) +
  theme_void()

ina_map












