
# Preamble ----------------------------------------------------------------

# Author       : Ihsan Fadilah
# Email        : ifadilah@oucru.org
# Project      : Prevalence of G6PD deficiency across Indonesia

# First created : 
# Last updated  : 

# Package dependencies
library(tidyverse)
library(ggbeeswarm)
library(broom)
library(raster)
library(here)
library(janitor)
library(sf)
library(rmapshaper)
library(rgdal)
library(cowplot)
library(extrafont); loadfonts()

# General plot formatting
theme_set(theme_minimal())
theme_update(
  text = element_text(size = 9.5, family = "Fira Code"), # Font
  plot.title = element_text(hjust = 0),      # Centre-align title
  plot.subtitle = element_text(hjust = 0),   # Centre-align subtitle
  legend.title = element_blank(),            # Remove legend title
  # legend.position = c(0.80, 0.15),         # Move legend to bottom right
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
         prev_male_0.7 = 100 * (male_int + male_def) / n_male,
         prev_female_0.7 = 100 * (fem_int + fem_def) / n_female,
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
         n_male_cat = case_when(
            n_male <= 10 ~ '[1, 10)',
            n_male <= 100 ~ '[10, 100)',
            n_male <= 1000 ~ '[100, 1000)',
            is.na(n_male) ~ NA_character_,
            TRUE ~ 'Check me!'
         ) |> factor(), 
         n_female_cat = case_when(
           n_female <= 10 ~ '[1, 10)',
           n_female <= 100 ~ '[10, 100)',
           n_female <= 1000 ~ '[100, 1000)',
           is.na(n_female) ~ NA_character_,
           TRUE ~ 'Check me!'
         ) |> factor(), 
         n_total_cat = case_when(
           n_total <= 10 ~ '[1, 10)',
           n_total <= 100 ~ '[10, 100)',
           n_total <= 1000 ~ '[100, 1000)',
           is.na(n_total) ~ NA_character_,
           TRUE ~ 'Check me!'
         ) |> factor(levels = c('[1, 10)', '[10, 100)', '[100, 1000)')), 
         prev_total = 100 * (total_def / n_total)
  )

# CI for males (assuming iid)
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
  dplyr::select(-prev_male) 
  
ci_male <- ci_male |> 
  mutate(estimate = 100 * estimate,
         conf.low = 100 * conf.low,
         conf.high = 100 * conf.high,
         sex = rep('Male', dim(ci_male)[1]) |> factor())

# <0.7 activity
ci_male_0.7 <- prevalence |> 
  drop_na(prev_male_0.7) |> 
  group_by(site_name) |> 
  nest() |> 
  mutate(out = map(
    .x = data,
    .f = ~summarise(.data = .x,
                    out = list(prop.test((male_int + male_def), n_male) |>
                                 tidy()))
  )) |> 
  ungroup() |> 
  unnest(out) |> 
  unnest(cols = c(data, out)) |> 
  dplyr::select(-prev_male_0.7) 

ci_male_0.7 <- ci_male_0.7 |> 
  mutate(estimate = 100 * estimate,
         conf.low = 100 * conf.low,
         conf.high = 100 * conf.high,
         sex = rep('Male', dim(ci_male_0.7)[1]) |> factor())

# CI for females (assuming iid)
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
  dplyr::select(-prev_female)

ci_female <- ci_female |> 
  mutate(estimate = 100 * estimate,
         conf.low = 100 * conf.low,
         conf.high = 100 * conf.high,
         sex = rep('Female', dim(ci_female)[1]) |> factor())

# <0.7 activity
ci_female_0.7 <- prevalence |> 
  drop_na(prev_female_0.7) |> 
  group_by(site_name) |> 
  nest() |> 
  mutate(out = map(
    .x = data,
    .f = ~summarise(.data = .x,
                    out = list(prop.test((fem_int + fem_def), n_female) |>
                                 tidy()))
  )) |> 
  ungroup() |> 
  unnest(out) |> 
  unnest(cols = c(data, out)) |> 
  dplyr::select(-prev_female_0.7) 

ci_female_0.7 <- ci_female_0.7 |> 
  mutate(estimate = 100 * estimate,
         conf.low = 100 * conf.low,
         conf.high = 100 * conf.high,
         sex = rep('Female', dim(ci_female_0.7)[1]) |> factor())

# Central tendencies by sex (weights approximated by sample sizes)
# Males
central_male <- ci_male |> 
  mutate(prevxn = estimate * n_male) |> 
  summarise(weighted_mean = sum(prevxn, na.rm = T) / sum(n_male, na.rm = T),
            median = median(estimate)) 

weighted_mean_male <- central_male |> pull(weighted_mean)
median_male <- central_male |> pull(median)

# <0.7
central_male_0.7 <- ci_male_0.7 |> 
  mutate(prevxn = estimate * n_male) |> 
  summarise(weighted_mean = sum(prevxn, na.rm = T) / sum(n_male, na.rm = T),
            median = median(estimate)) 

weighted_mean_male_0.7 <- central_male_0.7 |> pull(weighted_mean)
median_male_0.7 <- central_male_0.7 |> pull(median)

# Females
central_female <- ci_female |> 
  mutate(prevxn = estimate * n_female) |> 
  summarise(weighted_mean = sum(prevxn, na.rm = T) / sum(n_female, na.rm = T),
            median = median(estimate))

weighted_mean_female <- central_female |> pull(weighted_mean)
median_female <- central_female |> pull(median)

# <0.7
central_female_0.7 <- ci_female_0.7 |> 
  mutate(prevxn = estimate * n_female) |> 
  summarise(weighted_mean = sum(prevxn, na.rm = T) / sum(n_female, na.rm = T),
            median = median(estimate))

weighted_mean_female_0.7 <- central_female_0.7 |> pull(weighted_mean)
median_female_0.7 <- central_female_0.7 |> pull(median)

# Bees
est_male <- ci_male |> dplyr::select(sex, estimate, island)
est_female <- ci_female |> dplyr::select(sex, estimate, island)
estimate <- add_row(est_male, est_female)

# <0.7
est_male_0.7 <- ci_male_0.7 |> dplyr::select(sex, estimate, island)
est_female_0.7 <- ci_female_0.7 |> dplyr::select(sex, estimate, island)
estimate_0.7 <- add_row(est_male_0.7, est_female_0.7)

# Sample sizes
sample_male <- ci_male |>
  dplyr::select(sex, n_male) |>
  rename(n = n_male)
sample_female <- ci_female |>
  dplyr::select(sex, n_female) |> 
  rename(n = n_female)
sample <- add_row(sample_male, sample_female)

# Map
# Load the map data (run if complete rerunning is required)
# ina_shp <- readOGR(dsn = here('data', 'INDOPROV2017', 'INDOPROV2017.shp'),
#                    stringsAsFactors = F)
 
# Simplify the high-resolution map data to 0.1% of the original
# Run if complete rerunning is required
# ina_shp_compressed <- ms_simplify(ina_shp, keep = 0.001, keep_shapes = TRUE)

# Save the map object to cut processing time
# Run if complete rerunning is required
# save(ina_shp_compressed, file = here('data', 'ina_shp_compressed.RData'))

# Load the simplified map data
load(here('data', 'ina_shp_compressed.RData'))

# API and population-size data
population <- readxl::read_xlsx(
  path = here('data', 'Data Penduduk Malaria_2017-2022.xlsx'),
  sheet = 'Sheet1'
) |> 
  clean_names() |> 
  dplyr::select(province, district, population_2020, population_2021) |> 
  group_by(province) |> 
  nest() |> 
  mutate(
    pop_2020 = map(
      .x = data,
      .f = ~summarise(.data = .x, pop_2020 = sum(population_2020))),
    pop_2021 = map(
      .x = data,
      .f = ~summarise(.data = .x, pop_2021 = sum(population_2021)))
  ) |> 
  ungroup() |> 
  unnest(cols = c(pop_2020, pop_2021)) |> 
  drop_na(province) |> 
  rename(data_pop = data)

api <- readxl::read_xlsx(
  path = here('data', 'Data endemisitas 2022.xlsx'),
  sheet = 'Endemisitas'
) |> 
  clean_names() |> 
  dplyr::select(propinsi, kabupaten, positif_2020, positif_2021_per_13_mar22) |> 
  rename(province = propinsi,
         district = kabupaten,
         n_case_2020 = positif_2020,
         n_case_2021 = positif_2021_per_13_mar22) |> 
  group_by(province) |> 
  nest() |> 
  mutate(
    n_case_2020 = map(
      .x = data,
      .f = ~summarise(.data = .x, n_case_2020 = sum(n_case_2020))),
    n_case_2021 = map(
      .x = data,
      .f = ~summarise(.data = .x, n_case_2021 = sum(n_case_2021)))
  ) |> 
  ungroup() |> 
  unnest(cols = c(n_case_2020, n_case_2021)) |> 
  drop_na(n_case_2020) |> 
  filter(province != 'NASIONAL') |> 
  rename(data_api = data)

api_pop <- full_join(api, population, by = 'province') |> 
  mutate(
    api_2020 = 1000 * n_case_2020 / pop_2020,
    api_2021 = 1000 * n_case_2021 / pop_2021,
    api_wavg = ((api_2020 * pop_2020) + (api_2021 * pop_2021)) /
               (pop_2020 + pop_2021),
    api_avg = (api_2020 + api_2021) / 2,
    api_group = case_when(api_wavg >= 100 ~ 'High III [100, ∞)',
                          api_wavg >= 50 ~ 'High II [50, 100)',
                          api_wavg >= 5 ~ 'High I [5, 50)',
                          api_wavg >= 1 ~ 'Moderate [1, 5)',
                          TRUE ~ 'Low or elimination [0, 1)'),
    api_group = factor(api_group,
                       ordered = TRUE,
                       levels = c('High III [100, ∞)',
                                  'High II [50, 100)',
                                  'High I [5, 50)',
                                  'Moderate [1, 5)',
                                  'Low or elimination [0, 1)'))
  )

api_pop_for_map <- api_pop |>
  dplyr::select(province, api_wavg, api_group) |> 
  rename(ADMIN1 = province) |> 
  mutate(ADMIN1 = if_else(ADMIN1 == 'KEP. RIAU', 'KEPULAUAN RIAU', ADMIN1))

ina_shp_compressed@data <- full_join(ina_shp_compressed@data,
                                     api_pop_for_map,
                                     by = 'ADMIN1') |> 
  mutate(id = IDADMIN1)

shp_df <- broom::tidy(ina_shp_compressed, region = "id")
shp_df <- shp_df %>% left_join(ina_shp_compressed@data, by = c("id" = "id"))

# District-specific r(prevalence, api)
# district_prevalence <- prevalence |> 
#   dplyr::select(admin2, site_name,
#                 prev_male, prev_female, prev_male_0.7, prev_female_0.7) |> 
#   mutate(admin2 = str_to_lower(admin2)) |> 
#   rename(district = admin2)
# 
# district_api <- readxl::read_xlsx(
#   path = here('data', 'Data endemisitas 2022.xlsx'),
#   sheet = 'Endemisitas'
# ) |>
#   clean_names() |> 
#   dplyr::select(propinsi, kabupaten, api_2018) |> 
#   rename(district = kabupaten,
#          province = propinsi) |> 
#   mutate(district = str_to_lower(district),
#          province = str_to_lower(province))
#  
# # Manually calculate weighted-API for Banjarbaru & Banjarmasin (2018)
# api_banjar <- ((0.093897816 * 255597) + (0.008560801 * 700869)) /
#               (255597 + 700869)
# 
# district_prev_api <- full_join(district_prevalence,
#                                district_api,
#                                by = 'district') |> 
#   drop_na(site_name)
# 
# district_prev_api[13, 8] <- api_banjar
# district_prev_api[13, 7] <- 'kalimantan selatan'
# 
# x <- district_prev_api %>%
#   pivot_longer(
#     cols = starts_with("prev_"),
#     names_to = "sex",
#     names_prefix = "prev_",
#     values_to = "prev_def",
#     values_drop_na = FALSE
#   )
# 
# district_prev_api_def <- x |>
#   filter(sex %in% c('male', 'female')) |> 
#   mutate(sex = if_else(sex == 'male', 'Male', 'Female')) |> 
#   drop_na(prev_def)
# 
# plot(filter(district_prev_api_def, sex == 'Female')$api_2018, 
#     filter(district_prev_api_def, sex == 'Female')$prev_def)
# 
# district_prev_api_0.7 <- x |>
#   filter(sex %in% c('male_0.7', 'female_0.7')) |> 
#   mutate(sex = if_else(sex == 'male_0.7', 'Male', 'Female')) |> 
#   drop_na(prev_def)
# 
# cor(filter(district_prev_api_0.7, sex == 'Female')$api_2018, 
#     filter(district_prev_api_0.7, sex == 'Female')$prev_def)

## Plot: Points -----------------------------------------------------------
point_male <- ci_male |> 
  # drop_na(estimate) |> # 3 with NA
  filter(n_male >= 5) |> # Pund 0/2
  mutate(site_name = factor(site_name),
         site_name = fct_reorder(site_name, estimate)) |> 
  ggplot(aes(x = estimate, y = site_name)) +
  # geom_vline(xintercept = weighted_mean_male,
  #           linetype = 'dashed', colour = 'gray80', size = 0.4) +
  geom_point(colour = "#5b9877", alpha = 0.8, size = 1.8, shape = 19) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), colour = "#5b9877",
                alpha = 0.3, size = 1.5, linetype = 1, width = 0.5) +
  scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 10)) +
  labs(x = "\nPrevalence (%)",
       y = "")
# + facet_wrap(vars(factor(island)), scales = "free")

point_male_0.7 <- ci_male_0.7 |> 
  # drop_na(estimate) |> # 3 with NA
  filter(n_male >= 5) |> # Pund 0/2
  mutate(site_name = factor(site_name),
         site_name = fct_reorder(site_name, estimate)) |> 
  ggplot(aes(x = estimate, y = site_name)) +
  # geom_vline(xintercept = weighted_mean_male,
  #           linetype = 'dashed', colour = 'gray80', size = 0.4) +
  geom_point(colour = "#5b9877", alpha = 0.8, size = 1.8, shape = 19) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), colour = "#5b9877",
                alpha = 0.3, size = 1.5, linetype = 1, width = 0.5) +
  scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 10)) +
  labs(x = "\nPrevalence (%)",
       y = "")
# + facet_wrap(vars(factor(island)), scales = "free")

point_male
point_male_0.7

point_female <- ci_female |> 
  # drop_na(prev_female) |> # 7 with NA
  filter(n_female >= 5) |> # None excluded
  mutate(site_name = factor(site_name),
         site_name = fct_reorder(site_name, estimate)) |> 
  ggplot(aes(x = estimate, y = site_name)) +
  # geom_vline(xintercept = weighted_mean_female,
  #            linetype = 'dashed', colour = 'gray80', size = 0.4) +
  geom_point(colour = "#b37486", alpha = 0.8, size = 1.8, shape = 19) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), colour = "#b37486",
                alpha = 0.3, size = 1.5, linetype = 1, width = 0.5) +
  scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 10)) +
  labs(x = "\nPrevalence (%)",
       y = "")

point_female_0.7 <- ci_female_0.7 |> 
  # drop_na(prev_female) |> # 7 with NA
  filter(n_female >= 5) |> # None excluded
  mutate(site_name = factor(site_name),
         site_name = fct_reorder(site_name, estimate)) |> 
  ggplot(aes(x = estimate, y = site_name)) +
  # geom_vline(xintercept = weighted_mean_female,
  #            linetype = 'dashed', colour = 'gray80', size = 0.4) +
  geom_point(colour = "#b37486", alpha = 0.8, size = 1.8, shape = 19) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), colour = "#b37486",
                alpha = 0.3, size = 1.5, linetype = 1, width = 0.5) +
  scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 10)) +
  labs(x = "\nPrevalence (%)",
       y = "")

point_female
point_female_0.7


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
  labs(x = "\nYear conducted",
       y = "Number of prevalence studies\n")

trend

## Plot: Map --------------------------------------------------------------
ina_map <- shp_df |> 
  ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group, fill = api_group),
                 colour = 'white', size = 0.2, alpha = 1) +
    theme_void() +
    theme(legend.key.size = unit(1, 'lines'),
          legend.spacing.x = unit(0.3, 'lines'),
          plot.margin = margin(5, 5, 5, 5, "mm"),
          text = element_text(size = 9, family = "Fira Code"),
          legend.title = element_text(size = 7.5)) +
    labs(fill =
      '\nProvince-specific\nannual parasite incidence\nper 1000 population') +
    scale_fill_manual(values = c('#a84268',
                                 '#e3615f',
                                 '#fcb97d',
                                 '#9dbf9e'),
                      limits = c('High II [50, 100)',
                                 'High I [5, 50)',
                                 'Moderate [1, 5)',
                                 'Low or elimination [0, 1)'))

ina_map

map_female <- ina_map +
    geom_point(data = ci_female, pch = 20, alpha = 0.7,
               aes(x = long, y = lat, size = n_female_cat, colour = estimate)) +
    scale_colour_gradient(low = "#DB1F48", high = "#000000", na.value = NA,
                        limits = c(0, 50)) +
#    scale_size_continuous(range = c(3, 7),
#                          limits = c(0, 1000),
#                          breaks = seq(100, 1000, by = 200)) +
    # guides(size = guide_legend(title = "Sample size (n)")) +
    theme(legend.position = "right",
          legend.direction = "vertical",
          legend.key.size = unit(1, 'lines'),
          legend.spacing.x = unit(0.3, 'lines'),
          plot.margin = margin(5, 10, 5, 5, "mm"),
          text = element_text(size = 9, family = "Fira Code"),
          legend.title = element_text(size = 9, family = "Fira Code Bold")) +
    labs(size = '\nSample size',
         colour = '\nPrevalence (%)\n') +
  geom_rect(
    xmin = (119.779369 - 1.2),
    ymin = (-9.657382 - 1.2),
    xmax = (119.779369 + 1.2),
    ymax = (-9.657382 + 1.2),
    fill = NA, 
    colour = "black",
    size = 0.3
  )

map_female

inset_map_female <- ggdraw(map_female) +
  draw_plot(
    {
      map_female +
        coord_sf(
          xlim = c((119.779369 - 1.2), (119.779369 + 1.2)),
          ylim = c((-9.657382 - 1.2), (-9.657382 + 1.2)),
          expand = FALSE
        ) +
        theme(legend.position = 'none')
    },
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    x = 0.58, 
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    y = 0.60,
    # The width and height of the plot expressed as proportion of the entire
    # ggdraw object
    width = 0.38, 
    height = 0.38)

inset_map_female

map_male <- ina_map +
  geom_point(data = ci_male, pch = 20, alpha = 0.7,
             aes(x = long, y = lat, size = n_male_cat, colour = estimate)) +
  scale_colour_gradient(low = "#DB1F48", high = "#000000", na.value = NA,
                        limits = c(0, 50)) +
#  scale_size_continuous(range = c(3, 7),
#                        limits = c(0, 1000),
#                        breaks = seq(100, 1000, by = 200)) +
  # guides(size = guide_legend(title = "Sample size (n)")) +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.key.size = unit(1, 'lines'),
        legend.spacing.x = unit(0.3, 'lines'),
        plot.margin = margin(5, 10, 5, 5, "mm"),
        text = element_text(size = 9, family = "Fira Code"),
        legend.title = element_text(size = 9, family = "Fira Code Bold")) +
  labs(size = '\nSample size',
       colour = '\nPrevalence (%)\n') +
  geom_rect(
    xmin = (119.779369 - 1.2),
    ymin = (-9.657382 - 1.2),
    xmax = (119.779369 + 1.2),
    ymax = (-9.657382 + 1.2),
    fill = NA, 
    colour = "black",
    size = 0.3
  )

map_male

inset_map_male <- ggdraw(map_male) +
  draw_plot(
    {
      map_male +
        coord_sf(
          xlim = c((119.779369 - 1.2), (119.779369 + 1.2)),
          ylim = c((-9.657382 - 1.2), (-9.657382 + 1.2)),
          expand = FALSE
        ) +
        theme(legend.position = 'none')
    },
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    x = 0.58, 
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    y = 0.60,
    # The width and height of the plot expressed as proportion of the entire
    # ggdraw object
    width = 0.38, 
    height = 0.38)

inset_map_male

# <0.7
map_female_0.7 <- ina_map +
  geom_point(data = ci_female_0.7, pch = 20, alpha = 0.7,
             aes(x = long, y = lat, size = n_female_cat, colour = estimate)) +
  scale_colour_gradient(low = "#DB1F48", high = "#000000", na.value = NA,
                        limits = c(0, 50)) +
  # scale_size_continuous(range = c(3, 7),
  #                       limits = c(0, 1000),
  #                       breaks = seq(100, 1000, by = 200)) +
  # guides(size = guide_legend(title = "Sample size (n)")) +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.key.size = unit(1, 'lines'),
        legend.spacing.x = unit(0.3, 'lines'),
        plot.margin = margin(5, 10, 5, 5, "mm"),
        text = element_text(size = 9, family = "Fira Code"),
        legend.title = element_text(size = 9, family = "Fira Code Bold")) +
  labs(size = '\nSample size',
       colour = '\nPrevalence (%)\n') +
  geom_rect(
    xmin = (119.779369 - 1.2),
    ymin = (-9.657382 - 1.2),
    xmax = (119.779369 + 1.2),
    ymax = (-9.657382 + 1.2),
    fill = NA, 
    colour = "black",
    size = 0.3
  )

map_female_0.7

inset_map_female_0.7 <- ggdraw(map_female_0.7) +
  draw_plot(
    {
      map_female_0.7 +
        coord_sf(
          xlim = c((119.779369 - 1.2), (119.779369 + 1.2)),
          ylim = c((-9.657382 - 1.2), (-9.657382 + 1.2)),
          expand = FALSE
        ) +
        theme(legend.position = 'none')
    },
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    x = 0.58, 
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    y = 0.60,
    # The width and height of the plot expressed as proportion of the entire
    # ggdraw object
    width = 0.38, 
    height = 0.38)

inset_map_female_0.7

# map_male_0.7 <- ina_map +
#   geom_point(data = ci_male_0.7, pch = 20, alpha = 0.7,
#              aes(x = long, y = lat, size = n_male, colour = estimate)) +
#   scale_colour_gradient(low = "#DB1F48", high = "#000000", na.value = NA,
#                         limits = c(0, 50)) +
#   scale_size_continuous(range = c(3, 7),
#                         limits = c(0, 1000),
#                         breaks = seq(100, 1000, by = 200)) +
#   # guides(size = guide_legend(title = "Sample size (n)")) +
#   theme(legend.position = "right",
#         legend.direction = "vertical",
#         legend.key.size = unit(1, 'lines'),
#         legend.spacing.x = unit(0.3, 'lines'),
#         plot.margin = margin(5, 10, 5, 5, "mm"),
#         text = element_text(size = 9, family = "Fira Code"),
#         legend.title = element_text(size = 9, family = "Fira Code Bold")) +
#   labs(size = '\nSample size',
#        colour = '\nPrevalence (%)\n') +
#   geom_rect(
#     xmin = (119.779369 - 1.2),
#     ymin = (-9.657382 - 1.2),
#     xmax = (119.779369 + 1.2),
#     ymax = (-9.657382 + 1.2),
#     fill = NA, 
#     colour = "black",
#     size = 0.3
#   )
# 
# map_male_0.7
# 
# inset_map_male_0.7 <- ggdraw(map_male_0.7) +
#   draw_plot(
#     {
#       map_male_0.7 +
#         coord_sf(
#           xlim = c((119.779369 - 1.2), (119.779369 + 1.2)),
#           ylim = c((-9.657382 - 1.2), (-9.657382 + 1.2)),
#           expand = FALSE
#         ) +
#         theme(legend.position = 'none')
#     },
#     # The distance along a (0,1) x-axis to draw the left edge of the plot
#     x = 0.58, 
#     # The distance along a (0,1) y-axis to draw the bottom edge of the plot
#     y = 0.60,
#     # The width and height of the plot expressed as proportion of the entire
#     # ggdraw object
#     width = 0.38, 
#     height = 0.38)
# 
# inset_map_male_0.7

# Overall (male + female) prevalence
# Deficient
map_total <- ina_map +
  geom_point(data = prevalence, pch = 20, alpha = 0.7,
             aes(x = long, y = lat, size = n_total_cat, colour = prev_total)) +
  scale_colour_gradient(low = "#DB1F48", high = "#000000", na.value = NA,
                        limits = c(0, 50)) +
  scale_size_discrete(drop = FALSE) +
  #  scale_size_continuous(range = c(3, 7),
  #                        limits = c(0, 1000),
  #                        breaks = seq(100, 1000, by = 200)) +
  # guides(size = guide_legend(title = "Sample size (n)")) +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.key.size = unit(1, 'lines'),
        legend.spacing.x = unit(0.3, 'lines'),
        plot.margin = margin(5, 10, 5, 5, "mm"),
        text = element_text(size = 9, family = "Fira Code"),
        legend.title = element_text(size = 9, family = "Fira Code Bold")) +
  labs(size = '\nSample size',
       colour = '\nPrevalence (%)\n') +
  geom_rect(
    xmin = (119.779369 - 1.2),
    ymin = (-9.657382 - 1.2),
    xmax = (119.779369 + 1.2),
    ymax = (-9.657382 + 1.2),
    fill = NA, 
    colour = "black",
    size = 0.3
  )

map_total

inset_map_total <- ggdraw(map_total) +
  draw_plot(
    {
      map_total +
        coord_sf(
          xlim = c((119.779369 - 1.2), (119.779369 + 1.2)),
          ylim = c((-9.657382 - 1.2), (-9.657382 + 1.2)),
          expand = FALSE
        ) +
        theme(legend.position = 'none')
    },
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    x = 0.58, 
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    y = 0.60,
    # The width and height of the plot expressed as proportion of the entire
    # ggdraw object
    width = 0.38, 
    height = 0.38)

inset_map_total

# <0.7?

## Plot: Beeswarm ----------------------------------------------------------

est_by_sex <- estimate |> 
  ggplot() +
    # geom_hline(yintercept = weighted_mean_male,
    #            colour = "#5b9877", linetype = 2) +
    # geom_hline(yintercept = weighted_mean_female,
    #            colour = "#e7a29c", linetype = 2) +
    # geom_hline(yintercept = median_male,
    #            colour = "#5b9877", linetype = 3) +
    # geom_hline(yintercept = median_female,
    #            colour = "#e7a29c", linetype = 3) +
    geom_beeswarm(aes(x = sex, y = estimate, colour = sex),
                  alpha = 0.6, size = 3.1, cex = 2.3) +
    geom_boxplot(aes(x = sex, y = estimate),
                 alpha = 0, width = 0.5/2) +
    geom_point(aes(x = 'Male', y = weighted_mean_male),
               colour = 'maroon', shape = 18, size = 3.5) +
    geom_point(aes(x = 'Female', y = weighted_mean_female),
               colour = 'maroon', shape = 18, size = 3.5) +
    scale_color_manual(values = c("Male" = "#5b9877",
                                  "Female" = "#e7a29c")) +
    theme(legend.position = 'none',
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_y_continuous(limits = c(0, 50),
                       breaks = seq(0, 50, by = 10)) +
    labs(x = 'Sex',
         y = 'Prevalence (%)\n')

est_by_sex

# <0.7
est_by_sex_0.7 <- estimate_0.7 |> 
  ggplot() +
  # geom_hline(yintercept = weighted_mean_male_0.7,
  #            colour = "#5b9877", linetype = 2) +
  # geom_hline(yintercept = weighted_mean_female_0.7,
  #            colour = "#e7a29c", linetype = 2) +
  # geom_hline(yintercept = median_male_0.7,
  #            colour = "#5b9877", linetype = 3) +
  # geom_hline(yintercept = median_female_0.7,
  #            colour = "#e7a29c", linetype = 3) +
  geom_beeswarm(aes(x = sex, y = estimate, colour = sex),
                alpha = 0.6, size = 3.1, cex = 2.3) +
  geom_boxplot(aes(x = sex, y = estimate),
               alpha = 0, width = 0.5/2) +
  geom_point(aes(x = 'Male', y = weighted_mean_male_0.7),
             colour = 'maroon', shape = 18, size = 3.5) +
  geom_point(aes(x = 'Female', y = weighted_mean_female_0.7),
             colour = 'maroon', shape = 18, size = 3.5) +
  scale_color_manual(values = c("Male" = "#5b9877",
                                "Female" = "#e7a29c")) +
  theme(legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(0, 50),
                     breaks = seq(0, 50, by = 10)) +
  labs(x = 'Sex',
       y = 'Prevalence (%)\n')

est_by_sex_0.7

# In this plot, deficient/intermediate males?

# By island
est_by_sex_by_island <- est_by_sex + facet_wrap(~island) 
est_by_sex_0.7_by_island <- est_by_sex_0.7 + facet_wrap(~island) # <0.7

# Deficient-only
est_by_sex_by_island <- estimate |> 
  ggplot() +
  geom_beeswarm(aes(x = sex, y = estimate, colour = sex),
                alpha = 0.6, size = 3.1, cex = 2.3) +
  scale_color_manual(values = c("Male" = "#5b9877",
                                "Female" = "#e7a29c")) +
  theme(legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(0, 50),
                     breaks = seq(0, 50, by = 10)) +
  labs(x = 'Sex',
       y = 'Prevalence (%)\n') +
  facet_wrap(~island)

est_by_sex_0.7_by_island <- estimate_0.7 |>
  ggplot() +
  geom_beeswarm(aes(x = sex, y = estimate, colour = sex),
                alpha = 0.6, size = 3.1, cex = 2.3) +
  scale_color_manual(values = c("Male" = "#5b9877",
                                "Female" = "#e7a29c")) +
  theme(legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(0, 50),
                     breaks = seq(0, 50, by = 10)) +
  labs(x = 'Sex',
       y = 'Prevalence (%)\n') +
  facet_wrap(~island)

est_by_sex_by_island
est_by_sex_0.7_by_island

## Plot: Sample sizes ------------------------------------------------------

n_by_sex <- sample |> 
  ggplot() +
    geom_histogram(aes(x = n, fill = sex, colour = sex),
                   alpha = 0.7, bins = 35, position = 'identity') +
    scale_colour_manual(values = c("Male" = "#accbb9",
                                   "Female" = "#e7a29c")) +
    scale_fill_manual(values = c("Male" = "#accbb9",
                                 "Female" = "#e7a29c")) +
    theme(legend.position = c(0.82, 0.80)) +
    scale_x_continuous(breaks = seq(0, 1000, by = 100)) +
    scale_y_continuous(limits = c(0, 11),
                       breaks = seq(0, 12, by = 2)) +
    labs(x = '\nSample size  ',
         y = 'Number of studies\n')

n_by_sex

## Plot: Male percentage ----------------------------------------------------
mf_percentage <- prevalence |>
  mutate(n_male = if_else(is.na(n_male), 0, n_male),
         n_female = if_else(is.na(n_female), 0, n_female),
         n_total = n_male + n_female,
         mf_perc = 100 * (n_male / n_total),
         n_total_cat = case_when(
           n_total <= 10 ~ '[1, 10)',
           n_total <= 100 ~ '[10, 100)',
           n_total <= 1000 ~ '[100, 1000)',
           is.na(n_total) ~ NA_character_,
           TRUE ~ 'Check me!'
         ) |> factor(levels = c('[1, 10)', '[10, 100)', '[100, 1000)'))) |> 
  dplyr::select(site_name, n_male, n_female, n_total_cat, mf_perc)

mf_perc_plot <- mf_percentage |> 
  mutate(site_name = factor(site_name),
         site_name = fct_reorder(site_name, mf_perc)) |> 
  ggplot(aes(x = mf_perc, y = site_name)) +
    geom_point(aes(size = n_total_cat),
               colour = "gray50", alpha = 0.7, shape = 18) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
    labs(x = "\nProportion of males (%)",
         y = "") +
    scale_size_discrete(drop = FALSE, range = c(1, 4)) +
    theme(legend.position = 'bottom')

mf_perc_plot

# End session
xfun::session_info()












