# Faceted section plots (contour + points) of chlorophyll fluorescence
# by depth and latitude, one facet per 2016 sampling date, with
# interpolated regions masked below each profile's sampled depth.

library(dplyr)
library(tidyr)
library(ggplot2)
library(akima)   # irregular-grid interpolation for geom_contour_filled

darts <- read.csv("../../data/DaRTS_combined_data.csv")

# --- 1. Station -> latitude lookup --------------------------------------

station_lat <- data.frame(
  station = c(1, 2, 3, 4),
  lat     = c(43.903, 43.863, 43.811, 43.753)
)

darts2016 <- darts %>%
  filter(year == 2016) %>%
  left_join(station_lat, by = "station")

# --- 2. Bottom (max sampled) depth per station/date profile -------------
# This defines the mask boundary: anything interpolated deeper than a
# station's actual max sampled depth gets dropped before plotting.

bottom_depths <- darts2016 %>%
  filter(!is.na(fluor_mg_per_m3)) %>%
  group_by(date, station, lat) %>%
  summarize(max_depth = max(depth_m), .groups = "drop")

# --- 3. Interpolate one date onto a regular lat x depth grid, then mask -

interp_one_date <- function(df, bottom_df, n_lat = 100, n_depth = 100) {
  df <- df %>% filter(!is.na(fluor_mg_per_m3))
  if (nrow(df) < 4 || length(unique(df$lat)) < 2) return(NULL)
  
  interp_result <- akima::interp(
    x = df$lat, y = df$depth_m, z = df$fluor_mg_per_m3,
    xo = seq(min(df$lat), max(df$lat), length.out = n_lat),
    yo = seq(0, max(df$depth_m), length.out = n_depth),
    duplicate = "mean", linear = TRUE
  )
  
  grid_df <- expand.grid(lat = interp_result$x, depth_m = interp_result$y)
  grid_df$fluor <- as.vector(interp_result$z)
  
  # Bottom depth as a function of latitude, linearly interpolated between
  # the 4 station values, then used to mask the grid below each profile.
  bd <- bottom_df %>% arrange(lat)
  bottom_fun <- approxfun(bd$lat, bd$max_depth, rule = 2)
  grid_df$bottom <- bottom_fun(grid_df$lat)
  
  grid_df <- grid_df %>% filter(depth_m <= bottom)
  grid_df$date <- unique(df$date)
  grid_df
}

dates_2016 <- unique(darts2016$date)

grid_list <- lapply(dates_2016, function(d) {
  df_d <- darts2016 %>% filter(date == d)
  bd_d <- bottom_depths %>% filter(date == d)
  interp_one_date(df_d, bd_d)
})

grid_all <- bind_rows(grid_list)

# --- 4. Plot --------------------------------------------------------------

ggplot() +
  geom_contour_filled(
    data = grid_all,
    aes(x = lat, y = depth_m, z = fluor),
    alpha = 0.85
  ) +
  geom_point(
    data = darts2016 %>% filter(!is.na(fluor_mg_per_m3)),
    aes(x = lat, y = depth_m, color = fluor_mg_per_m3),
    size = 1.2, shape = 21, stroke = 0.3
  ) +
  scale_y_reverse() +
  scale_color_viridis_c(name = expression(Chl~fluor.~(mg~m^{-3}))) +
  facet_wrap(~ date) +
  labs(
    x = "Latitude (°N)",
    y = "Depth (m)",
    fill = expression(Chl~fluor.~(mg~m^{-3}))
  ) +
  theme_minimal()