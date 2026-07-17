suppressMessages({library(macrosheds); library(dplyr); library(tidyr)})
root <- "C:/Users/nic/cz_disturbance"
sites <- c("w1","GSWS10","WS-3","WS-5","MC06")

rd <- if (requireNamespace("arrow", quietly=TRUE)) arrow::read_feather else feather::read_feather
ws <- as.data.frame(rd(file.path(root, "ms_data/v2/watershed_summaries.feather")))

# lat/long/area from bundled site data
sd <- ms_site_data %>% filter(site_code %in% sites) %>%
  select(network, domain, site_code, latitude, longitude, ws_area_ha)

# dominant NLCD land-cover class per site
nlcd_names <- c(
  nlcd_forest_dec="Deciduous forest", nlcd_forest_evr="Evergreen forest",
  nlcd_forest_mix="Mixed forest", nlcd_shrub="Shrub/scrub", nlcd_shrub_dwr="Dwarf scrub",
  nlcd_grass="Grassland", nlcd_sedge="Sedge/herbaceous", nlcd_pasture="Pasture/hay",
  nlcd_crop="Cultivated crops", nlcd_wetland_wood="Woody wetland",
  nlcd_wetland_herb="Herbaceous wetland", nlcd_barren="Barren", nlcd_water="Open water",
  nlcd_dev_open="Developed, open", nlcd_dev_low="Developed, low",
  nlcd_dev_med="Developed, medium", nlcd_dev_hi="Developed, high",
  nlcd_moss="Moss", nlcd_lichens="Lichens", nlcd_ice_snow="Ice/snow")
nlcd_cols <- intersect(names(nlcd_names), names(ws))

lulc <- ws %>% filter(site_code %in% sites) %>%
  select(site_code, all_of(nlcd_cols)) %>%
  pivot_longer(-site_code, names_to = "cls", values_to = "pct") %>%
  group_by(site_code) %>% slice_max(pct, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(site_code,
            dominant_lulc = nlcd_names[cls],
            dominant_lulc_pct = round(pct, 1))

# aspect degrees -> 8-point cardinal
card <- c("N","NE","E","SE","S","SW","W","NW")
terr <- ws %>% filter(site_code %in% sites) %>%
  transmute(site_code,
            elevation_m = round(elev_mean),
            slope_pct   = round(slope_mean, 1),
            aspect_deg  = round(aspect_mean),
            aspect_card = card[floor((aspect_mean %% 360) / 45 + 0.5) %% 8 + 1])

meta <- sd %>% left_join(terr, "site_code") %>% left_join(lulc, "site_code") %>%
  transmute(network, domain, site_code,
            latitude = round(latitude, 4), longitude = round(longitude, 4),
            area_ha = round(ws_area_ha, 1), elevation_m,
            slope_pct, aspect_deg, aspect_card,
            dominant_lulc, dominant_lulc_pct) %>%
  arrange(match(site_code, c("w1","GSWS10","WS-3","WS-5","MC06")))

print(as.data.frame(meta), width = 200)
write.csv(meta, file.path(root, "data/site_metadata.csv"), row.names = FALSE)
cat("\nwrote data/site_metadata.csv\n")
