# =========================================================
# MAP OF SAMPLE SITES
# Europe-wide map + map with regional insets
# =========================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)# =========================================================
  # MAP OF SAMPLE SITES
  # Europe-wide map + map with regional insets
  # =========================================================
  
  suppressPackageStartupMessages({
    library(sf)
    library(dplyr)
    library(tibble)
    library(ggplot2)
    library(ggrepel)
    library(rnaturalearth)
    library(rnaturalearthdata)
    library(cowplot)
  })
  
  # =========================================================
  # 0. SETTINGS
  # =========================================================
  
  out_dir <- "site_map_outputs"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  forest_cols <- c(
    "Sonian"   = "#B8860B",
    "Rajhenav" = "#2E7D32",
    "Strødam"  = "#466C95",
    "Suserup"  = "#5FA8D3"
  )
  
  # =========================================================
  # 1. SITE COORDINATES
  # =========================================================
  # Coordinates converted from degrees-minutes-seconds to decimal degrees.
  # Latitude north and longitude east are positive.
  
  sites <- tribble(
    ~forest,    ~country,    ~lat,        ~lon,
    "Sonian",   "Belgium",   50.7517778,   4.4224722,
    "Rajhenav", "Slovenia",  45.67,  15.011,
    "Strødam",  "Denmark",   55.9577500,  12.2720000,
    "Suserup",  "Denmark",   55.3786111,  11.5605833
  ) %>%
    mutate(
      forest = factor(forest, levels = names(forest_cols)),
      label = as.character(forest)
    )
  
  sites_sf <- st_as_sf(
    sites,
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  )
  
  # ---------------------------------------------------------
  # Additional NAT-MAN sites not used in the present study
  # ---------------------------------------------------------
  
  other_natman <- tribble(
    ~site_name,                            ~country_code, ~lat,   ~lon,
    "Silkeborg Vesterskov, Knagerne",      "DK",          56.13,   9.53,
    "Møns Klinteskov, Kalsterbjerg",       "DK",          54.96,  12.54,
    "Velling Skov",                        "DK",          56.04,   9.50,
    "Kekes",                               "H",           47.80,  19.85,
    "Õserdõ",                              "H",           48.05,  20.43,
    "Krokar",                              "SI",          45.54,  14.78,
    "Utrecht, Amelisweerd",                "NL",          52.10,   5.18,
    "Veluwe, Dassenberg",                  "NL",          52.07,   5.88,
    "Veluwe, Drie",                        "NL",          52.07,   5.88,
    "Veluwe, Gortelse Bos",                "NL",          52.07,   5.88,
    "Utrecht, Oostbroek",                  "NL",          52.10,   5.18,
    "Veluwe, Speulderbos",                 "NL",          52.25,   5.72,
    "Veluwe, Weversbergen",                "NL",          52.07,   5.88,
    "Utrecht, Wulperhorst",                "NL",          52.10,   5.18
  )
  
  other_natman_sf <- st_as_sf(
    other_natman,
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  )
  
  # =========================================================
  # 2. BASE MAP DATA
  # =========================================================
  
  world <- ne_countries(
    scale = "medium",
    returnclass = "sf"
  )
  
  europe <- world %>%
    filter(
      continent == "Europe" |
        admin %in% c("Turkey", "Cyprus")
    )
  
  # =========================================================
  # 3. COMMON THEME
  # =========================================================
  
  theme_map <- function(base_size = 9) {
    theme_void(base_size = base_size) +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(
          face = "bold",
          hjust = 0,
          size = base_size + 2
        ),
        plot.subtitle = element_text(
          hjust = 0,
          size = base_size,
          colour = "grey30"
        ),
        plot.margin = margin(5, 5, 5, 5)
      )
  }
  # Manual label positions
  # Adjust label_lon and label_lat until the boxes sit where you want them.
  site_label_pos <- sites %>%
    mutate(
      label_lon = case_when(
        forest == "Sonian"   ~ 1.2,
        forest == "Rajhenav" ~ 13.95,
        forest == "Strødam"  ~ 14.5,
        forest == "Suserup"  ~ 9.00,
        TRUE ~ lon
      ),
      label_lat = case_when(
        forest == "Sonian"   ~ 50.75,
        forest == "Rajhenav" ~ 47.2,
        forest == "Strødam"  ~ 57.0,
        forest == "Suserup"  ~ 54.15,
        TRUE ~ lat
      )
    )
  # =========================================================
  # 4. EUROPE-WIDE MAP
  # =========================================================
  site_label_segments <- site_label_pos %>%
    mutate(
      x = lon,
      y = lat,
      xend = label_lon,
      yend = label_lat
    )
  
  p_europe_sites <- ggplot() +
    geom_sf(
      data = europe,
      fill = "grey92",
      colour = "grey65",
      linewidth = 0.25
    ) +
    geom_sf(
      data = other_natman_sf,
      shape = 21,
      size = 1,
      stroke = 0.30,
      fill = "black",
      colour = "grey80", alpha=0.75
    ) +
    geom_sf(
      data = sites_sf,
      aes(fill = forest),
      shape = 21,
      size = 1.5,
      stroke = 0.3,
      colour = "black"
    ) +
    geom_segment(
      data = site_label_segments,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      linewidth = 0.25,
      colour = "grey30"
    ) +
    geom_label(
      data = site_label_pos,
      aes(
        x = label_lon,
        y = label_lat,
        label = label,
        fill = forest
      ),
      colour = "black",
      size = 2.5,
      alpha = 0.99,
      label.padding = unit(0.12, "lines"),
      label.r = unit(0.12, "lines"),
      label.size = 0.25,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = forest_cols) +
    scale_colour_manual(values = forest_cols) +
    coord_sf(
      xlim = c(-5, 21),
      ylim = c(44, 58),
      expand = FALSE
    ) +
    labs(
      title = "European deadwood fungal survey sites",
      subtitle = "Forest reserves included in the historical and contemporary fruitbody inventories"
    ) +
    theme_map(base_size = 7.5)
  
  print(p_europe_sites)
  
  ggsave(
    filename = file.path(out_dir, "map_sample_sites_europe.png"),
    plot = p_europe_sites,
    width = 7.2/2,
    height = 5.2/2,
    dpi = 900,
    bg = "white"
  )
  
  # =========================================================
  # 5. INSET MAP FUNCTIONS
  # =========================================================
  
  make_zoom_map <- function(xlim, ylim, title_text, label_sites = TRUE) {
    p <- ggplot() +
      geom_sf(
        data = europe,
        fill = "grey94",
        colour = "grey65",
        linewidth = 0.22
      ) +
      geom_sf(
        data = other_natman_sf,
        shape = 21,
        size = 2.0,
        stroke = 0.25,
        fill = "black",
        colour = "black"
      ) +
      geom_sf(
        data = sites_sf,
        aes(fill = forest),
        shape = 21,
        size = 2.7,
        stroke = 0.35,
        colour = "black"
      ) +
      scale_fill_manual(values = forest_cols, drop = FALSE) +
      coord_sf(
        xlim = xlim,
        ylim = ylim,
        expand = FALSE
      ) +
      labs(title = title_text) +
      theme_void(base_size = 8) +
      theme(
        legend.position = "none",
        plot.title = element_text(
          face = "bold",
          size = 8.5,
          hjust = 0.02,
          margin = margin(0, 0, 2, 0)
        ),
        panel.border = element_rect(
          colour = "grey25",
          fill = NA,
          linewidth = 0.35
        ),
        plot.background = element_rect(
          fill = "white",
          colour = "grey25",
          linewidth = 0.35
        ),
        plot.margin = margin(3, 3, 3, 3)
      )
    
    if (label_sites) {
      p <- p +
        geom_text_repel(
          data = sites %>%
            filter(
              lon >= xlim[1],
              lon <= xlim[2],
              lat >= ylim[1],
              lat <= ylim[2]
            ),
          aes(
            x = lon,
            y = lat,
            label = label,
            colour = forest
          ),
          size = 2.6,
          min.segment.length = 0,
          segment.linewidth = 0.2,
          segment.colour = "grey35",
          box.padding = 0.25,
          point.padding = 0.2,
          seed = 42,
          show.legend = FALSE
        ) +
        scale_colour_manual(values = forest_cols, drop = FALSE)
    }
    
    p
  }
  
  p_belgium_inset <- make_zoom_map(
    xlim = c(2.4, 6.5),
    ylim = c(49.4, 51.6),
    title_text = "Belgium"
  )
  
  p_denmark_inset <- make_zoom_map(
    xlim = c(9.8, 13.4),
    ylim = c(54.7, 56.4),
    title_text = "Denmark"
  )
  
  p_slovenia_inset <- make_zoom_map(
    xlim = c(13.0, 16.8),
    ylim = c(44.8, 46.8),
    title_text = "Slovenia"
  )
  
  # =========================================================
  # 6. MAIN MAP PREPARED FOR INSETS
  # =========================================================
  
  p_main_for_insets <- ggplot() +
    geom_sf(
      data = europe,
      fill = "grey92",
      colour = "grey65",
      linewidth = 0.25
    ) +
    geom_sf(
      data = other_natman_sf,
      shape = 21,
      size = 2.2,
      stroke = 0.30,
      fill = "black",
      colour = "black"
    ) +
    geom_sf(
      data = sites_sf,
      aes(fill = forest),
      shape = 21,
      size = 3.0,
      stroke = 0.45,
      colour = "black"
    ) +
    scale_fill_manual(values = forest_cols) +
    coord_sf(
      xlim = c(-8, 28),
      ylim = c(41, 59),
      expand = FALSE
    ) +
    labs(
      title = "European deadwood fungal survey sites",
      subtitle = "Inset panels show the local position of the four forest reserves"
    ) +
    theme_map(base_size = 9)
  
  # =========================================================
  # 7. MAP WITH INSETS
  # =========================================================
  # The x, y, width and height arguments below are relative to the whole figure.
  # Adjust them slightly if an inset overlaps an important map feature.
  
  p_europe_sites_insets <- ggdraw(p_main_for_insets) +
    draw_plot(p_belgium_inset,  x = 0.05, y = 0.07, width = 0.265, height = 0.255) +
    draw_plot(p_denmark_inset,  x = 0.33, y = 0.07, width = 0.265, height = 0.255) +
    draw_plot(p_slovenia_inset, x = 0.675, y = 0.07, width = 0.265, height = 0.255)
  
  print(p_europe_sites_insets)
  
  ggsave(
    filename = file.path(out_dir, "map_sample_sites_europe_with_insets.png"),
    plot = p_europe_sites_insets,
    width = 7.2,
    height = 5.2,
    dpi = 600,
    bg = "white"
  )
  
  # =========================================================
  # 8. OPTIONAL: PRINT COORDINATE TABLE
  # =========================================================
  
  print(
    sites %>%
      select(forest, country, latitude_decimal = lat, longitude_decimal = lon)
  )
  
  print(
    other_natman %>%
      select(site_name, country_code, latitude_decimal = lat, longitude_decimal = lon)
  )
  
  
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(cowplot)
})

# =========================================================
# 0. SETTINGS
# =========================================================

out_dir <- "site_map_outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

forest_cols <- c(
  "Sonian"   = "#B8860B",
  "Rajhenav" = "#2E7D32",
  "Strødam"  = "#466C95",
  "Suserup"  = "#5FA8D3"
)

# =========================================================
# 1. SITE COORDINATES
# =========================================================
# Coordinates converted from degrees-minutes-seconds to decimal degrees.
# Latitude north and longitude east are positive.

sites <- tribble(
  ~forest,    ~country,    ~lat,        ~lon,
  "Sonian",   "Belgium",   50.7517778,   4.4224722,
  "Rajhenav", "Slovenia",  45.6565833,  15.0093611,
  "Strødam",  "Denmark",   55.9577500,  12.2720000,
  "Suserup",  "Denmark",   55.3786111,  11.5605833
) %>%
  mutate(
    forest = factor(forest, levels = names(forest_cols)),
    label = as.character(forest)
  )

sites_sf <- st_as_sf(
  sites,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

# =========================================================
# 2. BASE MAP DATA
# =========================================================

world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

europe <- world %>%
  filter(
    continent == "Europe" |
      admin %in% c("Turkey", "Cyprus")
  )

# =========================================================
# 3. COMMON THEME
# =========================================================

theme_map <- function(base_size = 9) {
  theme_void(base_size = base_size) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(
        face = "bold",
        hjust = 0,
        size = base_size + 2
      ),
      plot.subtitle = element_text(
        hjust = 0,
        size = base_size,
        colour = "grey30"
      ),
      plot.margin = margin(5, 5, 5, 5)
    )
}

# =========================================================
# 4. EUROPE-WIDE MAP
# =========================================================

p_europe_sites <- ggplot() +
  geom_sf(
    data = europe,
    fill = "grey92",
    colour = "grey65",
    linewidth = 0.25
  ) +
  geom_sf(
    data = sites_sf,
    aes(fill = forest),
    shape = 21,
    size = 2.5,
    stroke = 0.45,
    colour = "black"
  ) +
  geom_label_repel(
    data = sites,
    aes(
      x = lon,
      y = lat,
      label = label,
      fill = forest
    ),
    colour = "black",
    size = 2.5, alpha= 0.8,
    min.segment.length = 0,
    segment.linewidth = 0.25,
    segment.colour = "grey30",
    box.padding = 0.35,
    point.padding = 0.25,
    label.padding = unit(0.12, "lines"),
    label.r = unit(0.12, "lines"),
    label.size = 0.25,
    seed = 42,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = forest_cols) +
  scale_colour_manual(values = forest_cols) +
  coord_sf(
    xlim = c(-5, 18),
    ylim = c(44, 58),
    expand = FALSE
  ) +
  labs(
    title = "European deadwood fungal survey sites",
    subtitle = "Forest reserves included in the historical and contemporary fruitbody inventories"
  ) +
  theme_map(base_size = 7.5)

print(p_europe_sites)

ggsave(
  filename = file.path(out_dir, "map_sample_sites_europe.png"),
  plot = p_europe_sites,
  width = 7.2/2,
  height = 5.2/2,
  dpi = 600,
  bg = "white"
)

# =========================================================
# 5. INSET MAP FUNCTIONS
# =========================================================

make_zoom_map <- function(xlim, ylim, title_text, label_sites = TRUE) {
  p <- ggplot() +
    geom_sf(
      data = europe,
      fill = "grey94",
      colour = "grey65",
      linewidth = 0.22
    ) +
    geom_sf(
      data = sites_sf,
      aes(fill = forest),
      shape = 21,
      size = 2.7,
      stroke = 0.35,
      colour = "black"
    ) +
    scale_fill_manual(values = forest_cols, drop = FALSE) +
    coord_sf(
      xlim = xlim,
      ylim = ylim,
      expand = FALSE
    ) +
    labs(title = title_text) +
    theme_void(base_size = 8) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        face = "bold",
        size = 8.5,
        hjust = 0.02,
        margin = margin(0, 0, 2, 0)
      ),
      panel.border = element_rect(
        colour = "grey25",
        fill = NA,
        linewidth = 0.35
      ),
      plot.background = element_rect(
        fill = "white",
        colour = "grey25",
        linewidth = 0.35
      ),
      plot.margin = margin(3, 3, 3, 3)
    )
  
  if (label_sites) {
    p <- p +
      geom_text_repel(
        data = sites %>%
          filter(
            lon >= xlim[1],
            lon <= xlim[2],
            lat >= ylim[1],
            lat <= ylim[2]
          ),
        aes(
          x = lon,
          y = lat,
          label = label,
          colour = forest
        ),
        size = 2.6,
        min.segment.length = 0,
        segment.linewidth = 0.2,
        segment.colour = "grey35",
        box.padding = 0.25,
        point.padding = 0.2,
        seed = 42,
        show.legend = FALSE
      ) +
      scale_colour_manual(values = forest_cols, drop = FALSE)
  }
  
  p
}

p_belgium_inset <- make_zoom_map(
  xlim = c(2.4, 6.5),
  ylim = c(49.4, 51.6),
  title_text = "Belgium"
)

p_denmark_inset <- make_zoom_map(
  xlim = c(9.8, 13.4),
  ylim = c(54.7, 56.4),
  title_text = "Denmark"
)

p_slovenia_inset <- make_zoom_map(
  xlim = c(13.0, 16.8),
  ylim = c(44.8, 46.8),
  title_text = "Slovenia"
)

# =========================================================
# 6. MAIN MAP PREPARED FOR INSETS
# =========================================================

p_main_for_insets <- ggplot() +
  geom_sf(
    data = europe,
    fill = "grey92",
    colour = "grey65",
    linewidth = 0.25
  ) +
  geom_sf(
    data = sites_sf,
    aes(fill = forest),
    shape = 21,
    size = 3.0,
    stroke = 0.45,
    colour = "black"
  ) +
  scale_fill_manual(values = forest_cols) +
  coord_sf(
    xlim = c(-8, 28),
    ylim = c(41, 59),
    expand = FALSE
  ) +
  labs(
    title = "European deadwood fungal survey sites",
    subtitle = "Inset panels show the local position of the four forest reserves"
  ) +
  theme_map(base_size = 9)

# =========================================================
# 7. MAP WITH INSETS
# =========================================================
# The x, y, width and height arguments below are relative to the whole figure.
# Adjust them slightly if an inset overlaps an important map feature.

p_europe_sites_insets <- ggdraw(p_main_for_insets) +
  draw_plot(p_belgium_inset,  x = 0.05, y = 0.07, width = 0.265, height = 0.255) +
  draw_plot(p_denmark_inset,  x = 0.33, y = 0.07, width = 0.265, height = 0.255) +
  draw_plot(p_slovenia_inset, x = 0.675, y = 0.07, width = 0.265, height = 0.255)

print(p_europe_sites_insets)

ggsave(
  filename = file.path(out_dir, "map_sample_sites_europe_with_insets.png"),
  plot = p_europe_sites_insets,
  width = 7.2,
  height = 5.2,
  dpi = 600,
  bg = "white"
)

# =========================================================
# 8. OPTIONAL: PRINT COORDINATE TABLE
# =========================================================

print(
  sites %>%
    select(forest, country, latitude_decimal = lat, longitude_decimal = lon)
)

