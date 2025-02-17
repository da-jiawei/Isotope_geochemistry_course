library(tidyverse)
# library(sf)
# library(rnaturalearth)
# library(rnaturalearthdata)
# world = ne_countries(scale = "medium", returnclass = "sf")
theme = theme(axis.text.x = element_text(margin = margin(t = 0.1, unit = "cm")),
              axis.text.y = element_text(margin = margin(r = 0.1, unit = "cm")),
              axis.ticks.length=unit(0.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              text = element_text(color = "black", size = 12),
              axis.title = element_text(size = 15), 
              axis.text = element_text(color = "black", size = 12),
              plot.title = element_text(hjust = 0.1, vjust = -10),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
source("water_cycle/constants_coefficients.R")

ocean = read_csv("water_cycle/data/ocean_isotope.csv")
ocean = ocean[, 1:10]
ocean = ocean %>%
  mutate(dD = na_if(dD, "**")) %>%
  mutate(across(c(1:9), as.numeric))

Med = ocean %>% 
  filter(Longitude > -5 & Longitude < 40) %>%
  filter(Latitude > 30 & Latitude < 50) %>%
  drop_na(d18O, dD) %>%
  mutate(d.excess = dD - d18O * 8)
ggplot(Med, aes(x = d18O, y = d.excess, fill = Longitude)) +
  geom_point(shape = 21, size = 3) +
  scale_fill_viridis_c() +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression("d-excess (\u2030)"))


# inflow - evaporation ----
# dV_l / dt = F_in - E
# d(delta_l * V_l) / dt = F_in * delta_in - E * delta_E
# numerical model 
Medi_numerical = function(temp, w, h, n, duration){
  alpha_eq(temp)
  V_med = 3750000 # km3
  ocean_density = 1.03 # g/cm3
  M_med = 18 * ocean_density * 1e12 * V_med # mol
  F_in = 0.8 # sv: 1e6m3/s
  F_in = F_in * 1e6 * ocean_density * 1e3 * 18 * 3600 * 24 * 365 # mol/yr
  E = F_in / n
  d18_in = 0
  d2_in = 0
  diff18 = diffratio_18 * (1-w) + w
  diff2 = diffratio_2 * (1-w) + w
  dt = 1
  time = seq(0, duration, dt)
  results = data.frame("time" = time)
  results$d18m[1] = 1 # initial ocean water
  results$d2m[1] = 8 # initial ocean water
  results$M_med[1] = M_med
  for (i in 1:(length(time)-1)) {
    d_M = (F_in - E) * dt
    results$d18v[i] = (results$d18m[i] + 1000) / (alpha18_l_v * diff18 * (1 - h) + alpha18_l_v * h) - 1000
    results$d2v[i] = (results$d2m[i] + 1000) / (alpha2_l_v * diff2 * (1 - h) + alpha2_l_v * h) - 1000
    D_d18m = ((d18_in - results$d18m[i]) * (F_in / results$M_med[i]) - (results$d18v[i] - results$d18m[i]) * (E / results$M_med[i])) * dt
    D_d2m = ((d2_in - results$d2m[i]) * (F_in / results$M_med[i]) - (results$d2v[i] - results$d2m[i]) * (E / results$M_med[i])) * dt
    results$d18m[i+1] = results$d18m[i] + D_d18m
    results$d2m[i+1] = results$d2m[i] + D_d2m
    # results$d2m[i+1] = ((d2_in * F_in - results$d2v[i] * E) * dt - results$d2m[i] * d_M) / results$M_med[i]
    results$M_med[i+1] = results$M_med[i] + d_M
  }
  results$d18v[i+1] = (results$d18m[i+1] + 1000) / (alpha18_l_v * diff18 * (1 - h) + alpha18_l_v * h) - 1000
  results$d2v[i+1] = (results$d2m[i+1] + 1000) / (alpha2_l_v * diff2 * (1 - h) + alpha2_l_v * h) - 1000
  results$d.excess = results$d2m - 8*results$d18m
  ggplot(results, aes(x = d18m, y = d.excess)) +
    geom_point(data = Med, aes(x = d18O, y = d.excess), size = 2, shape = 21) +
    geom_line(size = 1) +
    theme_bw() +
    labs(x = expression(delta^"18"*"O (\u2030)"), y = expression("d-excess (\u2030)"))
  # ggplot(results) +
  #   geom_line(aes(x = time, y = d18m), color = "firebrick1") +
  #   geom_line(aes(x = time, y = d18v), color = "dodgerblue2") +
  #   theme_bw()
}
Medi_numerical(20, 0, 0.4, 11, 1e4)
