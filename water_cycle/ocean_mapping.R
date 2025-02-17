library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
world = ne_countries(scale = "medium", returnclass = "sf")
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

# ocean isotopes ----
#Schmidt, G.A., G. R. Bigg and E. J. Rohling. 1999. "Global Seawater Oxygen-18 Database - v1.22" 
#https://data.giss.nasa.gov/o18data/
ocean = read_csv("water_cycle/data/ocean_isotope.csv")
ocean = ocean[, 1:10]
ocean = ocean %>%
  mutate(dD = na_if(dD, "**")) %>%
  mutate(across(c(1:9), as.numeric))
m1 = lm(dD~d18O, data = ocean)
summary(m1)
ggplot(ocean, aes(x = d18O, y = dD, fill = abs(Latitude))) +
  geom_abline(slope = 8, intercept = 10) +
  geom_point(size = 3, shape = 21) +
  scale_fill_distiller(palette = "RdBu") +
  theme_bw() + theme +
  scale_x_continuous(limits = c(-10, 5)) +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"),
       fill = expression(paste("Latitude (", degree, "N/S)")))

ggplot(ocean, aes(x = Latitude, y = dD-8*d18O)) +
  geom_point(shape = 21, size =3) +
  theme_bw() + theme +
  labs(x = expression(paste("Latitude (", degree, "S-N)")),
       y = "d-excess (\u2030)")
ggplot(ocean, aes(x = Longitude, y = dD-8*d18O)) +
  geom_point(shape = 21, size =3) +
  theme_bw() + theme +
  labs(x = expression(paste("Longitude (", degree, "W-E)")),
       y = "d-excess (\u2030)")

ocean = ocean %>% drop_na(dD, d18O)
ggplot(data = world) +
  geom_sf(fill = "gray", color = NA) +
  coord_sf(expand = FALSE) +
  geom_point(data = ocean, aes(x = Longitude, y = Latitude, fill = dD-d18O*8), shape = 21, size = 2) +
  scale_fill_distiller(palette = "RdBu", limits = c(-12, 12)) +
  scale_x_continuous(breaks = seq(-180, 180, 30)) +
  scale_y_continuous(breaks = seq(-90, 90, 30)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray", linetype = "dotted"),
    panel.background = element_rect(fill = "white"),
    axis.text = element_text(size = 10),
    axis.title = element_blank()
  ) +
  labs(fill = "d-excess (\u2030)")


ocean.d18 = ocean %>% drop_na(d18O) %>% filter(d18O >=-9 & d18O <= 4)
ggplot(data = world) +
  geom_sf(fill = "gray", color = NA) +
  coord_sf(expand = FALSE) +
  geom_point(data = ocean.d18, aes(x = Longitude, y = Latitude, fill = d18O), shape = 21, size = 2) +
  scale_fill_distiller(palette = "RdBu", limits = c(-9, 4)) +
  scale_x_continuous(limits = c(-10,50), breaks = seq(0, 60, 10)) +
  scale_y_continuous(limits = c(30,70), breaks = seq(30, 60, 10)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray", linetype = "dotted"),
    panel.background = element_rect(fill = "white"),
    axis.text = element_text(size = 10),
    axis.title = element_blank()
  ) +
  labs(fill = expression(delta^"18"*"O (\u2030)"))

## Closed equilibrium system (ocean isotopes) ----
ocean = read_csv("data/ocean_isotope.csv")
ocean = ocean[, 1:10]
ocean = ocean %>%
  mutate(dD = na_if(dD, "**")) %>%
  mutate(across(c(1:9), as.numeric))
ctrl = function() {
  vars = list(
    d18 = 0,
    d2 = 0,
    f = 0.9,
    temp = 10,
    temp.a = 0
  )
}
ces = function(vars){
  list2env(vars, environment())
  R18_VSMOW <- 0.0020052 # ratio
  R2_VSMOW <- 0.00015576 # ratio
  R18 = (d18 * 1e-3 + 1) * R18_VSMOW
  R2 = (d2 * 1e-3 + 1) * R2_VSMOW
  alpha18_l_v = exp(-2.0667 * 10^-3 - 0.4156/(temp+273.15) + (1.137*10^3)/(temp+273.15)^2) # Majoube (1971)
  alpha2_l_v = exp(52.612*10^-3 - 76.248/(temp+273.15) + (24.844*10^3)/(temp+273.15)^2) # Majoube (1971)
  R18_l = R18 / (1/alpha18_l_v * (1-f) + f)
  R2_l = R2 / (1/alpha2_l_v * (1-f) + f)
  R18_v = R18_l / alpha18_l_v
  R2_v = R2_l / alpha2_l_v
  d18_l = (R18_l / R18_VSMOW - 1) * 1e3
  d2_l = (R2_l / R2_VSMOW - 1) * 1e3
  d18_v = (R18_v / R18_VSMOW - 1) * 1e3
  d2_v = (R2_v / R2_VSMOW - 1) * 1e3
  # condensation 
  alpha18_l_v = exp(-2.0667 * 10^-3 - 0.4156/(temp.a+273.15) + (1.137*10^3)/(temp.a+273.15)^2) # Majoube (1971)
  alpha2_l_v = exp(52.612*10^-3 - 76.248/(temp.a+273.15) + (24.844*10^3)/(temp.a+273.15)^2) # Majoube (1971)
  d18_c = (d18_v + 1e3) * alpha18_l_v - 1e3
  d2_c = (d2_v + 1e3) * alpha2_l_v - 1e3
  results = data.frame("temp" = rep(temp),
                       "f" = rep(f),
                       "d18l" = rep(d18_l), "d2l" = rep(d2_l),
                       "d18v" = rep(d18_v), "d2v" = rep(d2_v),
                       "d18c" = rep(d18_c), "d2c" = rep(d2_c))
}

vars = ctrl()
vars$f = seq(0.1, 0.9, 0.2)
sens.f = ces(vars)
m1 = lm(data = ocean, dD ~ d18O)
summary(m1) 
m2 = lm(data = sens.f, d2l ~ d18l)
summary(m2)
m3 = lm(data = sens.f, d2c ~ d18c)
summary(m3)
ggplot(sens.f) +
  # geom_abline(slope = 8, intercept = 10) +
  geom_point(data = ocean, aes(x = d18O, y = dD), shape = 21, size = 2, color = "coral") +
  geom_abline(slope = 7.23, intercept = -1.54,  color = "firebrick1", linetype = "dashed") +
  annotate(geom = "point", x = 0, y = 0, shape = 23, size = 3, fill = "green") +
  geom_point(aes(x = d18l, y = d2l, fill = f), shape = 21, size = 3) +
  geom_point(aes(x = d18v, y = d2v, fill = f), shape = 22, size = 3) +
  geom_point(aes(x = d18c, y = d2c, fill = f), shape = 23, size = 3) +
  # geom_line(size = 1, linetype = "dashed", color = "red") +
  # geom_line(aes(x = d18v, y = d2v),
  #           size = 1, linetype = "dashed", color = "blue") +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"), 
       fill = expression(italic(f)))

vars = ctrl()
vars$temp = seq(0, 30, 5)
sens.t = ces(vars)
ggplot(sens.t) +
  # geom_abline(slope = 8, intercept = 10) +
  geom_point(data = ocean, aes(x = d18O, y = dD), size = 2, color = "coral") +
  geom_abline(slope = 7.23, intercept = -1.54,  color = "firebrick1", linetype = "dashed") +
  annotate(geom = "point", x = 0, y = 0, shape = 23, size = 3, fill = "green") +
  geom_point(aes(x = d18l, y = d2l, fill = temp), shape = 21, size = 3) +
  geom_point(aes(x = d18v, y = d2v, fill = temp), shape = 22, size = 3) +
  geom_point(aes(x = d18c, y = d2c, fill = temp), shape = 23, size = 3) +
  # geom_line(size = 1, linetype = "dashed", color = "red") +
  # geom_line(aes(x = d18v, y = d2v),
  #           size = 1, linetype = "dashed", color = "blue") +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"), 
       fill = expression(paste("T (", degree, "C)")))


## Mediterranean ----
Med = ocean %>% 
  filter(Longitude > -5 & Longitude < 40) %>%
  filter(Latitude > 30 & Latitude < 50) %>%
  drop_na(d18O, dD) %>%
  mutate(d.excess = dD - d18O * 8)
ggplot(Med, aes(x = d18O, y = d.excess)) +
  geom_point(shape = 21, size = 3) +
  scale_fill_distiller(palette = "RdBu") +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression("d-excess (\u2030)"))

## Baltic ----
Baltic = ocean %>% 
  filter(Longitude > 0 & Longitude < 30) %>%
  filter(Latitude > 50 & Latitude < 60) %>%
  drop_na(d18O, dD) %>%
  mutate(d.excess = dD - d18O * 8)
ggplot(Baltic, aes(x = d18O, y = d.excess)) +
  geom_point(shape = 21) +
  scale_fill_distiller(palette = "RdBu")
