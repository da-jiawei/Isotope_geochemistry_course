library(tidyverse)
library(readxl)
library(ggpubr)
source('constants_coefficients.R')
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

# read GNIP data ----
dat = read_xlsx("water_cycle/data/GNIP.xlsx")
dat2 = dat[, c(12:15, 35, 39)]
dat.d18 = dat2 %>%
  filter(`Measurand Symbol` == "O18") %>%
  group_by(`Sample Site Name`) %>%
  summarise(d18 = mean(`Measurand Amount`))
names(dat.d18) = c("site", "d18")
dat.d2 = dat2 %>%
  filter(`Measurand Symbol` == "H2") %>%
  group_by(`Sample Site Name`) %>%
  summarise(d18 = mean(`Measurand Amount`))
names(dat.d2) = c("site", "d2")
GNIP = merge(dat.d18, dat.d2, by = "site")

ggplot(GNIP, aes(x = d18, y = d2)) +
  geom_point(color = "gray") +
  # annotate("point", x = 0, y = 0, shape = 22, fill = "white", size = 4) +
  geom_abline(slope = 8, intercept = 10) +
  theme_bw() + theme +
  ggtitle("GMWL") +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"))

## Rayleigh distillation ----
# Rf = Ro * f ^ (alpha - 1) # alpha = Ref/R (effluent flux)
# set the initial vapor in equilibrium with ocean water
ctrl = function() {
  vars = list(
    d18o = -10,
    d2o = -10 * 8 + 10,
    f = 0.5,
    temp = 10
  )
}
ray = function(vars) {
  list2env(vars, environment())
  alpha_eq(temp)
  d18_f = (d18o + 1e3) * f ^ (alpha18_l_v - 1) - 1e3
  d2_f = (d2o + 1e3) * f ^ (alpha2_l_v - 1) - 1e3
  d18_p = (d18_f + 1e3) * alpha18_l_v - 1e3
  d2_p = (d2_f + 1e3) * alpha2_l_v -1e3
  results = data.frame("temp" = rep(temp), "f" = rep(f),
                       "d18f" = rep(d18_f), "d2f" = rep(d2_f),
                       "d18p" = rep(d18_p), "d2p" = rep(d2_p))
}

# fraction of remaining cloud
vars = ctrl()
vars$f = seq(0.01, 1, 0.01)
sens.f = ray(vars)
p1 = ggplot(sens.f) +
  geom_line(aes(x = f, y = d18p), color = "dodgerblue") +
  geom_line(aes(x = f, y = d18f), color = "firebrick1") +
  annotate("text", x = 0.8, y = 2, label = expression(delta^"18"*"O"[rain])) +
  annotate("text", x = 0.8, y = -9, label = expression(delta^"18"*"O"[vapor])) +
  theme_bw() +
  labs(x = expression(italic(f)), 
       y = expression(delta^"18"*"O (\u2030)")) +
  scale_x_continuous(limits = c(0,1))
vars$f = seq(0.1, 1, 0.1)
sens.f = ray(vars)
p2 = ggplot(sens.f) +
  geom_point(data = GNIP, aes(x = d18, y = d2), color = "gray") +
  geom_abline(intercept = 10, slope = 8) +
  geom_point(aes(x = d18p, y = d2p, fill = f), shape = 23, size = 3) +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme_bw() +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"),
       fill = expression(italic(f))) +
  guides(color = "none")
ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv")

# d-excess
ggplot(sens.f) +
  geom_point(data = GNIP, aes(x = d18, y = d2), color = "gray") +
  geom_abline(intercept = 10, slope = 8) +
  geom_point(aes(x = d18p, y = d2p, fill = f), shape = 23, size = 3) +
  geom_point(aes(x = d18f, y = d2f, fill = f), shape = 21, size = 3) +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  annotate(geom = "point", x = 0, y = 0, shape = 22, fill = "white", size = 4) +
  theme_bw() +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"),
       fill = expression(italic(f))) +
  guides(color = "none")
sens.f = sens.f %>%
  mutate(d_excess.p = d2p - 8*d18p,
         d_excess.f = d2f - 8*d18f)


# temperature
vars$f = seq(0.01, 1, 0.01)
vars$d18o = 0
vars$d2o = 0
for (i in 1:5) {
  vars$temp = (i-1) * 5
  sens.f = ray(vars)
  if(i == 1){
    sens = sens.f
  } else {
    sens = rbind(sens, sens.f)
  }
}

ggplot(sens, aes(x = d18p, y = d2p, color = temp)) +
  geom_point(data = GNIP, aes(x = d18, y = d2), color = "gray") +
  geom_abline(intercept = 10, slope = 8) +
  geom_path() +
  annotate(geom = "point", x = vars$d18, y = vars$d2, shape = 22, fill = "green", size = 4) +
  scale_color_distiller(palette = "RdBu") +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"),
       color = expression(paste("T (", degree, "C)"))) +
  scale_x_continuous(limits = c(-25, 5)) +
  scale_y_continuous(limits = c(-200, 10))

sens$d.excess = sens$d2p - 8*sens$d18p
ggplot(sens, aes(x = f, y = d.excess, group = temp, color = temp)) +
  geom_path() + 
  theme_bw() + theme +
  ggtitle("rain") +
  labs(x = expression(italic(f)[r]),
       y = "d-excess (\u2030)")

# Rayleigh distillation - numerical solution ----
temp.o = 20
temp.f = 0
dt = 0.1
length = (temp.o - temp.f) / dt + 1
results = data.frame(matrix(nrow = length, ncol = 6))
names(results) = c("temp", "conc", "d18f", "d2f", "d18p", "d2p")
for (i in 1:length) {
  if(i == 1) {
    results$temp[i] = temp.o
    Es = (1.0007 + 3.46e-8) * (611.21 * exp(17.502 * temp.o / (240.97 + temp.o)))
    results$conc[i] = Es / (8.3144 * (temp.o + 273.15))
    results$d18f[i] = -10
    results$d2f[i] = results$d18f[i] * 8 + 10
  } else {
    results$temp[i] = results$temp[i-1] + dt
    Es = (1.0007 + 3.46e-8) * (611.21 * exp(17.502 * results$temp[i] / (240.97 + results$temp[i])))
    results$conc[i] = Es / (8.3144 * (results$temp[i] + 273.15))
    f = results$conc[i] / results$conc[i-1]
    alpha_eq(results$temp[i])
    results$d18f[i] = (results$d18f[i-1] + 1000) * f ^ (alpha18_l_v - 1) - 1000
    results$d2f[i] = (results$d2f[i-1] + 1000) * f ^ (alpha2_l_v - 1) - 1000
    results$d18p[i] = results$d18f[i] / alpha18_l_v
    results$d2p[i] = results$d2f[i] / alpha2_l_v
  }
}

ggplot(results, aes(x = d18p, y = d2p)) +
  geom_point(data = GNIP, aes(x = d18, y = d2), color = "gray") +
  geom_abline(intercept = 10, slope = 8) +
  geom_point()



