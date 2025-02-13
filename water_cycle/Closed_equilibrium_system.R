library(tidyverse)
library(ggpubr)
source('constants_coefficients.R')
# Closed equilibrium system ----
ctrl = function() {
  vars = list(
    d18 = 0,
    d2 = 0,
    f = 0.9,
    temp = 10
  )
}

ces = function(vars){
  list2env(vars, environment())
  R18 = (d18 * 1e-3 + 1) * R18_VSMOW
  R2 = (d2 * 1e-3 + 1) * R2_VSMOW
  alpha_eq(temp) # evaporation of ocean water
  R18_l = R18 / (1/alpha18_l_v * (1-f) + f)
  R2_l = R2 / (1/alpha2_l_v * (1-f) + f)
  R18_v = R18_l / alpha18_l_v
  R2_v = R2_l / alpha2_l_v
  d18_l = (R18_l / R18_VSMOW - 1) * 1e3
  d2_l = (R2_l / R2_VSMOW - 1) * 1e3
  d18_v = (R18_v / R18_VSMOW - 1) * 1e3
  d2_v = (R2_v / R2_VSMOW - 1) * 1e3
  results = data.frame("temp" = rep(temp), "f" = rep(f),
                       "d18l" = rep(d18_l), "d2l" = rep(d2_l),
                       "d18v" = rep(d18_v), "d2v" = rep(d2_v))
}

# sensitivity test ---
# fraction
vars = ctrl()
vars$f = seq(0,1,0.1)
sens.f = ces(vars)
p1 = ggplot(sens.f) +
  geom_line(aes(x = f, y = d18l), color = "dodgerblue") +
  geom_line(aes(x = f, y = d18v), color = "firebrick1") +
  annotate("text", x = 0.9, y = -8, label = expression(delta^"18"*"O"[vapor])) +
  annotate("text", x = 0.9, y = 3, label = expression(delta^"18"*"O"[liquid])) +
  theme_bw() +
  labs(x = expression(italic(f)), 
       y = expression(delta^"18"*"O (\u2030)"))

p2 = ggplot(sens.f) +
  geom_point(aes(x = d18l, y = d2l, fill = f), shape = 21, size = 3) +
  geom_point(aes(x = d18v, y = d2v, fill = f), shape = 22, size = 3) +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme_bw() +
  labs(x = expression(delta^"18"*"O (\u2030)"), 
       y = expression(delta*"D (\u2030)"),
       fill = expression(italic(f)))

ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv")

# temperature
for (i in 1:4) {
  vars$temp = (i-1)*10
  if (i == 1) {
    sens = ces(vars)
  } else {
    results = ces(vars)
    sens = rbind(sens, results)
  }
}

p1 = ggplot(sens, aes(color = temp)) +
  geom_path(aes(x = f, y = d18l)) +
  geom_path(aes(x = f, y = d18v)) +
  scale_color_distiller(palette = "RdBu") +
  annotate("text", x = 0.9, y = -6, label = expression(delta^"18"*"O"[vapor])) +
  annotate("text", x = 0.9, y = 4, label = expression(delta^"18"*"O"[liquid])) +
  theme_bw() +
  labs(x = expression(italic(f)), 
       y = expression(delta^"18"*"O (\u2030)"),
       color = expression(paste("T (", degree, "C)")))

p2 = ggplot(sens) +
  geom_path(aes(x = d18l, y = d2l, color = temp)) +
  geom_path(aes(x = d18v, y = d2v, color = temp)) +
  scale_color_distiller(palette = "RdBu") +
  theme_bw() +
  labs(x = expression(delta^"18"*"O (\u2030)"), 
       y = expression(delta*"D (\u2030)"),
       color = expression(paste("T (", degree, "C)")))

ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv", common.legend = TRUE)
