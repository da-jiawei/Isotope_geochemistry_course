library(tidyverse)
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

# Evaporation from the ocean ----
ctrl = function() {
  vars = list(
    # SMOW
    d18l = 0,
    d2l = 0,
    RH = 0.8,
    temp = 20, # ocean temperature
    w = 0.7 # turbulence factor (1 - pure diffusion, 0 - pure turbulence)
  )
}

# assuming the isotopic compositions of the evaporation flux and atmospheric vapor are the same
cg = function(vars) {
  list2env(vars, environment())
  alpha_eq(temp)
  alpha18_diff = w * diffratio_18 + (1-w) 
  alpha2_diff = w * diffratio_2 + (1-w) 
  R18l = (d18l/1000 + 1) * R18_VSMOW
  R2l = (d2l/1000 + 1) * R2_VSMOW
  R18e = R18l / (alpha18_l_v * (alpha18_diff * (1 - RH)) + RH * alpha18_l_v)
  R2e = R2l / (alpha2_l_v * (alpha2_diff * (1 - RH)) + RH * alpha2_l_v)
  d18e = (R18e / R18_VSMOW - 1) * 1e3
  d2e = (R2e / R2_VSMOW - 1) * 1e3
  d_excess = d2e - 8*d18e
  results = data.frame("RH" = rep(RH), "temp" = rep(temp), "w" = rep(w),
                       "d18dv" = rep(d18e), "d2dv" = rep(d2e), "d_excess" = rep(d_excess))
}

# Relative humidity 
vars = ctrl()
vars$RH = seq(0.1, 1, 0.2)
sens.h = cg(vars)
p1 = ggplot(sens.h) +
  geom_point(data = GNIP, aes(x = d18, y = d2), color = "gray") +
  geom_abline(slope = 8, intercept = 10) +
  geom_point(aes(x = d18dv, y = d2dv, fill = RH*100), shape = 21, size = 3) +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  annotate("point", x = 0, y = 0, shape = 23, size = 4, fill = "white") +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"),
       fill = "RH (%)")
p2 = ggplot(sens.h, aes(x = RH*100, y = d_excess)) +
  geom_point(shape = 21, size = 3) +
  theme_bw() + theme +
  labs(x = "RH (%)", y = "d-excess (\u2030)")
ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv", common.legend = TRUE)

# temperature
vars = ctrl()
vars$temp = seq(0,30,5)
sens.t = cg(vars)
p1 = ggplot(sens.t) +
  geom_point(data = GNIP, aes(x = d18, y = d2), color = "gray") +
  geom_abline(slope = 8, intercept = 10) +
  geom_point(aes(x = d18dv, y = d2dv, fill = temp), shape = 21, size = 3) +
  scale_fill_distiller(palette = "RdBu") +
  annotate("point", x = 0, y = 0, shape = 23, size = 3) +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"),
       fill = expression(paste("T (", degree, "C)")))
p2 = ggplot(sens.t, aes(x = temp, y = d_excess)) +
  geom_point(shape = 21, size = 3) +
  theme_bw() + theme +
  labs(x = expression(paste("T (", degree, "C)")), y = "d-excess (\u2030)")
ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv", common.legend = TRUE)

# w - turbulence
vars = ctrl()
vars$w = seq(0.1, 1, 0.2)
sens.n = cg(vars)
p1 = ggplot(sens.n) +
  geom_point(data = GNIP, aes(x = d18, y = d2), color = "gray") +
  geom_abline(slope = 8, intercept = 10) +
  geom_point(aes(x = d18dv, y = d2dv, fill = w), shape = 21, size = 3) +
  scale_fill_distiller(palette = "RdBu") +
  annotate("point", x = 0, y = 0, shape = 23, size = 3) +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"),
       fill = expression(italic(w)))
p2 = ggplot(sens.n, aes(x = w, y = d_excess)) +
  geom_point(shape = 21, size = 3) +
  theme_bw() + theme +
  labs(x = expression(italic(w)), y = "d-excess (\u2030)")
ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv", common.legend = TRUE)





