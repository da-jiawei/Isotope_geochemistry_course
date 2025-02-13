library(tidyverse)
library(ggpubr)
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

## Saturated water vapor pressure - Tetens formula  ----
T_air = 20
E_s = (1.0007 + 3.46e-8) * (611.21 * exp(17.502 * T_air / (240.97 + T_air))) # unit: Pa

Ps = function(temp) {
  (1.0007 + 3.46e-8) * (611.21 * exp(17.502 * temp / (240.97 + temp)))
}

dat = data.frame(temp = seq(0, 30, 1), E_s = NA)
for (i in 1:nrow(dat)) {
  dat$E_s[i] = Ps(dat$temp[i])
}

ggplot(dat, aes(x = temp, y = E_s)) +
  geom_line() +
  theme_bw() + theme +
  labs(x = expression(paste("T (", degree, "C)")),
       y = "Vapor Pressure (Pa)")

## Craig - Gordon model ----
# sensitivity test
ctrl = function(){
  vars = list(
    n = 0.5,
    theta = 0.8, 
    RH = 0.7,
    temp = 20
  )
}
cg = function(vars){
  list2env(vars, environment())
  d18_s = 0
  d2_s = 0
  alpha18_v_l = 1 / exp(-2.0667 * 10^-3 - 0.4156/(temp+273.15) + (1.137*10^3)/(temp+273.15)^2) # Majoube (1971)
  alpha2_v_l = 1 / exp(52.612*10^-3 - 76.248/(temp+273.15) + (24.844*10^3)/(temp+273.15)^2) # Majoube (1971)
  d18_v = d18_s - 1e3*(1-alpha18_v_l)
  d2_v = d2_s - 1e3*(1-alpha2_v_l)
  Ck18 = 28.5 # unit: permil - Merlivat (1978)
  Ck2 = 25.1 # unit: permil - Merlivat (1978)
  eps18 = (1 - alpha18_v_l) * 10^3
  eps2 = (1 - alpha2_v_l) * 10^3
  d18_e = (d18_s - RH * d18_v - eps18 - (1 - RH) * theta * n * Ck18) / (1 - RH)
  d2_e = (d2_s - RH * d2_v - eps2 - (1 - RH) * theta * n * Ck2) / (1 - RH)
  results = data.frame("temp" = rep(temp),
                       "RH" = rep(RH),
                       "n" = rep(n),
                       "theta" = rep(theta),
                       "d18_s" = rep(d18_s), "d2_s" = rep(d2_s),
                       "d18_v" = rep(d18_v), "d2_v" = rep(d2_v),
                       "d18_e" = rep(d18_e), "d2_e" = rep(d2_e))
}

# temperature 
vars = ctrl()
vars$temp = seq(0, 30, 5)
sens.t = cg(vars)
p1 = ggplot(sens.t, aes(x = d18_e, y = d2_e, fill = temp)) +
  geom_abline(intercept = 10, slope = 8) +
  geom_point(aes(x = d18_s, y = d2_s), shape = 22, size = 3, fill = "white") +
  geom_point(aes(x = d18_v, y = d2_v), shape = 23, size = 3, fill = "white") +
  geom_point(shape = 21, size = 3) + 
  scale_fill_distiller(palette = "RdBu") +
  theme_bw() + theme +
  scale_y_continuous(limits = c(-300, 100)) +
  scale_x_continuous(limits = c(-50, 10)) +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"),
       fill = expression(paste("T (", degree, "C)")))
p1
# humidity
vars = ctrl()
vars$RH = seq(0, 0.9, 0.1)
sens.h = cg(vars)
p2 = ggplot(sens.h, aes(x = d18_e, y = d2_e, fill = RH*100)) +
  geom_abline(intercept = 10, slope = 8) +
  geom_point(aes(x = d18_s, y = d2_s), shape = 22, size = 3, fill = "white") +
  geom_point(aes(x = d18_v, y = d2_v), shape = 23, size = 3, fill = "white") +
  geom_point(shape = 21, size = 3) + 
  scale_fill_distiller(palette = "RdBu") +
  theme_bw() + theme +
  scale_y_continuous(limits = c(-300, 100)) +
  scale_x_continuous(limits = c(-50, 10)) +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"),
       fill = "RH (%)")
p2
# theta (rm/r) 
# 1 for a small water body whose evaporation flux does not perturb the ambient moisture significantly
vars = ctrl()
vars$theta = seq(0.5, 1, 0.1)
sens.th = cg(vars)
p3 = ggplot(sens.th, aes(x = d18_e, y = d2_e, fill = theta)) +
  geom_abline(intercept = 10, slope = 8) +
  geom_point(aes(x = d18_s, y = d2_s), shape = 22, size = 3, fill = "white") +
  geom_point(aes(x = d18_v, y = d2_v), shape = 23, size = 3, fill = "white") +
  geom_point(shape = 21, size = 3) + 
  scale_fill_distiller(palette = "RdBu") +
  theme_bw() + theme +
  scale_y_continuous(limits = c(-300, 100)) +
  scale_x_continuous(limits = c(-50, 10)) +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"),
       fill = expression(theta)) 
p3
# n
# 1 for stagnant, 0.5 for fully turbulent wind conditions
vars = ctrl()
vars$n = seq(0.5, 1, 0.1)
sens.n = cg(vars)
p4 = ggplot(sens.n, aes(x = d18_e, y = d2_e, fill = n)) +
  geom_abline(intercept = 10, slope = 8) +
  geom_point(aes(x = d18_s, y = d2_s), shape = 22, size = 3, fill = "white") +
  geom_point(aes(x = d18_v, y = d2_v), shape = 23, size = 3, fill = "white") +
  geom_point(shape = 21, size = 3) + 
  scale_fill_distiller(palette = "RdBu") +
  theme_bw() + theme +
  scale_y_continuous(limits = c(-300, 100)) +
  scale_x_continuous(limits = c(-50, 10)) +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"),
       fill = "n") 
p4
ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, align = "hv")
ggsave("figures/craig_gordon_sens.jpg", width = 7.6, height = 6)






