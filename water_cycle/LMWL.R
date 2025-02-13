library(tidyverse)
library(ggpubr)
library(lubridate)
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
data_grooming = function(dat){
  dat = dat[, c(20, 35, 39)] %>% drop_na()
  names(dat) = c("date", "param", "value")
  d18 = dat %>% filter(param == "O18")
  d18 = d18[, c(1,3)]
  names(d18) = c("date", "d18")
  d2 = dat %>% filter(param == "H2")
  d2 = d2[, c(1,3)]
  names(d2) = c("date", "d2")
  dat = left_join(d18, d2, by = "date") %>%
    mutate(month = month(date),
           d.excess = d2 - 8 * d18)
  return(dat)
}
dat = read_xlsx("data/GNIP/Waco.xlsx", sheet = 1)
dat = data_grooming(dat)
p1 = ggplot(dat, aes(x = d18, y = d2, fill = month)) +
  geom_abline(slope = 8, intercept = 10) +
  geom_point(shape = 21, size = 3) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE) +
  scale_fill_distiller(palette = "RdBu") +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"))
p2 = ggplot(dat, aes(x = d18, y = d.excess, fill = month)) +
  geom_point(shape = 21, size = 3) +
  scale_fill_distiller(palette = "RdBu") +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression("d-excess (\u2030)"))
ggarrange(p1, p2, nrow = 1, ncol = 2, common.legend = TRUE)

# isolated water body model
ctrl = function(){
  vars = list(
    d18o = -5,
    temp = 0,
    f = 0.5,
    w = 0.4,
    h = 0.7
  )
}

iwb = function(vars){
  list2env(vars, environment())
  d2o = d18o * 8 + 10
  R18o = (d18o/1000 + 1) * R18_VSMOW
  R2o = (d2o/1000 + 1) * R2_VSMOW
  alpha_eq(temp)
  # assume ambient vapor in isotope equilibrium with rainfall
  R18a = R18o / alpha18_l_v 
  R2a = R2o / alpha2_l_v
  diff18 = diffratio_18 * (1 - w) + w
  diff2 = diffratio_2 * (1 - w) + w
  u18 = 1 / (alpha18_l_v * diff18 * (1 - h)) - 1
  u2 = 1 / (alpha2_l_v * diff2 * (1 - h)) - 1
  R18ss = (h * alpha18_l_v * R18a) / (1 - alpha18_l_v * diff18 * (1 - h))
  R2ss = (h * alpha2_l_v * R2a) / (1 - alpha2_l_v * diff2 * (1 - h))
  R18f = f ^ u18 * (R18o - R18ss) + R18ss
  R2f = f ^ u2 * (R2o - R2ss) + R2ss
  R18v = (R18f - (R18a * alpha18_l_v * h)) / (alpha18_l_v * diff18 * (1 - h))
  R2v = (R2f - (R2a * alpha2_l_v * h)) / (alpha2_l_v * diff2 * (1 - h))
  d18a = (R18a / R18_VSMOW - 1) * 1000
  d18f = (R18f / R18_VSMOW - 1) * 1000
  d18v = (R18v / R18_VSMOW - 1) * 1000
  d2a = (R2a / R2_VSMOW - 1) * 1000
  d2f = (R2f / R2_VSMOW - 1) * 1000
  d2v = (R2v / R2_VSMOW - 1) * 1000
  results = data.frame("temp" = rep(temp), "f" = rep(f), "w" = rep(w), "h" = rep(h),
                       "d18a" = rep(d18a), "d2a" = rep(d2a),
                       "d18f" = rep(d18f), "d2f" = rep(d2f), "d18v" = rep(d18v), "d2v" = rep(d2v))
}

# temperature
vars = ctrl()
dat_s = iwb(vars)
for (i in 1:5) {
  vars$temp = 5 * (i-1)
  vars$f = seq(0.1, 1, 0.1)
  results = iwb(vars)
  dat_s = rbind(dat_s, results)
}
ggplot(dat_s) +
  geom_abline(slope = 8, intercept = 10) +
  geom_point(data = dat, aes(x = d18, y = d2), color = "gray") +
  geom_smooth(data = dat, aes(x = d18, y = d2), method = "lm", linetype = "dashed", se = FALSE) +
  geom_point(aes(x = d18f, y = d2f, fill = temp), shape = 21, size = 3) +
  scale_fill_distiller(palette = "RdBu") +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"),
       fill = expression(paste("T (", degree, "C)")))
# humidity
vars = ctrl()
dat_s = iwb(vars)
for (i in 1:4) {
  vars$h = 0.5 + i/10
  vars$f = seq(0.1, 1, 0.1)
  results = iwb(vars)
  dat_s = rbind(dat_s, results)
}
ggplot(dat_s) +
  geom_abline(slope = 8, intercept = 10) +
  geom_point(data = dat, aes(x = d18, y = d2), color = "gray") +
  geom_smooth(data = dat, aes(x = d18, y = d2), method = "lm", linetype = "dashed", se = FALSE) +
  geom_point(aes(x = d18f, y = d2f, fill = h*100), shape = 21, size = 3) +
  scale_fill_distiller(palette = "RdBu") +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"),
       fill = "RH (%)")
# turbulence
vars = ctrl()
dat_s = iwb(vars)
for (i in 1:5) {
  vars$w = 0.3 + i/10
  vars$f = seq(0.1, 1, 0.1)
  results = iwb(vars)
  dat_s = rbind(dat_s, results)
}
ggplot(dat_s) +
  geom_abline(slope = 8, intercept = 10) +
  geom_point(data = dat, aes(x = d18, y = d2), color = "gray") +
  geom_smooth(data = dat, aes(x = d18, y = d2), method = "lm", linetype = "dashed", se = FALSE) +
  geom_point(aes(x = d18f, y = d2f, fill = w), shape = 21, size = 3) +
  scale_fill_distiller(palette = "RdBu") +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"))
