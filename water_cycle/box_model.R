library(tidyverse)

# photosynthesis ----
delta_in = -3
temp = 22
alpha18_l_v = exp(-2.0667 * 10^-3 - 0.4156/(temp+273.15) + (1.137*10^3)/(temp+273.15)^2) # Majoube (1971)

dt = 1 # time step
time = seq(0, 1000, dt) # s

results = data.frame("time" = time)
results$delta_leaf[1] = -10
results$M_leaf[1] = 1 # mass of Oxygen in leaf - unit:gram
F_in = results$M_leaf[1] / 100  # influx - g/s
F_out = F_in

for (i in 1:(length(time)-1)) {
  results$M_leaf[i+1]  = results$M_leaf[i] + F_in - F_out
  results$delta_out[i] = (1000 + results$delta_leaf[i]) / alpha18_l_v - 1000
  Delta_leaf = ((delta_in - results$delta_leaf[i]) * F_in / results$M_leaf[i] -
    (results$delta_out[i] - results$delta_leaf[i]) * F_out / results$M_leaf[i]) * dt
  results$delta_leaf[i+1] = results$delta_leaf[i] + Delta_leaf
}
results$delta_out[i+1] = (1000 + results$delta_leaf[i+1]) / alpha18_l_v - 1000

ggplot(results) +
  geom_line(aes(x = time, y = delta_leaf), color = "dodgerblue") +
  geom_line(aes(x = time, y = delta_out), color = "firebrick1") +
  theme_bw() +
  annotate("text", x = 800, y = results$delta_out[i+1] - 2, label = expression(delta^"18"*"O"[out])) +
  annotate("text", x = 800, y = results$delta_leaf[i+1] - 2, label = expression(delta^"18"*"O"[leaf])) +
  labs(y = expression(delta^"18"*"O (\u2030)"))


# Rayleigh ----
temp = 0
alpha18_l_v = exp(-2.0667 * 10^-3 - 0.4156/(temp+273.15) + (1.137*10^3)/(temp+273.15)^2) # Majoube (1971)

dt = 1 # time step
time = seq(0, 1000, dt) # s

results = data.frame("time" = time)
results$delta_v[1] = -10
results$M_v[1] = 1 # mass of water vapor in cloud - unit:gram
results$f[1] = 1 # fraction of remaining vapor
F_out = results$M_v[1] / 1000  # influx - g/s

for (i in 1:(length(time)-1)) {
  results$M_v[i+1]  = results$M_v[i] - F_out
  results$f[i+1] = results$M_v[i+1] / results$M_v[1]
  results$delta_p[i] = (1000 + results$delta_v[i]) * alpha18_l_v - 1000
  Delta_v = - dt * (results$delta_p[i] - results$delta_v[i]) * F_out / results$M_v[i]
  results$delta_v[i+1] = results$delta_v[i] + Delta_v
}
results$delta_p[i+1] = (1000 + results$delta_v[i+1]) * alpha18_l_v - 1000

ggplot(results) +
  geom_line(aes(x = time, y = delta_v), color = "dodgerblue") +
  geom_line(aes(x = time, y = delta_p), color = "firebrick1") +
  theme_bw() +
  annotate("text", x = 50, y = results$delta_p[1] + 4, label = expression(delta^"18"*"O"[p])) +
  annotate("text", x = 50, y = results$delta_v[1] + 4, label = expression(delta^"18"*"O"[v])) +
  labs(y = expression(delta^"18"*"O (\u2030)"))

ggplot(results) +
  geom_line(aes(x = f, y = delta_v), color = "dodgerblue") +
  geom_line(aes(x = f, y = delta_p), color = "firebrick1") +
  theme_bw() +
  annotate("text", x = 0.9, y = results$delta_p[1] + 4, label = expression(delta^"18"*"O"[p])) +
  annotate("text", x = 0.9, y = results$delta_v[1] + 4, label = expression(delta^"18"*"O"[v])) +
  labs(y = expression(delta^"18"*"O (\u2030)"))

