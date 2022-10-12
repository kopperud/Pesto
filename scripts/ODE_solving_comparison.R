# ode_testing
library(deSolve)
library(patchwork)
library(dplyr)
library(ggplot2)
library(tidyr)

source("scripts/ODE.branch.backwards.rk4.R")
source("scripts/ODE.branch.backwards.R")

E.init <- c(0.0, 0.0)
D.init <- c(1.0, 0.0)
out_rk4 <- branch.prob.backwards.ode(lambda, mu, eta, 1.0, D.init, E.init)

out_eulers <- list()
step_sizes <- c(100, 1000, 10000, 100000, 1000000)
for (i in seq_along(step_sizes)){
  out_eulers[[i]] <- branch.prob.backwards(lambda, mu, eta, 1.0, D.init, E.init, step_sizes[i]) %>%
    (function(x) c("E1" = x$E[1], "E2" = x$E[2], "D1" = x$D[1], "D2" = x$D[2], "STEP" = step_sizes[i]))
}
library(dplyr)
bind_rows(out_eulers)

euler_summary <- out_eulers %>%
  bind_rows()

p1 <- euler_summary %>%
  tidyr::gather("variable", "value", -STEP) %>%
  ggplot(aes(x = STEP, y = value)) +
  facet_wrap(~ variable, scales = "free") +
  geom_point() +
  geom_line(alpha = 0.5) +
  scale_x_log10() +
  theme_classic()

euler_diff <- euler_summary
for (var in c("E1", "E2", "D1", "D2")){
  euler_diff[[var]] <- euler_diff[[var]] - tail(out_rk4[[var]], n = 1)
}

p2 <- euler_diff %>%
  tidyr::gather("variable", "value", -STEP) %>%
  filter(STEP > 500) %>%
  ggplot(aes(x = STEP, y = value)) +
  facet_wrap(~ variable, scales = "free") +
  geom_point() +
  geom_line(alpha = 0.5) +
  scale_x_log10() +
  theme_classic() +
  geom_hline(yintercept = 0, color = "red")

p <- p1 | p2 

ggsave("figures/ode_solver_euler_versus_rungekutta4.pdf", width = 450, height = 180, units = "mm")

