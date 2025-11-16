library(bgms)
library(easybgm)

dat <- df_clean %>%
  mutate(
    across(
      all_of(ids),
      ~ case_when(
        .x %in% c(1, 2) ~ -1,
        .x == 3          ~  0,
        .x %in% c(4, 5) ~  1,
        TRUE            ~ NA_real_   
      )
    )
  )


net_bgms0 <- bgm(dat[, c(ids)] |> round(1),
                variable_type = c(
                rep("blume-capel", length(ids))),
                baseline_category = c(rep(0,18)), 
                iter = 50000, warmup = 2500, seed = 111)

summary(net_bgms0)

plot_network(net_bgms0, dashed = T,  groups = c("core", "peri", "peri", 
                                               "core", "core", "peri",
                                               "core", "core", "peri", 
                                               "peri", "peri", "core",
                                               "peri", "peri", "core", 
                                               "peri", "peri", "core"))

alpha_beta0 <- bc_alpha_beta_table(net_bgms0,
                                  var_names = ids)

bc_alpha_beta_plot(alpha_beta0,
                   core_alpha = 2,
                   core_beta  = 1.0)



net_bgms <- bgm(dat[, c(ids, pol_dims)] |> round(1),
                variable_type = c(
                  rep("blume-capel", length(ids)),
                  rep("ordinal", 3)), baseline_category = c(rep(0,18), rep(0,3)),
                iter = 5000, warmup = 2500, seed = 111)
                

plot_network(net_bgms, dashed = T,  groups = c("core", "peri", "peri", 
                                               "core", "core", "peri",
                                               "core", "core", "peri", 
                                               "peri", "peri", "core",
                                               "peri", "peri", "core", 
                                               "peri", "peri", "core"))

# compute energies
energies <- energy_from_bgms(
  bgms_fit          = net_bgms,
  data              = dat[, c(ids,pol_dims)],
  reference_category = 0,
  na_as_zero        = TRUE
)

summary(energies)
hist(energies)

dat_bind <- cbind(dat, energies)

lmfit <- lm(energies ~ RWA + I(RWA^2) + 
              SDO + I(SDO^2) +
              POP + I(POP^2), dat_bind)

summary(lmfit)

sjPlot::plot_model(lmfit, "pred", show.data = T, title = "Predicted dissonance scores (energy)")
