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
                iter = 50000, warmup = 2500, seed = 111, na_action = "impute")

summary(net_bgms0)

summary(net_bgms0)$main

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
                

# compute energies
E <- energy_from_bgms(
  bgms_fit = net_bgms0,
  data =
    dat |> 
    select(any_of(ids)) |> 
    mutate(across(
      everything(),
      ~ case_when(
        .x == -1 ~ 0,   
        .x == 0 ~ 1,   
        .x == 1 ~ 2,    
        TRUE ~ NA_real_
      )
    )),
  reference_category = 1,
  na_as_zero = F
)


summary(E)
hist(E)

dat_bind <- cbind(dat, E)

lmfit <- lm(E ~ RWA + I(RWA^2) + 
              SDO + I(SDO^2) +
              POP + I(POP^2), dat_bind |> mutate(across(
                ids), round(1)))

summary(lmfit)

sjPlot::plot_model(lmfit, "pred", show.data = T, title = "Predicted dissonance scores (energy)")

plot_network(net_bgms, dashed = T,  groups = c("core", "peri", "peri", 
                                               "core", "core", "peri",
                                               "core", "core", "peri", 
                                               "peri", "peri", "core",
                                               "peri", "peri", "core", 
                                               "peri", "peri", "core"))


