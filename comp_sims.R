source('preamble.R')





# t_series <- function(a1, a2, len = 100) {
#     X <- numeric(len)
#     for (t in 3:len) {
#         X[t] <- a1 * X[t-1] + a2 * X[t-2] + rnorm(1)
#     }
#     return(X)
# }
# plot(1:100, t_series(0.5, 0.25), type = 'l', xlab = 'time', ylab = 'abundance')
# lines(1:100, t_series(0.5, 0.25), col = 'red')



# 3 species
# 1:2 decline with X
# 3 increases with X
# 1 and 2 compete
# 2 and 3 are in mutualism with greater benefits in poor environment


sym_3sp_comp <- function(b2_change_fun, b3s, 
                         phases = 1, abundances_0 = rep(1000, 3), 
                         b1s = c(0.1,-0.1), phase_len = 100, r_sd = 0.05, K = 2000,
                         spp_names = c('competitor', 'host', 'symbiont')) {
    time_len <- as.integer(phases * phase_len)
    X <- u_curve(time_len, phases = phases, c_max = 1, c_min = -1)
    n_spp <- 3
    pop_sizes <- matrix(integer(time_len*n_spp), ncol = n_spp)
    
    pop_sizes[1,] <- abundances_0  # Initial abundances
    
    for (t in 2:time_len) {
        # Mutualism parameter: How sp 3 affects 2, and vice versa
        b2s <- b2_change_fun(X[t-1])
        # Changing rs_t based on X[t], plus adding random noise
        # Also including effects based on other species abundances
        rs_t <- b1s[c(1,1,2)] * X[t] + c(0, b2s) * (pop_sizes[t-1,c(1,3,2)]/K) + 
            c(b3s, 0) * (pop_sizes[t-1,c(2,1,3)]/K) + rnorm(n_spp, sd = r_sd)
        pop_sizes_t <- log_growth(pop_sizes[t-1,], rs_t, K)
        if (any(is.na(pop_sizes_t) | abs(pop_sizes_t) > (10 * K))) {
            warning(paste("Population reached very high/low values.",
                    "You should play around with parameters."))
        }
        pop_sizes[t,] <- pop_sizes_t
    }
    pop_sizes_df <- as_data_frame(pop_sizes)
    colnames(pop_sizes_df) <- spp_names
    pop_sizes_df <- pop_sizes_df %>%
        mutate(generation = 1:time_len) %>%
        select(generation, everything())
    return(pop_sizes_df)
}


sym_3sp_comp_lin <- function(b2_change_fun, b3s, abundances_0 = rep(1000, 3), 
                         b1s = c(0.1,-0.1), time_len = 100, r_sd = 0.05, K = 2000,
                         spp_names = c('competitor', 'host', 'symbiont')) {
    X <- seq(1, -1, length.out = time_len)
    n_spp <- 3
    pop_sizes <- matrix(integer(time_len*n_spp), ncol = n_spp)
    
    pop_sizes[1,] <- abundances_0  # Initial abundances
    
    for (t in 2:time_len) {
        # Mutualism parameter: How sp 3 affects 2, and vice versa
        b2s <- b2_change_fun(X[t-1])
        # Changing rs_t based on X[t], plus adding random noise
        # Also including effects based on other species abundances
        rs_t <- b1s[c(1,1,2)] * X[t] + c(0, b2s) * (pop_sizes[t-1,c(1,3,2)]/K) + 
            c(b3s, 0) * (pop_sizes[t-1,c(2,1,3)]/K) + rnorm(n_spp, sd = r_sd)
        pop_sizes_t <- log_growth(pop_sizes[t-1,], rs_t, K)
        if (any(is.na(pop_sizes_t) | abs(pop_sizes_t) > (10 * K))) {
            warning(paste("Population reached very high/low values.",
                          "You should play around with parameters."))
        }
        pop_sizes[t,] <- pop_sizes_t
    }
    pop_sizes_df <- as_data_frame(pop_sizes)
    colnames(pop_sizes_df) <- spp_names
    pop_sizes_df <- pop_sizes_df %>%
        mutate(generation = 1:time_len) %>%
        select(generation, everything())
    return(pop_sizes_df)
}






# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================

# Linearly decreasing X[t] and competition

# ---------------
# Setting simulation parameters
# ---------------
n_sims <- 100
time_len <- 100
# Effect of X on sp 1:2, 3
b1s <- c(0.15, -0.15)
r_sd <- 0.05
b2_phi <- 0.25
b2_change_cont <- function(X_t, phi = b2_phi){ 
    b3_2 <- (1 - X_t)/2 * phi
    b2_3<- (X_t + 1)/2 * phi
    return(c(b3_2, b2_3))
}

# Effect of sp 2 on 1, and vice versa
b3s <- c(-0.15, -0.2)

sim_df <- do_sims(n_sims, sym_3sp_comp_lin, 2, b2_change_fun = b2_change_cont, b3s = b3s, 
                  b1s = b1s, time_len = time_len, r_sd = r_sd) %>%
    gather(species, abundance, -generation, -sim, factor_key = TRUE) %>% 
    group_by(sim, species) %>% 
    mutate(extinct = any(abundance == 0)) %>%
    ungroup



sim_df %>% 
    ggplot(aes(generation, abundance)) +
    geom_text(data = data_frame(species = unique(sim_df$species), 
                                generation = rep(0,3), 
                                abundance = rep(2100, 3)), 
              aes(label = species),
              size = 4, vjust = 0, hjust = 0, fontface = 'bold') + 
    geom_line(aes(group = sim, color = factor(extinct)), alpha = 0.3) +
    stat_summary(geom = 'line', fun.y = mean) +
    theme_lan() +
    facet_grid(species ~ .) + 
    coord_cartesian(ylim = c(0, 2250)) +
    theme(strip.text = element_blank(), legend.position = 'none') +
    scale_color_manual(values = c('dodgerblue', 'red'))




sim_df %>% 
    filter(generation == time_len) %>% 
    group_by(species) %>% 
    summarize(survived = mean(abundance > 0))







# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================



# Cyclical X[t] and competition


# ---------------
# Setting simulation parameters
# ---------------
n_phases <- 4
n_sims <- 100
phase_len <- 100
# Effect of X on sp 1:2, 3
b1s <- c(1, -1) * 0.15
r_sd <- 0.05
b2_phi <- 0.15
b2_change_cont <- function(X_t, phi = b2_phi){ 
    b3_2 <- (1 - X_t)/2 * phi
    b2_3<- (X_t + 1)/2 * phi
    return(c(b3_2, b2_3))
}

# (function(X_t, phi = 1){b3_2 <- (1.5 - X_t)/2 * phi; b2_3<- (X_t + 1)/2 * phi; return(c(b3_2, b2_3))})(1)


# Effect of sp 2 on 1, and vice versa
# b3s <- c(-1, -2.5) / 1e3
b3s <- c(-0.1, -0.2)

sim_df <- do_sims(n_sims, sym_3sp_comp, 2, b2_change_fun = b2_change_cont, b3s = b3s, 
                  phases = n_phases, b1s = b1s, phase_len = phase_len, r_sd = r_sd) %>%
    gather(species, abundance, -generation, -sim, factor_key = TRUE) %>% 
    group_by(sim, species) %>% 
    mutate(extinct = any(abundance == 0)) %>%
    ungroup



sim_df %>% 
    ggplot(aes(generation, abundance)) +
    geom_text(data = data_frame(species = unique(sim_df$species), 
                                generation = rep(0,3), 
                                abundance = rep(2100, 3)), 
              aes(label = species),
              size = 4, vjust = 0, hjust = 0, fontface = 'bold') + 
    geom_line(aes(group = sim, color = factor(extinct)), alpha = 0.3) +
    stat_summary(geom = 'line', fun.y = mean) +
    theme_lan() +
    facet_grid(species ~ .) + 
    coord_cartesian(ylim = c(0, 2250)) +
    theme(strip.text = element_blank(), legend.position = 'none') +
    scale_color_manual(values = c('dodgerblue', 'red'))






sim_df %>% 
    filter(generation == (n_phases * phase_len)) %>% 
    group_by(species) %>% 
    summarize(survived = mean(abundance > 0))




