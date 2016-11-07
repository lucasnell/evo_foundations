

source('preamble.R')




# 3 species
# 1:2 decline with X
# 3 increases with X
# when X reaches a threshold, r for 2 affects r for 3, and vice versa


sym_3sp_dr <- function(b2_change_fun, phases = 1, abundances_0 = rep(1000, 3), 
                       b1s = c(0.1,-0.1), phase_len = 100, r_sd = 0.05, K = 2000,
                       spp_names = c('independent', 'host', 'symbiont')) {
    time_len <- as.integer(phases * phase_len)
    X <- u_curve(time_len, phases = phases, c_max = 1, c_min = -1)
    n_spp <- 3
    pop_sizes <- matrix(integer(time_len*n_spp), ncol = n_spp)
    
    pop_sizes[1,] <- abundances_0  # Initial abundances
    
    for (t in 2:time_len) {
        # How sp 3 affects 2, and vice versa
        b2s <- b2_change_fun(X[t-1])
        # Changing rs_t based on X[t], plus adding random noise
        # Linear deterministic portion based on X[t] and other sp abundances
        lin_det <- b1s[c(1,1,2)] * X[t] + c(0, b2s) * (pop_sizes[t-1,c(1,3,2)]/ K)
        rs_t <- lin_det + rnorm(n_spp, sd = r_sd)
        pop_sizes_t <- log_growth(pop_sizes[t-1,], rs_t, K)
        if (any(is.na(pop_sizes_t))) {
            warning("Population size is NA.")
        } else if (any(pop_sizes_t > K)) {
            warning("Population exceeded K. Maybe reduce r.")
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




# =========================================================
# =========================================================
# Simulation set 1: sp 2 and sp 3 equally affected by X, but in different directions
# Cyclical X for 8 phases
# Relatively low b1
# =========================================================
# =========================================================

# ---------------
# Setting simulation parameters
# ---------------
n_phases <- 8
n_sims <- 100
phase_len <- 100
# Effect of X on sp 1:2, 3
b1s <- c(0.1, -0.1)
r_sd <- abs(b1s[1]) * 1.75
b2_change_cont <- function(X_t, phi = abs(b1s[1]) * 0.5){ 
    #rel_t, phases_done, phi = abs(b1s[1]) * 0.25) {
    # beta1 <- abs((2 * rel_t) - 1)
    b3_2 <- (1 - X_t)/2 * phi # (1 - beta1) * phi
    b2_3<- (X_t + 1)/2 * phi  # beta1 * phi
    return(c(b3_2, b2_3))
}


# 80,000 simulated generations took 2.72 secs
sim_df <- do_sims(n_sims, sym_3sp_dr, 8, b2_change_fun = b2_change_cont, 
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
    ylim(c(0, 2250)) +
    theme(strip.text = element_blank(), legend.position = 'none') +
    scale_color_manual(values = c('dodgerblue', 'red'))


sim_df %>% 
    filter(generation == (n_phases * phase_len)) %>% 
    group_by(species) %>% 
    summarize(survived = mean(abundance > 0))




# Focusing on 1 simulation where independent went extinct
sim_df %>%
    filter(sim == 3) %T>%
    {p_limits <<- c(c_max = max(c(min(.$abundance), 0 + diff(range(.$abundance))/10)), 
                    c_min = 0); .} %>% 
    ggplot(aes(generation, abundance, color = species)) +
    stat_function(fun = u_curve_1x,
                  args = list(phases = n_phases, time_len = n_phases * phase_len,
                              c_max = p_limits['c_max'],
                              c_min = p_limits['c_min']),
                  linetype = 2, color = 'gray80', size = 0.5, n = 1001) +
    geom_line(size = 0.75) +
    theme_lan()










# =========================================================
# =========================================================
# Simulation set 2: sp 2 and sp 3 equally affected by X, but in different directions
# Cyclical X for 1 phase
# =========================================================
# =========================================================

# ---------------
# Setting simulation parameters
# ---------------
n_phases <- 1
n_sims <- 100
phase_len <- 200
# Effect of X on sp 1:2, 3
b1s <- c(0.25, -0.25)
r_sd <- 0.2
b2_change_cont <- function(X_t, phi = max(b1s) * 0.5){ 
    #rel_t, phases_done, phi = abs(b1s[1]) * 0.25) {
    # beta1 <- abs((2 * rel_t) - 1)
    b3_2 <- (1 - X_t)/2 * phi # (1 - beta1) * phi
    b2_3<- (X_t + 1)/2 * phi  # beta1 * phi
    return(c(b3_2, b2_3))
}

sim_df <- do_sims(n_sims, sym_3sp_dr, 7, b2_change_fun = b2_change_cont, 
                  phases = n_phases, b1s = b1s, phase_len = phase_len, 
                  r_sd = r_sd) %>%
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
    ylim(c(0, 2250)) +
    theme(strip.text = element_blank(), legend.position = 'none') +
    scale_color_manual(values = c('dodgerblue', 'red'))


sim_df %>% filter()



sim_df %>% 
    filter(generation == (n_phases * phase_len)) %>% 
    group_by(species) %>% 
    summarize(survived = mean(abundance > 0))



# # Focusing on 1 simulation where independent went extinct
# # sim_df %>% filter(generation == (n_phases * phase_len), abundance == 0)
# 
# sim_df %>%
#     filter(sim == 26) %T>%
#     {p_limits <<- c(c_max = max(.$abundance), c_min = 0); .} %>% 
#     ggplot(aes(generation, abundance, color = species)) +
#     stat_function(fun = u_curve_1x,
#                   args = list(phases = n_phases, time_len = n_phases * phase_len,
#                               c_max = p_limits['c_max'],
#                               c_min = p_limits['c_min']),
#                   linetype = 2, color = 'gray80', size = 0.5, n = 1001) +
#     geom_line(size = 0.75) +
#     theme_lan() + 
#     theme(legend.position = c(0.8, 0.2))





# =========================================================
# =========================================================
# Simulation set 3: sp 2 and sp 3 equally affected by X, but in different directions
# Cyclical X for 8 phases
# Parasitism when environment is good for 2:3
# =========================================================
# =========================================================

# ---------------
# Setting simulation parameters
# ---------------
n_phases <- 8
n_sims <- 100
phase_len <- 100
# Effect of X on sp 1:2, 3
b1s <- c(0.1, -0.1)
r_sd <- abs(b1s[1]) * 1.75
b2_change_cont <- function(X_t, phi = 0.15){ 
    b3_2 <- (0.5 - X_t)/2 * phi
    b2_3<- (X_t + 0.5)/2 * phi
    return(c(b3_2, b2_3))
}


# 80,000 simulated generations took 2.72 secs
sim_df <- do_sims(n_sims, sym_3sp_dr, 8, b2_change_fun = b2_change_cont, 
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
    ylim(c(0, 2250)) +
    theme(strip.text = element_blank(), legend.position = 'none') +
    scale_color_manual(values = c('dodgerblue', 'red'))


sim_df %>% 
    filter(generation == (n_phases * phase_len)) %>% 
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
# =======================================================================================










# # For figure 1:
# ggplot(data = data_frame(x = seq(200), y = u_curve(200, 2, c_min = -1)), aes(x, y)) +
#     geom_line() +
#     xlab('Time') +
#     ylab('X') +
#     geom_segment(x = 0, xend = 100, y = 0, yend = 0, linetype = 2) +
#     geom_segment(data = data_frame(x = c(0, 100), xend = c(0, 100),
#                                    y = rep(-0.05,2), yend = rep(0.05,2)),
#                  aes(xend = xend, yend = yend), linetype = 1) +
#     annotate('text', x = 150, y = 1, hjust = 0.5, vjust = 1, size = 5,
#              label = 'italic(L) == 100', parse = TRUE) +
#     annotate('text', x = 50, y = 0.05, hjust = 0.5, vjust = 0, size = 5,
#              label = 'one phase') +
#     theme_lan()
