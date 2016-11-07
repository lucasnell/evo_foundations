library(MASS)
library(tidyverse)
library(magrittr)
library(parallel)
library(grid)
library(car)



# Custom ggplot2 theme
theme_lan <- function(base_size = 10, base_family = 'Helvetica') {
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(
            strip.text = element_text(face = 'bold', size = 11),
            strip.background = element_rect(color = NA, fill = NA),
            panel.grid = element_blank()
        )
}


inv_logit <- function(x){
    exp(x)/(1 + exp(x))
}
logit <- function(x){
    log(x/(1 - x))
}



# Discrete-time logistic growth for time t+1
log_growth <- function(N_t, r, K){
    new_N <- N_t + r * N_t * (1 - N_t / K)
    # new_N <- N_t * exp(r * (1 - N_t / K))
    # To avoid negative or very low abundances
    new_N <- ifelse(new_N < 0.5, 0, new_N)
    return(new_N)
    # new_N <- ifelse(new_N < 0, 0, new_N)
    # return(matrix(as.integer(new_N), nrow = 1))
}




# Sine curve with limits and set number of phases
# First function is for one x value, and is used inside the 2nd one and for plotting
u_curve_1x <- function(x, time_len, phases = 1, c_max = 1, c_min = 0) {
    phase_len <- time_len / phases
    (c_max - c_min) * 0.5 * {cos(x * (2*pi/phase_len)) + 1} + c_min
}
u_curve <- function(time_len, phases = 1, c_max = 1, c_min = 0) {
    sapply(1:time_len, u_curve_1x, time_len = time_len, phases = phases, 
           c_max = c_max, c_min = c_min)
}




# Change a covariance matrix
change_cov_mat <- function(cov_matrix_0, i1, i2, new_cov){
    if (i1 == i2) { stop("Cannot change covariance with itself.") }
    new_cov_matrix <- cov_matrix_0
    new_cov_matrix[i1, i2] <- new_cov
    new_cov_matrix[i2, i1] <- new_cov
    return(new_cov_matrix)
}







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
        if (any(pop_sizes_t > K)) {
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

# This assumes that each simulation returns df of same nrows
do_sims <- function(n_sims, sim_fun, seed = NULL, ...) {
    s0 <- Sys.time()
    if (!is.null(seed)) {set.seed(seed)}
    sim_list <- lapply(1:n_sims, function(i) { sim_fun(...) })
    per_sim <- nrow(sim_list[[1]])
    sim_df <- bind_rows(sim_list)
    sim_df <- mutate(sim_df, sim = rep(1:n_sims, each = per_sim))
    s1 = Sys.time()
    cat(format(n_sims * per_sim, big.mark = ',', scientific = FALSE), 
        'simulated generations took', sprintf("%.2f", s1 - s0), 'secs')
    return(sim_df)
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






sym_3sp_dk <- function(b2_change_fun, phases = 1, abundances_0 = rep(1000, 3), 
                       b1s = c(-0.1,0.1), phase_len = 100, r_sd = 0.05, K_0 = 2000,
                       spp_names = c('independent', 'host', 'symbiont')) {
    time_len <- as.integer(phases * phase_len)
    X <- u_curve(time_len, phases = phases, c_max = 1, c_min = -1)
    n_spp <- 3
    pop_sizes <- matrix(integer(time_len*n_spp), ncol = n_spp)
    
    pop_sizes[1,] <- abundances_0  # Initial abundances
    
    for (t in 2:time_len) {
        # How sp 3 affects 2, and vice versa
        b2s <- b2_change_fun(X[t-1])
        # Changing Ks_t based on other sp abundances
        Ks_t <- K_0 + c(0, b2s) * pop_sizes[t-1,c(1,3,2)]
        # Changing rs_t based on X[t] and random noise
        rs_t <- b1s[c(1,1,2)] * X[t] + rnorm(n_spp, sd = r_sd)
        pop_sizes_t <- log_growth(pop_sizes[t-1,], rs_t, Ks_t)
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
# Simulation set 4: sp 2 and sp 3 equally affected by X, but in different directions
# Cyclical X for multiple phases
# =========================================================
# =========================================================



# ---------------
# Setting simulation parameters
# ---------------
n_phases <- 6
n_sims <- 100
phase_len <- 200
# Effect of X on sp 1:2, 3
b1s <- c(0.2, -0.2)
r_sd <- 0.2
b2_change_cont <- function(X_t, phi = max(b1s) * 1){ 
    b3_2 <- (1 - X_t)/2 * phi
    b2_3<- (X_t + 1)/2 * phi
    return(c(b3_2, b2_3))
    # return(c(0,0))
}



# 80,000 simulated generations took 2.72 secs
sim_df <- do_sims(n_sims, sym_3sp_dk, 8, b2_change_fun = b2_change_cont, phases = n_phases,
                  b1s = b1s, phase_len = phase_len, r_sd = r_sd) %>%
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
    # ylim(c(0, 2250)) +
    coord_cartesian(ylim = c(0, 4000)) +
    theme(strip.text = element_blank(), legend.position = 'none') +
    scale_color_manual(values = c('dodgerblue', 'red'))


sim_df %>% 
    filter(generation == (n_phases * phase_len)) %>% 
    group_by(species) %>% 
    summarize(survived = mean(abundance > 0))




# Focusing on 1 randomly chosen simulation
set.seed(4)
sim_df %>%
    filter(sim == sample(1:n_sims, 1)) %T>%
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
