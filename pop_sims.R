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
    # new_N <- N_t + r * N_t * (1 - N_t / K)
    new_N <- N_t * exp(r * (1 - N_t / K))
    # To avoid negative or very low abundances
    new_N <- ifelse(new_N < 0.5, 0, new_N)
    # Also rounding bc a partial individual makes no sense
    return(new_N)
}




# Sine curve with limits and set number of phases
# First function is for one x value, and is used inside the 2nd one and for plotting
u_curve_1x <- function(x, time_len, phases = 1, c_max = 1, c_min = 0) {
    phase_len <- time_len / phases
    (c_max - c_min) * 0.5 * {cos(phases * (x/phase_len) * pi) + 1} + c_min
}
u_curve <- function(time_len, phases = 1, c_max = 1, c_min = 0) {
    sapply(1:time_len, u_curve_1x, time_len, phases, c_max, c_min)
}

c_max = 1
c_min = -2
phases = 2
phase_len = 100

curve(u_curve_1x(x, time_len = phase_len*phases, phases = phases), 0, phase_len*phases, ylab='')



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
                       b1s = c(-0.1,0.1), phase_len = 100, r_sd = 0.05, K = 2000,
                       spp_names = c('independent', 'host', 'symbiont')) {
    time_len <- as.integer(phases * phase_len)
    X <- u_curve(time_len, phases = phases, c_max = 1, c_min = -1)
    n_spp <- 3
    pop_sizes <- matrix(integer(time_len*n_spp), ncol = n_spp)
    
    pop_sizes[1,] <- abundances_0  # Initial abundances
    
    for (t in 2:time_len) {
        # How sp 3 affects 2, and vice versa
        rel_t <- (t %% phase_len) / phase_len
        b2s <- b2_change_fun(rel_t, t %/% phase_len)
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
n_phases <- 12
n_sims <- 100
phase_len <- 100
# Effect of X on sp 1:2, 3
b1s <- c(0.05, -0.05)
r_sd <- 0.05
b2_change_disc <- function(rel_t, phases_done, b1_X = abs(b1s[1])) {
    i <- max(which(c(-Inf, 0.25, 0.5, 0.75) <= (rel_t + 0 * phases_done)))
    # Effect of species 3 on 2
    b3_2 <- (c(0, 0.25, 0.5, 0) * b1_X)[i]
    # Effect of species 2 on 3
    b2_3 <- (c(0.5, 0, 0, 0.25) * b1_X)[i]
    return(c(b3_2, b2_3))
}
b2_change_cont <- function(rel_t, phases_done, lambda = 0.025) {
    beta1 <- abs((2 * rel_t) - 1)
    b3_2 <- (1 - beta1) * lambda
    b2_3<- beta1 * lambda
    return(c(b3_2, b2_3))
}



# 80,000 simulated generations took 2.72 secs
sim_df <- do_sims(n_sims, sym_3sp_dr, 8, b2_change_fun = b2_change_cont, 
                  phases = n_phases, b1s = b1s, phase_len = phase_len, r_sd = r_sd) %>%
    gather(species, abundance, -generation, -sim, factor_key = TRUE)


sim_df %>% 
    ggplot(aes(generation, abundance)) +
    geom_line(aes(group = sim), alpha = 0.3, color = 'dodgerblue') +
    stat_summary(geom = 'line', fun.y = median) +
    theme_lan() +
    facet_grid( ~ species)


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
phase_len <- 100
# Effect of X on sp 1:2, 3
b1s <- c(0.5, -0.5)
r_sd <- 0.1
b2_change_cont <- function(rel_t, phases_done, lambda = max(b1s)/2) {
    beta1 <- abs((2 * rel_t) - 1)
    b3_2 <- (1 - beta1) * lambda
    b2_3<- beta1 * lambda
    return(c(b3_2, b2_3))
}

sim_df <- do_sims(n_sims, sym_3sp_dr, 7, b2_change_fun = b2_change_cont, 
                  phases = n_phases, b1s = b1s, phase_len = phase_len, 
                  r_sd = r_sd) %>%
    gather(species, abundance, -generation, -sim, factor_key = TRUE)

sim_df %>% 
    ggplot(aes(generation, abundance)) +
    geom_line(aes(group = sim), alpha = 0.3, color = 'dodgerblue') +
    stat_summary(geom = 'line', fun.y = median) +
    theme_lan() +
    facet_grid( ~ species)


sim_df %>% 
    filter(generation == (n_phases * phase_len)) %>% 
    group_by(species) %>% 
    summarize(survived = mean(abundance > 0))



# Focusing on 1 randomly chosen simulation
set.seed(7)
sim_df %>%
    filter(sim == sample(1:n_sims, 1)) %T>%
    {p_limits <<- c(c_max = max(c(min(.$abundance), diff(range(.$abundance))/5)), 
                    c_min = 0); .} %>% 
    ggplot(aes(generation, abundance, color = species)) +
    stat_function(fun = u_curve_1x,
                  args = list(phases = n_phases, time_len = n_phases * phase_len,
                              c_max = p_limits['c_max'],
                              c_min = p_limits['c_min']),
                  linetype = 2, color = 'gray80', size = 0.5, n = 1001) +
    geom_line(size = 0.75) +
    theme_lan()












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
        rel_t <- (t %% phase_len) / phase_len
        b2s <- b2_change_fun(rel_t, t %/% phase_len)
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
n_phases <- 12
n_sims <- 100
phase_len <- 100
# Effect of X on sp 1:2, 3
b1s <- c(0.05, -0.05)
r_sd <- 0.05
b2_change_cont <- function(rel_t, phases_done, lambda = 1) {
    beta1 <- abs((2 * rel_t) - 1)
    b3_2 <- (1 - beta1) * lambda
    b2_3<- beta1 * lambda
    return(c(b3_2, b2_3))
}




# 80,000 simulated generations took 2.72 secs
sim_df <- do_sims(n_sims, sym_3sp_dk, 8, b2_change_fun = b2_change_cont, phases = n_phases,
                  b1s = b1s, phase_len = phase_len, r_sd = r_sd) %>%
    gather(species, abundance, -generation, -sim, factor_key = TRUE)

sim_df %>% 
    ggplot(aes(generation, abundance)) +
    geom_line(aes(group = sim), alpha = 0.3, color = 'dodgerblue') +
    stat_summary(geom = 'line', fun.y = median) +
    theme_lan() +
    facet_grid( ~ species)


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



# Discrete-time logistic growth for time t+1
log_growth <- function(N_t, r, K){
    # new_N <- N_t + r * N_t * (1 - N_t / K)
    new_N <- N_t * exp(r * (1 - N_t / K))
    # To avoid negative abundances
    new_N <- ifelse(new_N < 0, 0, new_N)
    # Also rounding bc a partial individual makes no sense
    # return(matrix(as.integer(new_N), nrow = 1))
    return(new_N)
}

sapply(seq(1000, 2000, 1000), log_growth, N_t = 100, r = 0.1)



log_growth(100, 0.1, 000)


