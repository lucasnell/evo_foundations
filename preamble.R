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
# `other_spp` is `alpha_12 * N_2[t]`, indicator of how sp 2 affects sp 1 and N of sp 2 at
# time t; it can be of length > 1
log_growth <- function(N_t, r, K, other_spp = 0){
    new_N <- N_t + r * N_t * {(K - N_t + other_spp) / K}
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


