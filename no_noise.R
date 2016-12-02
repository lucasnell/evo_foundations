# Simulations with no noise

source('preamble.R')



# Discrete-time logistic growth for time t+1
# `other_spp` is `alpha_12 * N_2[t]`, indicator of how sp 2 affects sp 1 and N of sp 2 at
# time t; it can be of length > 1
log_growth <- function(N_t, r, K, other_spp = 0){
    new_N <- N_t + r * N_t * {(K - N_t + other_spp) / K}
    # new_N <- N_t * exp(r * (1 - N_t / K))
    # To avoid negative and very low abundances
    new_N <- ifelse(new_N < 0.5, 0, new_N)
    return(new_N)
}

u_curve <- function(phase_len, phases = 1, c_max = 1, c_min = -1, mean_zero = FALSE) {
    time_len <- phase_len * phases
    x <- 1:(time_len)
    out <- (c_max - c_min) * 0.5 * {cos(x * (2*pi/phase_len)) + 1} + c_min
    if (mean_zero) {
        # Above function can make mean(X) *very* close to, but still != 0
        # Subtracting the mean multiple times takes care of this.
        while (mean(out) != 0) { out <- out - mean(out) }
    }
    return(out)
}

triangle_fun <- function(phase_len, phases = 1, c_max = 1, c_min = -1, mean_zero = FALSE) {
    m <- abs((c_max - c_min) / (phase_len / 2)) * c(1, -1)
    b <- c(c_min, (c_max + diff(c(c_min, c_max))))
    one_tri <- function(t){
        x <- t - 1
        i <- ifelse(x <= (phase_len/2), 1, 2)
        m[i] * x + b[i]
    }
    out <- rep(one_tri(1:phase_len), phases)
    if (mean_zero) {
        # Above function can make mean(X) *very* close to, but still != 0
        # Subtracting the mean multiple times takes care of this.
        while (mean(out) != 0) { out <- out - mean(out) }
    }
    return(out)
}



# 3 species
# 1:2 decline with X
# 3 increases with X
# 1 and 2 compete
# 2 and 3 are in mutualism with greater benefits in poor environment



sym_3sp_comp_1p <- function(abundances_0, X, b_sym_fun, alphas, phase_len, b_Xs, K) {
    
    n_spp <- 3
    pop_sizes <- matrix(integer((phase_len + 1) * n_spp), ncol = n_spp)
    pop_sizes[1,] <- abundances_0  # Initial abundances
    
    for (t in 1:phase_len) {
        # Symbiosis parameter: How sp 3 affects 2, and vice versa
        b_syms <- b_sym_fun(X[t])
        # Changing rs_t based on X[t]
        rs_t <- b_Xs[c(1,1,2)] * X[t]
        # Effects based on other species' abundances
        other_spp <- c(0, b_syms) * pop_sizes[t, c(1,3,2)] + 
            c(alphas, 0) * pop_sizes[t, c(2,1,3)]
        pop_sizes_t <- log_growth(pop_sizes[t,], rs_t, K, other_spp)
        pop_sizes[(t+1),] <- pop_sizes_t
    }
    return(pop_sizes[2:(phase_len+1),])
}





sym_3sp_comp_equil <- function(b_sym_fun, alphas, phase_len = 100, 
                               abundances_0 = rep(1000, 3), b_Xs = c(0.1,-0.1), K = 2000,
                               spp_names = c('competitor', 'host', 'symbiont'),
                               t_max = 60) {
    
    X <- u_curve(phase_len, phases = 1, mean_zero = TRUE)
    n_spp <- 3

    t0 <- Sys.time()
    
    pops_1 <- sym_3sp_comp_1p(abundances_0, X, b_sym_fun, alphas, phase_len, b_Xs, K)
    pops_2 <- sym_3sp_comp_1p(as.numeric(pops_1[phase_len,]), X, b_sym_fun, 
                              alphas, phase_len, b_Xs, K)
    not_equil <- any(c(apply(pops_1, 2, min) != apply(pops_2, 2, min), 
                       apply(pops_1, 2, max) != apply(pops_2, 2, max)))
    t_diff <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
    
    while (not_equil) {
        pops_1 <- pops_2
        pops_2 <- sym_3sp_comp_1p(as.numeric(pops_1[phase_len,]), X, b_sym_fun, 
                              alphas, phase_len, b_Xs, K)
        not_equil <- any(c(apply(pops_1, 2, min) != apply(pops_2, 2, min), 
                           apply(pops_1, 2, max) != apply(pops_2, 2, max)))
        t_diff <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
        if (t_diff >= t_max) {
            stop('Maximum time reached without equilibrium.')
        }
    }
    
    pop_sizes_df <- data_frame(min = apply(pops_2, 2, min), 
                               max = apply(pops_2, 2, max),
                               species = spp_names)
    
    return(pop_sizes_df)
}




sym_3sp_comp <- function(b_sym_fun, alphas, phase_len = 100, phases = 1,
                         abundances_0 = rep(1000, 3), b_Xs = c(0.1,-0.1), K = 2000,
                         spp_names = c('competitor', 'host', 'symbiont')) {
    
    X <- u_curve(phase_len, phases = 1, mean_zero = TRUE)
    n_spp <- length(spp_names)
    
    if (n_spp != length(abundances_0)) {
        stop(paste("The abundances_0 parameter isn't the same length as the ",
                   "spp_names parameter."))
    }
    
    time_len <- as.integer(phase_len * phases)
    
    pop_sizes <- matrix(integer(time_len*n_spp), ncol = n_spp)
    
    pop_sizes[1:phase_len,] <- sym_3sp_comp_1p(abundances_0, X, b_sym_fun, alphas, 
                                               phase_len, b_Xs, K)
    for (p in 2:phases){
        new_inds <- c((p - 1) * phase_len + 1, p * phase_len)
        old_inds <- new_inds[2] - phase_len
        pop_sizes[new_inds[1]:new_inds[2],] <- sym_3sp_comp_1p(
            as.numeric(pop_sizes[old_inds,]), X, b_sym_fun, 
            alphas, phase_len, b_Xs, K
        )
    }
    
    pop_sizes_df <- as_data_frame(pop_sizes)
    colnames(pop_sizes_df) <- spp_names
    pop_sizes_df <- pop_sizes_df %>%
        mutate(generation = 1:time_len) %>%
        select(generation, everything())
    
    return(pop_sizes_df)
}

phases = 3
phase_len = 100


rm(phases, phase_len, p, start, end, new_inds, old_inds)




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


# expand.grid(b_c = seq()) %>% as.tbl

phase_len <- 100
# Effect of X on sp 1:2, 3
b_Xs <- c(1, -1) * 0.1
b_sym_phi <- 0.15
# How symbiotic effects change with X_t
b_sym_change_cont <- function(X_t, phi = 0.15){ 
    b3_2 <- (1 - X_t)/2 * phi
    b2_3<- (X_t + 1)/2 * phi
    return(c(b3_2, b2_3))
}

# Effect of sp 2 on 1, and vice versa
alphas <- c(-0.15, -0.2)

# system.time(sym_3sp_comp_equil(b_sym_change_cont, alphas = alphas, b_Xs = b_Xs))
#  user  system elapsed 
# 0.344   0.003   0.348

# alpha_21 is competitive effect of sp 2 on sp 1
# It's assumed they all have the same abs(b_X)
equil_sims <- function(b_sym_fun, b_sym_phis, alphas_21, alphas_12, b_Xs) {
    par_df <- expand.grid(b_sym_phis, alphas_21, alphas_12, b_Xs)
    one_sim <- function(i) {
        b_sym_phi_i <- par_df$b_sym_phis[i]
        alphas_i <- c(par_df$alphas_21[i], par_df$alphas_12[i])
        b_Xs_i <- par_df$b_Xs[i] * c(1, -1)
        b_sym_fun_i <- function(..., phi = b_sym_phi_i) {
            b_sym_fun(..., phi = phi)
        }
        out_df <- sym_3sp_comp_equil(b_sym_fun_i, alphas_i, b_Xs = b_Xs_i)
        
    }
}

b_sym_phis = c(0.1, 0.2, 0.5)
alphas_21 = -0.1
alphas_12 = -0.15
b_Xs = 0.3
par_df <- expand.grid(b_sym = b_sym_phis, a_21 = alphas_21, a_12 = alphas_12, b_X = b_Xs)

out_df <- sym_3sp_comp_equil(b_sym_change_cont, alphas = c(alphas_21, alphas_12), 
                             b_Xs = c(1, -1) * b_Xs)
out_df






