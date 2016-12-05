# Simulations with no noise
# DEPRECATED

source('preamble.R')



# Logistic reponse between X and Y
resp_fun <- function(X_prop, y_min = 0, y_max = 1, shape = 100, rate = 10){
    {(y_max - y_min) / (1 + shape * exp(-rate * X_prop))} + y_min
}


u_curve <- function(phase_len, phases = 1, y_max = 1, y_min = 0, exact_mean = NULL) {
    time_len <- phase_len * phases
    x <- 1:(time_len)
    out <- (y_max - y_min) * 0.5 * {cos(x * (2*pi/phase_len)) + 1} + y_min
    if (!is.null(exact_mean)) {
        # Above function can make mean(X) *very* close to, but still != 0
        # Subtracting the mean multiple times takes care of this.
        while (mean(out) != 0) { out <- out - mean(out) }
        out <- out + exact_mean
    }
    return(out)
}

triangle_fun <- function(phase_len, phases = 1, y_max = 1, y_min = 0, exact_mean = NULL) {
    m <- abs((y_max - y_min) / (phase_len / 2)) * c(1, -1)
    b <- c(y_min, (y_max + diff(c(y_min, y_max))))
    one_tri <- function(t){
        x <- t - 1
        i <- ifelse(x <= (phase_len/2), 1, 2)
        m[i] * x + b[i]
    }
    out <- rep(one_tri(1:phase_len), phases)
    if (!is.null(exact_mean)) {
        # Above function can make mean(X) *very* close to, but still != 0
        # Subtracting the mean multiple times takes care of this.
        while (mean(out) != 0) { out <- out - mean(out) }
        out <- out + exact_mean
    }
    return(out)
}

# Linearly changing X variable
lin_fun <- function(time_len, start = 1, end = 0) {
    slope <- (end - start) / (time_len - 1)
    out <- slope * (1:time_len - 1) + start
    return(out)
}






# 3 species
# 1:2 decline with X
# 3 increases with X
# 1 and 2 compete
# 2 and 3 are in mutualism with greater benefits in poor environment



sym_3sp_1phase <- function(abundances_0, X, a_sym_fun, alphas, phase_len, b_Xs, K) {
    
    a_mat <- alphas  # Working version of alphas matrix
    n_spp <- 3
    pop_sizes <- matrix(integer((phase_len + 1) * n_spp), ncol = n_spp)
    pop_sizes[1,] <- abundances_0  # Initial abundances
    
    for (t in 1:phase_len) {
        # Symbiosis parameter: How sp 3 affects 2, and vice versa
        a_syms <- a_sym_fun(X[t])
        # Updating alphas matrix for symbiotic effects
        a_mat[3,2] <- a_syms[1]
        a_mat[2,3] <- a_syms[2]
        # All interspecific effects
        other_spp <- a_mat[cbind(1:3, c(2,3,1))] * pop_sizes[t, c(2,3,1)] +
            a_mat[cbind(1:3, c(3,1,2))] * pop_sizes[t, c(3,2,1)]
        # Changing rs_t based on X[t]
        rs_t <- b_Xs * X[t]
        pop_sizes_t <- log_growth(pop_sizes[t,], rs_t, K, other_spp)
        pop_sizes[(t+1),] <- pop_sizes_t
    }
    if (any(is.na(pop_sizes))) {
        warning("Population size is NA.")
    }
    return(pop_sizes[2:(phase_len+1),])
}





sym_3sp_comp_equil <- function(a_sym_fun, alphas, phase_len = 100, 
                               abundances_0 = rep(1000, 3), b_Xs = c(0.1,-0.1), K = 2000,
                               spp_names = c('competitor', 'host', 'symbiont'),
                               t_max = 60) {
    
    # X <- u_curve(phase_len, phases = 1, exact_mean = 0.5)
    X <- lin_fun(phase_len, start = 1, end = 0)
    n_spp <- 3

    t0 <- Sys.time()
    
    pops_1 <- sym_3sp_1phase(abundances_0, X, a_sym_fun, alphas, phase_len, b_Xs, K)
    if (any(is.na(pops_1))) {
        warning('First phase returned NAs.')
        return(pops_1)
    }
    pops_2 <- sym_3sp_1phase(as.numeric(pops_1[phase_len,]), X, a_sym_fun, 
                              alphas, phase_len, b_Xs, K)
    not_equil <- any(c(apply(pops_1, 2, min) != apply(pops_2, 2, min), 
                       apply(pops_1, 2, max) != apply(pops_2, 2, max)))
    t_diff <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
    p_i <- 2
    
    while (not_equil) {
        pops_1 <- pops_2
        pops_2 <- sym_3sp_1phase(as.numeric(pops_1[phase_len,]), X, a_sym_fun, 
                              alphas, phase_len, b_Xs, K)
        not_equil <- any(c(apply(pops_1, 2, min) != apply(pops_2, 2, min), 
                           apply(pops_1, 2, max) != apply(pops_2, 2, max)))
        p_i <- p_i + 1
        if (any(is.na(pops_2))) {
            warning('Phase ', p_i, ' returned NAs.')
            break
        }
        t_diff <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
        if (t_diff >= t_max) {
            stop('Maximum time reached without equilibrium at phase ', p_i)
        }
    }
    
    pop_sizes_df <- data_frame(min = apply(pops_2, 2, min, na.rm = TRUE), 
                               max = apply(pops_2, 2, max, na.rm = TRUE),
                               species = spp_names,
                               phases = p_i)
    
    return(pop_sizes_df)
}





sym_3sp_comp <- function(a_sym_fun, alphas, phase_len = 100, phases = 1,
                         abundances_0 = rep(1000, 3), b_Xs = c(0.1,-0.1), K = 2000,
                         spp_names = c('competitor', 'host', 'symbiont'),
                         X_curve = FALSE) {
    
    time_len <- as.integer(phase_len * phases)
    
    if (X_curve){
        X <- u_curve(phase_len, phases = phases, exact_mean = 0.5)
    } else {
        X <- lin_fun(time_len, start = 1, end = 0)
    }
    
    n_spp <- length(spp_names)
    
    if (n_spp != length(abundances_0)) {
        stop(paste("The abundances_0 parameter isn't the same length as the ",
                   "spp_names parameter."))
    }
    
    pop_sizes <- matrix(integer(time_len*n_spp), ncol = n_spp)
    
    pop_sizes[1:phase_len,] <- sym_3sp_1phase(abundances_0, X[1:phase_len], 
                                               a_sym_fun, alphas, 
                                               phase_len, b_Xs, K)
    for (p in 2:phases){
        S <- (p - 1) * phase_len + 1
        E <- p * phase_len
        E_0 <- E - phase_len
        X_p <- X[S:E]
        pop_sizes[S:E,] <- sym_3sp_1phase(as.numeric(pop_sizes[E_0,]), X_p, a_sym_fun,
                                           alphas, phase_len, b_Xs, K)
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





# # Cyclical X[t] and competition
# 
# 
# # expand.grid(b_c = seq()) %>% as.tbl
# 
# phase_len <- 100
# # Effect of X on sp 1:2, 3
# b_Xs <- c(1, -1) * 0.1
# a_sym_phi <- 0.15
# # How symbiotic effects change with X_t
# a_sym_change_cont <- function(X_t, phi){ 
#     b3_2 <- (1 - X_t)/2 * phi
#     b2_3<- (X_t + 1)/2 * phi
#     return(c(b3_2, b2_3))
# }
# 
# 
# 
# 
# # Effect of sp 2 on 1, and vice versa
# alphas <- c(-0.15, -0.2)
# 
# # system.time(sym_3sp_comp_equil(a_sym_change_cont, alphas = alphas, b_Xs = b_Xs))
# #  user  system elapsed 
# # 0.344   0.003   0.348
# 
# # alpha_21 is competitive effect of sp 2 on sp 1
# # It's assumed they all have the same abs(b_X)
# equil_sims <- function(a_sym_fun, a_sym_phis, alphas_21, alphas_12, b_Xs) {
#     par_df <- expand.grid(a_sym = a_sym_phis, a_21 = alphas_21, a_12 = alphas_12,
#                           b_X = b_Xs)
#     one_sim <- function(i) {
#         a_sym_i <- par_df$a_sym[i]
#         alphas_i <- c(par_df$a_21[i], par_df$a_12[i])
#         b_Xs_i <- par_df$b_X[i] * c(1, -1)
#         a_sym_fun_i <- function(x, bs = a_sym_i) { a_sym_fun(x, bs) }
#         out_df <- sym_3sp_comp_equil(a_sym_fun_i, alphas_i, b_Xs = b_Xs_i)
#         out_df$a_sym <- a_sym_i
#         out_df$a_21 <- alphas_i[1]
#         out_df$a_12 <- alphas_i[2]
#         out_df$b_X <- abs(b_Xs_i[1])
#         return(out_df)
#     }
#     bind_rows(lapply(1:nrow(par_df), one_sim))
# }
# 
# a_sym_phis = seq(2, 6, 2)
alphas <- matrix(numeric(9), nrow = 3)
alphas[2,1] = -0.1
alphas[1,2] = -0.15
b_Xs = 0.1
# # par_df <- expand.grid(a_sym = a_sym_phis, a_21 = alphas_21, a_12 = alphas_12, 
# #                       b_X = b_Xs)
# 
# # out_df <- sym_3sp_comp_equil(a_sym, alphas = c(alphas_21, alphas_12),
# #                              b_Xs = c(1, -1) * b_Xs)
# # out_df
# 
# eq_df <- equil_sims(a_sym_change_cont, a_sym_phis, alphas_21, alphas_12, b_Xs)
# eq_df


# a_sym <- function(x, p = -0.5){a_sym_change_cont(x, p)}
a_sym <- function(x, p = 0.5){
    # c(resp_fun(x) * p, resp_fun(x, 0, 1, shape = 200) * p)
    c(p, p)
}

df <- sym_3sp_comp(a_sym, 
                   # alphas = matrix(numeric(9), nrow = 3),
                   alphas = alphas,
                   b_Xs = c(1, 1, 0.1) * 0.005, 
                   phases = 2, phase_len = 100)
df %>% gather(species, abundance, competitor:symbiont) %>%
    ggplot(aes(generation, abundance, color = species)) +
    geom_line() +
    theme_lan()# +
    # ylim(c(0, 2000))



# warnings()


system.time({x <- sym_3sp_comp(a_sym, 
             # alphas = matrix(numeric(9), nrow = 3),
             alphas = alphas,
             b_Xs = c(1, 1, 0.1) * 0.005, 
             phases = 3, phase_len = 100)})









# abundances_0 = rep(1000, 3)
# X = u_curve(phase_len, phases = 1, mean_zero = TRUE)
# a_sym_fun = function(x, p = -1){a_sym_change_cont(x, p)}
# alphas = c(alphas_21, alphas_12)
# b_Xs = c(1, -1) * b_Xs
# phase_len = 100
# K = 2000
# 
# n_spp <- 3
# pop_sizes <- matrix(integer((phase_len + 1) * n_spp), ncol = n_spp)
# pop_sizes[1,] <- abundances_0  # Initial abundances
# 
# for (t in 1:phase_len) {
#     # Symbiosis parameter: How sp 3 affects 2, and vice versa
#     a_syms <- a_sym_fun(X[t])
#     # Changing rs_t based on X[t]
#     rs_t <- b_Xs[c(1,1,2)] * X[t]
#     # Effects based on other species' abundances
#     other_spp <- c(0, a_syms) * pop_sizes[t, c(1,3,2)] + 
#         c(alphas, 0) * pop_sizes[t, c(2,1,3)]
#     pop_sizes_t <- log_growth(pop_sizes[t,], rs_t, K, other_spp)
#     pop_sizes[(t+1),] <- pop_sizes_t
# }
# 
# as_data_frame(pop_sizes) %>% mutate(generation = 0:100) %>% 
#     gather(species, abundance, V1:V3) %>% 
#     ggplot(aes(generation, abundance, color = species)) +
#     geom_line(aes(linetype = species)) +
#     theme_lan() +
#     ylim(c(0, 2000))


