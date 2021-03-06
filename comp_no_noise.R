# Competition in fluctuating environment, no noise

source('preamble.R')



# Logistic reponse between X and Y
resp_fun <- function(X_prop, y_start = 0, y_end = 1, shape = 100, rate = 10){
    {(y_end - y_start) / (1 + shape * exp(-rate * X_prop))} + y_start
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


# Linearly changing X variable
lin_fun <- function(time_len, start = 1, end = 0) {
    slope <- (end - start) / (time_len - 1)
    out <- slope * (1:time_len - 1) + start
    return(out)
}






sym_3sp_1phase <- function(abundances_0, X, a_sym_fun, alphas, b_Xs, K, r_is) {
    
    a_mat <- alphas  # Working version of alphas matrix
    n_spp <- nrow(alphas)
    phase_len <- length(X)
    pop_sizes <- matrix(integer((phase_len + 1) * n_spp), ncol = n_spp)
    pop_sizes[1,] <- abundances_0  # Initial abundances
    
    for (t in 1:phase_len) {
        # Symbiosis parameter: How sp 3 affects 2, and vice versa
        a_syms <- a_sym_fun(X[t])
        # Updating alphas matrix for symbiotic effects
        a_mat[3,2] <- a_syms[1]
        a_mat[2,3] <- a_syms[2]
        # All interspecific effects
        other_spp <- a_mat[cbind(c(2,3,1), 1:3)] * pop_sizes[t, c(2,3,1)] +
            a_mat[cbind(c(3,1,2), 1:3)] * pop_sizes[t, c(3,1,2)]
        # Changing rs_t based on X[t]
        rs_t <- b_Xs * X[t] + r_is
        pop_sizes_t <- log_growth(pop_sizes[t,], rs_t, K, other_spp)
        pop_sizes[(t+1),] <- pop_sizes_t
    }
    if (any(is.na(pop_sizes))) {
        warning("Population size is NA.")
    }
    return(pop_sizes[2:(phase_len+1),])
}





sym_3sp_comp_equil <- function(a_sym_fun, alphas, phase_len = 100, 
                               abundances_0 = rep(1000, 3), b_Xs = c(0.1,0.1,-0.1), 
                               K = 2000,
                               spp_names = c('competitor', 'host', 'endosymbiont'),
                               r_is = rep(0.1,3),
                               t_max = 60, X_curve = TRUE, leeway = 0) {
    
    if (X_curve){
        X <- u_curve(phase_len, phases = 1)
    } else {
        X <- lin_fun(phase_len, start = 1, end = 0)
    }
    n_spp <- length(spp_names)
    
    t0 <- Sys.time()
    
    pops_1 <- sym_3sp_1phase(abundances_0, X, a_sym_fun, alphas, b_Xs, K, r_is)
    if (any(is.na(pops_1))) {
        warning('First phase returned NAs.')
        return(pops_1)
    }
    pops_2 <- sym_3sp_1phase(as.numeric(pops_1[phase_len,]), X, a_sym_fun, 
                             alphas, b_Xs, K, r_is)
    not_equil <- any(c(abs(apply(pops_1, 2, min) - apply(pops_2, 2, min)) > leeway, 
                       abs(apply(pops_1, 2, max) - apply(pops_2, 2, max)) > leeway))
    t_diff <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
    p_i <- 2
    
    while (not_equil) {
        pops_1 <- pops_2
        pops_2 <- sym_3sp_1phase(as.numeric(pops_1[phase_len,]), X, a_sym_fun, 
                                 alphas, b_Xs, K, r_is)
        not_equil <- any(c(abs(apply(pops_1, 2, min) - apply(pops_2, 2, min)) > leeway, 
                           abs(apply(pops_1, 2, max) - apply(pops_2, 2, max)) > leeway))
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
    
    pop_sizes_mat <- cbind(abund_lo = round(apply(pops_2, 2, min, na.rm = TRUE), 2), 
                          abund_hi = round(apply(pops_2, 2, max, na.rm = TRUE), 2),
                          phases = p_i)
    pop_df <- as_data_frame(pop_sizes_mat)
    pop_df$species <- spp_names
    return(pop_df)
}




sym_3sp_comp <- function(a_sym_fun, alphas, phase_len = 100, phases = 1,
                         abundances_0 = rep(1000, 3), b_Xs = c(0.1,0.1,-0.1), K = 2000,
                         spp_names = c('competitor', 'host', 'endosymbiont'),
                         r_is = c(0.1, 0.1, 0.1),
                         X_curve = FALSE) {
    
    time_len <- as.integer(phase_len * phases)
    
    if (X_curve){
        X <- u_curve(phase_len, phases = phases)
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
                                              a_sym_fun, alphas, b_Xs, K, r_is)
    for (p in 2:phases){
        S <- (p - 1) * phase_len + 1
        E <- p * phase_len
        E_0 <- E - phase_len
        X_p <- X[S:E]
        pop_sizes[S:E,] <- sym_3sp_1phase(as.numeric(pop_sizes[E_0,]), X_p, a_sym_fun,
                                          alphas, b_Xs, K, r_is)
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





alphas <- matrix(numeric(9), nrow = 3)
alphas[2,1] <- -0.1
alphas[1,2] <- -0.2



# a_sym <- function(x, p = -0.5){a_sym_change_cont(x, p)}
# Effect of species 3 on 2, and vice versa
a_sym <- function(x, phi = 1, y_max = 0.1, y_min = -0.1){
    c(resp_fun(x, y_max, y_min, shape = 200) * phi, resp_fun(x, y_min, y_max) * phi)
}

comp_df <- sym_3sp_comp(a_sym, alphas = alphas,
                        b_Xs = c(0.1, 0.1, 0.1),
                        phases = 3, X_curve = TRUE)%>% 
    gather(species, abundance, competitor:endosymbiont)


one_comp_p <- comp_df %>%
    ggplot(aes(generation, abundance, color = species)) +
    geom_line() +
    theme_lan() +
    scale_y_continuous('Abundance', breaks = c(1000, 1500, 2000)) +
    xlab('Time') +
    scale_color_manual(values = gg_colors) +
    scale_fill_manual(values = gg_colors) +
    theme(legend.position = 'none') +
    geom_text_repel(
        data = comp_df %>% filter(generation == 120),
        aes(y = abundance, label = species),
        size = 3, segment.color = NA, fontface = 'bold',
        max.iter = 1e3L, nudge_x = c(0,0,-65),
        box.padding = unit(0.5, 'lines')) +
    geom_text(data = data_frame(generation = rep(200,2), abundance = c(1300,1200),
                                lab = c('alpha[s * ",min" ] == -0.1', 
                                        'alpha[s * ",max" ] == 0.1')),
              aes(label = lab), parse = TRUE, color = 'black')
# one_comp_p
# ggsave('fig4.pdf', one_comp_p, width = 6, height = 3, units = 'in')

# one_equil <- sym_3sp_comp_equil(a_sym, alphas, b_Xs = c(1, 1, 0.1) * 0.1,
#                                r_is = c(0.1, 0.1, 0.1))
# one_equil


equil_sims <- function(a_sym_fun, sym_maxs, sym_mins, b_Xs, b_X3s, ...) {
    par_df <- expand.grid(sym_max = sym_maxs, sym_min = sym_mins, b_X3 = b_X3s)
    par_df <- filter(par_df, sym_max >= sym_min)
    one_sim <- function(i) {
        sym_max_i <- par_df$sym_max[i]
        sym_min_i <- par_df$sym_min[i]
        b_Xs_i <- b_Xs
        b_Xs_i[3] <- par_df$b_X3[i]
        a_sym_fun_i <- function(x) { a_sym_fun(x, y_max = sym_max_i, y_min = sym_min_i) }
        out_df <- sym_3sp_comp_equil(a_sym_fun_i, b_Xs = b_Xs_i, ...)
        out_df$sym_max <- sym_max_i
        out_df$sym_min <- sym_min_i
        out_df$b_X3 <- b_Xs_i[3]
        return(out_df)
    }
    bind_rows(lapply(1:nrow(par_df), one_sim))
}


equil_df <- equil_sims(
    a_sym, 
    sym_maxs = seq(0.2, 0, -0.1), 
    sym_mins = seq(0.1, -0.1, length.out = 101), 
    b_Xs = rep(0.1,3), 
    b_X3s = c(0),#, -0.1, -0.15), 
    alphas = alphas, leeway = 0.01) %>% 
    mutate(abund_mid = mapply(function(x,y) median(c(x, y)), abund_hi, abund_lo))


plot_df <- equil_df %>% 
    mutate(b_X3 = factor(b_X3, levels = sort(unique(b_X3), decreasing = TRUE), 
                         labels = paste('beta[3] ==', 
                                        sort(unique(b_X3), decreasing = TRUE)))#,
           # sym_max = factor(sym_max, levels = sort(unique(sym_max)), 
           #               labels = paste('alpha[s * ",max"] ==', 
           #                              sort(unique(sym_max))))
           ) %>% 
    mutate_at(vars(abund_lo, abund_hi, abund_mid), funs(. / 2000))


equil_p <- plot_df %>% 
    ggplot(aes(sym_min, abund_mid, fill = species, color = species)) + 
    geom_ribbon(aes(ymin = abund_lo, ymax = abund_hi), alpha = 0.4, color = NA) +
    geom_line() +
    ylab(expression('Relative long-term abundance (' * N[infinity] / K * ')')) +
    scale_x_continuous(expression(paste('Minimum symbiotic effect (', 
                                        alpha[s * ',min'], ')')),
                       breaks = c(-0.1, 0, 0.1)) +
    facet_grid(. ~ sym_max, labeller = label_parsed, scales = 'free_x') +
    theme_lan() +
    ggtitle(expression("Maximum symbiotic effect (" * alpha[s * ",max"] * ")")) +
    scale_color_manual(values = gg_colors) +
    scale_fill_manual(values = gg_colors) +
    theme(legend.position = 'none', 
          plot.title = element_text(size = 11, face = 'plain', hjust = 0.5),
          axis.title = element_text(size = 11, face = 'plain', hjust = 0.5)) +
    geom_text_repel(
        data = plot_df %>% filter(sym_max == 0, sym_min == 0),
        aes(y = abund_hi, label = species),
        size = 3, segment.color = NA, fontface = 'bold',
        max.iter = 1e3L, nudge_x = -0.005, nudge_y = 0.01,
        box.padding = unit(0.5, 'lines'))

# equil_p
# ggsave('fig5.pdf', equil_p, width = 6, height = 4, units = 'in')

# When were median host abundances > median competitor abundances?
equil_df %>% 
    group_by(sym_max) %>% 
    summarize(big = sym_min[
        which(abund_mid[species == 'host'] >= abund_mid[species == 'competitor'])] %>% 
            tail(., 1))


# Vary b_X[3] from 0 to -b_X[2]
# Vary phi






# # Response curves for symbiotic effects (Figure 2)
# plot_df <- expand.grid(x = seq(0,1,length.out = 101), g = c(1,2), y_m = c(0, -0.1)) %>%
#     as_data_frame %>%
#     rowwise %>%
#     mutate(y = a_sym(x, p = 1, y_min = y_m, y_max = 0.1)[g]) %>%
#     ungroup %>%
#     mutate(g = factor(g, labels = c('host', 'endosymbiont')))
# resp_p <- plot_df %>%
#     ggplot(aes(x, y, color = factor(g), linetype = factor(y_m))) +
#     # geom_hline(yintercept = 0, linetype = 3, color = 'gray80') +
#     geom_line(size = 0.75) +
#     theme_lan() +
#     theme(legend.position = 'none') +
#     scale_color_manual(values = gg_colors[3:2]) +
#     scale_linetype_manual(values = c(2,1)) +
#     scale_y_continuous(expression('Symbiotic effect (' * alpha[s] * ')'),
#                        limits = c(-0.1, 0.12)) +
#     scale_x_continuous(expression('Environmental variable (' * X[t] * ')'),
#                        breaks = c(0,0.5,1)) +
#     geom_text_repel(data = plot_df %>% group_by(g) %>% filter(x == max(x), y_m == 0),
#         aes(label = g),
#         size = 4, segment.color = NA,
#         max.iter = 1e3L,
#         nudge_x = c(-0.03, -0.3), nudge_y = c(-0.05, 0.025),
#         box.padding = unit(1, 'lines'), fontface = 'bold')
# resp_p
# # ggsave('fig2.pdf', resp_p, width = 4, height = 3, units = 'in')
