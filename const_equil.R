# Constant environment, no noise, looking for equilibrium values

# Simulations with no noise

source('preamble.R')


# Assuming constant environment, return equilibrium abundance (where dN/dt = 0)
equil_dens <- function(alphas, Ks) {
    k1 <- Ks[1]
    k2 <- Ks[2]
    k3 <- Ks[3]
    a21 <- alphas[2,1]
    a31 <- alphas[3,1]
    a12 <- alphas[1,2]
    a32 <- alphas[3,2]
    a13 <- alphas[1,3]
    a23 <- alphas[2,3]
    denom <- (-1 + a12 * a21 + a13 * a31 + a12 * a23 * a31 + a13 * a21 * a32 + a23 * a32)
    n1 <- - (k1 - a23 * a32 * k1 + a21 * k2 + a23 * a31 * k2 + a31 * k3 + 
                 a21 * a32 * k3)
    n2 <- - (a12 * k1 + a13 * a32 * k1 + k2 - a13 * a31 * k2 + a12 * a31 * k3 + 
                 a32 * k3)
    n3 <- - (a13 * k1 + a12 * a23 * k1 + a13 * a21 * k2 + a23 * k2 + k3 - 
                 a12 * a21 * k3)
    return(c(n1, n2, n3) / denom)
}

alpha_sim <- function(alphas_original, symbioses, competitions, Ks = rep(2000, 3),
                      spp_names = c('competitor', 'host', 'symbiont')) {
    par_df <- expand.grid(s = abs(symbioses), c = -abs(competitions))
    one_sim <- function(i, a_o, Ks) {
        s <- par_df$s[i]
        c <- par_df$c[i]
        alphas_sim <- a_o
        alphas_sim[3,2] <- s
        alphas_sim[2,3] <- s
        alphas_sim[1,2] <- c
        dens <- equil_dens(alphas_sim, Ks)
        matrix(c(dens, s, c), nrow = 1)
    }
    sim_list <- lapply(1:nrow(par_df), one_sim, Ks = Ks, a_o = alphas_original)
    sim_df <- as_data_frame(do.call(rbind, sim_list))
    names(sim_df) <- c(spp_names, 'a_sym', 'a_comp')
    return(sim_df)
}



alphas <- matrix(numeric(9), nrow = 3)
alphas[2,1] <- -0.1
a_comps <- -c(0.15, 0.5)

equil_sim_df <- alpha_sim(alphas, seq(0.0, 0.9, length.out = 1000 + 1), a_comps)



plot_df <- equil_sim_df %>% 
    gather(species, abundance, competitor:symbiont) %>% 
    mutate(a_comp = paste('alpha[12] / alpha[21] == ', a_comp / -0.1))


plot_b <- plot_df %>%
    ggplot(aes(a_sym, abundance / 2000, color = species)) +
    geom_line(size = 0.75) +
    facet_wrap( ~ a_comp, labeller = label_parsed, ncol = 1) +
    theme_lan() +
    ggtitle('B') +
    ylab(expression('Relative equilibrium abundance (' * N[e] ~ '/' ~ K * ')')) +
    xlab(expression('Symbiotic effect (' * alpha[s] * ')')) +
    scale_color_brewer(palette = 'Dark2') +
    theme(legend.position = 'none')



plot_a <- plot_df %>%
    filter(a_sym <= 0.41, a_sym >= 0.035) %>%
    ggplot(aes(a_sym, abundance / 2000, color = species)) +
    geom_vline(data = data_frame(
        y = c(0.05, 0.4), 
        a_comp = paste('alpha[12] / alpha[21] == ', unique(equil_sim_df$a_comp) / -0.1)), 
               aes(xintercept = y), linetype = 2) +
    geom_line(size = 0.75) +
    facet_wrap( ~ a_comp, labeller = label_parsed, ncol = 1) +
    theme_lan() +
    ggtitle('A') +
    ylab(expression('Relative equilibrium abundance (' * N[e] ~ '/' ~ K * ')')) +
    xlab(expression('Symbiotic effect (' * alpha[s] * ')')) +
    scale_color_brewer(palette = 'Dark2') +
    theme(legend.position = 'none') +
    geom_text_repel(
        data = plot_df %>% filter(a_sym == quantile(a_sym, 0.3), 
                                  a_comp==paste('alpha[12] / alpha[21] == ', 1.5)),
        aes(label = species),
        size = 4, segment.color = NA,
        max.iter = 1e3L, nudge_x = c(0,0,-0.1), nudge_y = c(-0.05,0,0),
        box.padding = unit(1, 'lines'), fontface = 'bold')




grid.newpage()
grid.draw(cbind(ggplotGrob(plot_a), ggplotGrob(plot_b), size = 'last'))



# Simulations for when symbiont isn't at higher population density than others
# df <- alpha_sim(alphas, seq(0.0, 0.9, length.out = 1000 + 1), a_comps,
#                           Ks = c(2000, 2000, 1000))
# df %>% filter(a_comp == -0.5)
# 
# df %>% 
#     gather(species, abundance, competitor:symbiont) %>% 
#     mutate(a_comp = paste('alpha[12] / alpha[21] == ', a_comp / -0.1)) %>%
#     # filter(a_sym <= 0.31, a_sym >= 0.035) %>%
#     ggplot(aes(a_sym, abundance / 2000, color = species, linetype = species)) +
#     geom_vline(data = data_frame(
#         y = c(0.05, 0.3), 
#         a_comp = paste('alpha[12] / alpha[21] == ', unique(equil_sim_df$a_comp) / -0.1)), 
#         aes(xintercept = y), linetype = 2) +
#     geom_line(size = 0.75) +
#     facet_wrap( ~ a_comp, labeller = label_parsed, ncol = 1) +
#     scale_linetype_manual(values = c(3, 1, 2)) +
#     theme_lan()


alphas_ <- alphas
alphas_[1,2] <- -0.5
alphas_[2,1] <- -0.1
alphas_[3,2] <- 0.5
alphas_[2,3] <- 0.5

equil_dens(alphas_, Ks = c(2000, 2000, 1000))


# For where competitor and host have equal equilibrium pop. abundances:
alpha_sim(alphas, seq(0.0, 0.5, length.out = 10000 + 1), a_comps) %>%
    group_by(a_comp) %>% 
    mutate(diff = abs(competitor - host)) %>% 
    filter(diff == min(diff))
# 0.044 and 0.288
