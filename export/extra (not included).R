# PC corelation -----------------------------------------------------------
sites <- readxl::read_excel(path, "sites", skip = 2) %>% 
    mutate(
        # site = as.numeric(str_extract(site, "[:digit:]+"))
        # log_km = log10(км)
    ) %>% 
    select(-zone, -km)

pca <- sites %>% 
    mutate(id = paste0(site, "_", plot), .keep = "unused") %>% 
    column_to_rownames("id") %>% 
    # select(-mean_height:-grass_shan, -moss.n.spec625:-shrub_n.spec625) %>% 
    mutate_at(c("Cu", "Pb", "Cd"), log10) %>% 
    mutate_all(~scale(.x)[,1]) %>% 
    prcomp(., center = TRUE, scale. = TRUE)

eig <- round(pca$sdev / sum(pca$sdev) * 100, 1 )

f.loads <- pca$rotation %>% 
    as.data.frame() %>% 
    rownames_to_column("factors") %>% 
    select(1:4) %>% 
    as_tibble() %>% 
    arrange(PC1) %>% 
    mutate_if(is.numeric, ~round(.x, 3))

readxl::read_excel(path, "sites", skip = 1) %>% 
    slice(1) %>% 
    select(5:ncol(.)) %>% 
    t %>% 
    as.data.frame() %>% 
    rownames_to_column("factor_rus") %>% 
    rename(factors = V1) %>% 
    left_join(f.loads, by = "factors") %>% 
    # f.loads %>% 
    mutate(
        PC1 = formattable::color_tile("darkred","green")(PC1),
        PC2 = formattable::color_tile("darkred","green")(PC2),
        PC3 = formattable::color_tile("darkred","green")(PC3)) %>% 
    kable("html", escape = FALSE) %>%
    kable_styling(full_width = FALSE)

pca <- pca$x %>% 
    as.data.frame() %>% 
    rownames_to_column("id") %>% 
    separate(id, c("site", "plot")) %>% 
    as_tibble() %>% 
    mutate(plot = as.numeric(plot)) %>% 
    left_join(filter(div, year == "2009"), by = c("site", "plot")) %>% 
    filter(!is.na(abu)) %>% 
    transmute(site, plot, PC1, PC12 = PC1^2, abu, nsp, nsp100, shan) 

fits <- expand_grid(
    barren = c("incl", "excl"), 
    f = c("abu", "nsp", "nsp100", "shan")) %>% 
    mutate(
        f = case_when(
            barren == "excl" ~ paste(f, "~ PC1"), 
            TRUE ~ paste(f,"~ PC1 + PC12"))
    ) %>% 
    # slice(1, 5) %>%
    split(1:nrow(.)) %>% 
    lapply(function(a){
        if(a$barren == "excl"){
            d <- filter(pca, site != "K01S")
        } else {
            d <- pca
        }
        lm(a$f, data = d)
    })

tabs <- fits %>% 
    lapply(function(a){
        fits <- summary(a)
        tibble(
            barren = nrow(a$model),# a$barren, 
            index = capture.output(a$terms)[1],
            inter = round(fits$coefficients[1,1], 1),
            est = round(fits$coefficients[-1,1], 1),
            se = round(fits$coefficients[-1,2], 1),
            t.val = round(fits$coefficients[-1,3], 1),
            p.val = round(fits$coefficients[-1,4], 3),
            r.adj = round(fits$adj.r.squared, 3)
        )
    }) %>% 
    map_df(rbind) %>% 
    mutate(
        barren = case_when(barren < 29 ~ "excl", TRUE ~ "incl"),
        index = str_replace_all(index, "PC12|PC1|~| +|\\+", ""))

sq.models <- map(fits[1:4], 
                 ~tibble(PC1 = seq(-4, 4, by = 0.1), PC12 = PC1^2) %>% 
                     mutate(val = stats::predict.lm(object = .x, newdata = .))) %>% 
    `names<-`(c("abu", "nsp", "nsp100", "shan")) %>% 
    map_df(rbind, .id = "beetles")

plots$PCA_models <- pca %>% 
    filter(site != "K01S") %>% 
    pivot_longer(names_to = "beetles", values_to = "val", -1:-4) %>% 
    ggplot(aes(PC1, val)) + 
    geom_point() + 
    geom_line(data = sq.models, color = "tomato") +
    geom_abline(
        aes(slope = est, intercept = inter),
        data = mutate(filter(tabs, barren == "excl"), beetles = index),
        color = "cyan") +
    facet_wrap(~beetles, scales = "free") + 
    labs(y = NULL) +
    theme_bw()


if(export){
    ggsave("export/PCA_models.png", plots$PCA_models, 
           height = 6, width = 8, dpi = 300)
    
    
    
    writexl::write_xlsx(tabs, "export/PCA, models.xlsx")
    
} else { 
    plots$PCA_models
}

way1 <- sites %>% 
    left_join(
        transmute(pc, year, zone, site = as.numeric(site), 
                  plot = as.numeric(plot), Axis.1, Axis.2),
        ., 
        by = c("site", "plot")) %>% 
    select(-zone) %>% 
    unite("id", site, plot, sep = "_")

cor.axis.1 <- way1 %>% 
    select(-Axis.2) %>% 
    pivot_longer(names_to = "param", values_to = "val", -1:-3) %>% 
    split(.$param) %>% 
    # `[`(1:2) %>% 
    lapply(function(a){
        b = cor.test(a$Axis.1, a$val)
        tibble(
            Axis.1 = mean(range(a$Axis.1)), 
            val = mean(range(a$val)), 
            r = b$estimate, 
            p = b$p.value) %>% 
            mutate(
                r = round(r, 2), p = round(p, 3), 
                pp = case_when(p <= 0.001 ~ "***", 
                               p<= 0.01 ~ "**", 
                               p <= 0.05 ~ "*",
                               TRUE ~ ""),
                pp = paste0(r, pp))
    }) %>% 
    map_df(rbind, .id = "param")

plots$cor.axis.1 <- way1 %>% 
    select(-Axis.2) %>% 
    pivot_longer(names_to = "param", values_to = "val", -1:-3) %>% 
    ggplot(aes(val, Axis.1, color = year)) + 
    geom_point() + 
    geom_smooth(method = "lm", formula = 'y ~ x') + 
    geom_label(data = cor.axis.1, mapping = aes(
        val, Axis.1, label = pp),
        color = "black"
    ) + 
    facet_wrap(~param, scales = "free") + 
    theme_bw() + 
    theme(legend.position = "bottom") + 
    labs(x = NULL, y = NULL, subtitle = "Axis 1, 33%")

cor.axis.2 <- way1 %>% 
    select(-Axis.1) %>% 
    pivot_longer(names_to = "param", values_to = "val", -1:-3) %>% 
    split(.$param) %>% 
    # `[`(1:2) %>% 
    lapply(function(a){
        b = cor.test(a$Axis.2, a$val)
        tibble(
            Axis.1 = mean(range(a$Axis.2)), 
            val = mean(range(a$val)), 
            r = b$estimate, 
            p = b$p.value) %>% 
            mutate(
                r = round(r, 2), p = round(p, 3), 
                pp = case_when(p <= 0.001 ~ "***", 
                               p<= 0.01 ~ "**", 
                               p <= 0.05 ~ "*",
                               TRUE ~ ""),
                pp = paste0(r, pp))
    }) %>% 
    map_df(rbind, .id = "param")

plots$cor.axis.2 <- way1 %>% 
    select(-Axis.1) %>% 
    pivot_longer(names_to = "param", values_to = "val", -1:-3) %>% 
    ggplot(aes(val, Axis.2, color = year)) + 
    geom_point() + 
    geom_smooth(method = "lm", formula = 'y ~ x', 
                se = FALSE, ) + 
    geom_label(data = cor.axis.2, mapping = aes(
        val, Axis.1, label = pp),
        color = "black"
    ) + 
    facet_wrap(~param, scales = "free") + 
    theme_bw() + 
    theme(legend.position = "bottom") + 
    labs(x = NULL, y = NULL, subtitle = "Axis 2, 16%")

if(export){
    ggsave("export/corr.axis.1.png", plots$cor.axis.1, 
           height = 9, width = 11, dpi = 300)
    ggsave("export/corr.axis.2.png", plots$cor.axis.2, 
           height = 9, width = 11, dpi = 300)
    lst(axis.1 = cor.axis.1, axis.2 = cor.axis.2) %>% 
        map(~select(.x, Параметр = param, pp)) %>% 
        map_dfr(rbind, .id = "Axis") %>% 
        pivot_wider(names_from = Axis, values_from = pp) %>% 
        writexl::write_xlsx("export/PC, correlation.xlsx")
    
} else { 
    plots$cor.axis.1
    plots$cor.axis.2
}

# PCA  -----------------------------------------------------------
sites.labs <- readxl::read_excel(path, "sites", skip = 2) %>% 
    slice(1) %>% 
    select(5:23) %>% 
    t %>% 
    as.data.frame() %>% 
    select(param = 1) %>% 
    rownames_to_column("param.labs") %>% 
    filter(substr(param.labs, 1, 3) != "...") %>% 
    as_tibble()
sites <- readxl::read_excel(path, "sites", skip = 3) %>% 
    select_at(c("zone", "site", "plot", "km", sites.labs$param)) %>% 
    mutate_at((ncol(.)-2):ncol(.), log10) %>% 
    mutate(logkm = log10(km), .after = km)

result <- sites %>% 
    select(-km, -zone) %>% 
    mutate(id = paste0(site, "_", plot), .keep = "unused") %>% 
    column_to_rownames("id") %>% 
    mutate_at(c("Cu", "Pb", "Cd"), log10) %>%
    mutate_all(~scale(.x)[,1]) %>% 
    prcomp(., center = TRUE, scale. = TRUE)
eig.s <- round(result$sdev / sum(result$sdev) * 100, 1 )
f.loads <- result$rotation %>% 
    as.data.frame() %>% 
    rownames_to_column("factors") %>% 
    select(1:4) %>% 
    as_tibble() %>% 
    mutate_if(is.numeric, ~round(.x, 3))
coords <- result$x %>%
    as.data.frame() %>% 
    rownames_to_column("id") %>% 
    separate(id, c("site", "plot")) %>% 
    as_tibble() %>% 
    mutate(plot = as.numeric(plot)) %>% 
    left_join(filter(div, year == "2009"), by = c("site", "plot")) %>% 
    filter(!is.na(abu)) %>% 
    transmute(site, plot, PC1, PC12 = PC1^2, abu, nsp, nsp100, shan) 

fits <- paste(c("abu", "nsp", "nsp100", "shan"), "~ PC1") %>% 
    map(~lm(.x, data = coords)) %>% 
    `names<-`(c("abu", "nsp", "nsp100", "shan"))
linear_coef <- map_dfr(fits, coefficients, .id = "param")
linear_models <- lapply(fits, summary)
linear_models <- map_dfr(linear_models, ~.x %>% 
                             pluck(coefficients) %>% 
                             as.data.frame() %>% 
                             rownames_to_column("coefficient") %>% 
                             mutate(r2adj = .x$adj.r.squared), 
                         .id = "parameter"
) %>% 
    mutate_at(c(3:5, 7), ~round(.x, 2)) %>% 
    mutate_at(6, ~round(.x, 4))

fits2 <- paste(c("abu", "nsp", "nsp100", "shan"), "~ PC1 + PC12") %>% 
    map(~lm(.x, data = coords)) %>% 
    `names<-`(c("abu", "nsp", "nsp100", "shan"))
nonlinear_models <- lapply(fits2, summary)
nonlinear_models <- map_dfr(nonlinear_models, 
                            ~.x %>% 
                                pluck(coefficients) %>% 
                                as.data.frame() %>% 
                                rownames_to_column("coefficient") %>% 
                                mutate(r2adj = .x$adj.r.squared), 
                            .id = "parameter"
) %>% 
    mutate_at(c(3:5, 7), ~round(.x, 2)) %>% 
    mutate_at(6, ~round(.x, 4))
nonlinear <- tibble(
    PC1 = seq(min(coords$PC1), max(coords$PC1), by = 0.1), 
    PC12 = PC1^2)
nonlinear <- map_dfc(
    fits2, ~predict(.x, newdata = nonlinear)) %>% 
    cbind(nonlinear, .) %>% 
    select(-PC12) %>% 
    pivot_longer(names_to = "param", values_to = "val", -PC1)
plot <- coords %>% 
    pivot_longer(names_to = "param", values_to = "val", -1:-4) %>% 
    ggplot(aes(PC1, val)) + 
    geom_point() + 
    geom_line(data = nonlinear, color = "tomato") +
    geom_abline(
        aes(slope = PC1, intercept = `(Intercept)`),
        data = linear_coef ,
        color = "blue") +
    facet_wrap(~param, scales = "free") + 
    labs(y = NULL) +
    theme_bw()
plot

# lm_km.log <- 
sites %>% 
    mutate(
        Cu = log10(Cu), 
        Pb = log10(Pb), 
        Cd = log10(Cd), 
        # km = km
        km = log10(km)
    ) %>% 
    select(-site, -plot) %>% 
    pivot_longer(values_to = "val", names_to = "var", -c("km", "zone")) %>% 
    split(.$var) %>% 
    # `[`(1:2) %>% 
    map(~summary(lm(val ~ km, data = .x))) %>% 
    map(~as.data.frame(.x$coefficients) %>% mutate(r2adj = .x$adj.r.squared)) %>%  
    map_dfr(~filter(.x, rownames(.x) != "(Intercept)"), rbind, .id = "Variables") %>% 
    as_tibble() %>% 
    mutate_at(c(2:4, 6), ~ round(.x, 2)) %>% 
    mutate_at(5, ~ round(.x, 4))



if(export){
    ggsave("export/PCA_models4zone.png", pca4$plot,
           height = 6, width = 8, dpi = 300)
    ggsave("export/PCA_models3zone.png", pca3$plot,
           height = 6, width = 8, dpi = 300)
    
    
    writexl::write_xlsx(lst(
        data4zones = pca4$coords, 
        data3zones = pca3$coords, 
        factor_loads4zones = pca4$f.loads,
        factor_loads3zones = pca3$f.loads, 
        linear4zone = pca4$linear_models,
        linear3zone = pca3$linear_models,
        nonlinear4zone = pca4$nonlinear_models,
        nonlinear3zone = pca3$nonlinear_models
    ), "export/PCA, models.xlsx")
    
} else { 
    plots$PCA_models
}

coords <- result$x %>%
    as.data.frame() %>% 
    rownames_to_column("id") %>% 
    separate(id, c("site", "plot")) %>% 
    as_tibble() %>% 
    mutate(plot = as.numeric(plot)) %>% 
    left_join(filter(div, year == "2009"), by = c("site", "plot")) %>% 
    filter(!is.na(abu)) %>% 
    transmute(site, plot, PC2, PC22 = PC1^2, abu, nsp, nsp100, shan) 

fits <- paste(c("abu", "nsp", "nsp100", "shan"), "~ PC2") %>% 
    map(~lm(.x, data = coords)) %>% 
    `names<-`(c("abu", "nsp", "nsp100", "shan"))
linear_coef <- map_dfr(fits, coefficients, .id = "param")
linear_models <- lapply(fits, summary)
linear_models <- map_dfr(linear_models, ~.x %>% 
                             pluck(coefficients) %>% 
                             as.data.frame() %>% 
                             rownames_to_column("coefficient") %>% 
                             mutate(r2adj = .x$adj.r.squared), 
                         .id = "parameter"
) %>% 
    mutate_at(c(3:5, 7), ~round(.x, 2)) %>% 
    mutate_at(6, ~round(.x, 4))

fits2 <- paste(c("abu", "nsp", "nsp100", "shan"), "~ PC2 + PC22") %>% 
    map(~lm(.x, data = coords)) %>% 
    `names<-`(c("abu", "nsp", "nsp100", "shan"))
nonlinear_models <- lapply(fits2, summary)
nonlinear_models <- map_dfr(nonlinear_models, 
                            ~.x %>% 
                                pluck(coefficients) %>% 
                                as.data.frame() %>% 
                                rownames_to_column("coefficient") %>% 
                                mutate(r2adj = .x$adj.r.squared), 
                            .id = "parameter"
) %>% 
    mutate_at(c(3:5, 7), ~round(.x, 2)) %>% 
    mutate_at(6, ~round(.x, 4))
nonlinear <- tibble(
    PC2 = seq(min(coords$PC2), max(coords$PC2), by = 0.1), 
    PC22 = PC2^2)
nonlinear <- map_dfc(
    fits2, ~predict(.x, newdata = nonlinear)) %>% 
    cbind(nonlinear, .) %>% 
    select(-PC22) %>% 
    pivot_longer(names_to = "param", values_to = "val", -PC2)

list(
    linear_models = filter(linear_models, coefficient != "(Intercept)"),
    nonlinear_models = filter(nonlinear_models, coefficient != "(Intercept)")
) %>% writexl::write_xlsx("diversity~PCA.xls")





plot <- coords %>% 
    pivot_longer(names_to = "param", values_to = "val", -1:-4) %>% 
    ggplot(aes(PC2, val)) + 
    geom_point() + 
    geom_line(data = nonlinear, color = "tomato") +
    geom_abline(
        aes(slope = PC2, intercept = `(Intercept)`),
        data = linear_coef ,
        color = "blue") +
    facet_wrap(~param, scales = "free") + 
    labs(y = NULL) +
    theme_bw()

# p <- 
result$x %>%
    as.data.frame() %>% 
    rownames_to_column("id") %>% 
    separate(id, c("site", "plot")) %>% 
    as_tibble() %>% 
    group_by(site) %>% 
    summarise(PC1 = mean(PC1), PC2 = mean(PC2)) %>% # writexl::write_xlsx("coords.xlsx")
    # mutate(plot = as.numeric(plot)) %>% 
    # left_join(filter(div, year == "2009"), by = c("site", "plot")) %>% 
    # filter(!is.na(abu)) %>% 
    mutate(site = str_remove_all(site, "K|S|N"), 
           site = as.numeric(site)) %>% 
    ggplot(aes(PC1, PC2, label = site, color = site)) + 
    geom_text() +
    scale_color_gradient(low = "red", high = "darkgreen") + 
    # stat_ellipse() +
    # geom_path() + 
    theme_bw() + 
    labs(x = paste("PC1, ", eig.s[1], "%"), x = paste("PC2, ", eig.s[2], "%")) + 
    geom_vline(xintercept = c(-4, -2, 0, 2), color = 'black', linetype = "dashed") + 
    theme(legend.position = "none")
p
plotly::ggplotly(p)
# library(plotly)



way1 <- sites %>% 
    left_join(
        transmute(pc, year, zone, site = as.numeric(site), 
                  plot = as.numeric(plot), Axis.1, Axis.2),
        ., 
        by = c("site", "plot")) %>% 
    select(-zone) %>% 
    unite("id", site, plot, sep = "_")

cor.axis.1 <- way1 %>% 
    select(-Axis.2) %>% 
    pivot_longer(names_to = "param", values_to = "val", -1:-3) %>% 
    split(.$param) %>% 
    # `[`(1:2) %>% 
    lapply(function(a){
        b = cor.test(a$Axis.1, a$val)
        tibble(
            Axis.1 = mean(range(a$Axis.1)), 
            val = mean(range(a$val)), 
            r = b$estimate, 
            p = b$p.value) %>% 
            mutate(
                r = round(r, 2), p = round(p, 3), 
                pp = case_when(p <= 0.001 ~ "***", 
                               p<= 0.01 ~ "**", 
                               p <= 0.05 ~ "*",
                               TRUE ~ ""),
                pp = paste0(r, pp))
    }) %>% 
    map_df(rbind, .id = "param")

plots$cor.axis.1 <- way1 %>% 
    select(-Axis.2) %>% 
    pivot_longer(names_to = "param", values_to = "val", -1:-3) %>% 
    ggplot(aes(val, Axis.1, color = year)) + 
    geom_point() + 
    geom_smooth(method = "lm", formula = 'y ~ x') + 
    geom_label(data = cor.axis.1, mapping = aes(
        val, Axis.1, label = pp),
        color = "black"
    ) + 
    facet_wrap(~param, scales = "free") + 
    theme_bw() + 
    theme(legend.position = "bottom") + 
    labs(x = NULL, y = NULL, subtitle = "Axis 1, 33%")

cor.axis.2 <- way1 %>% 
    select(-Axis.1) %>% 
    pivot_longer(names_to = "param", values_to = "val", -1:-3) %>% 
    split(.$param) %>% 
    # `[`(1:2) %>% 
    lapply(function(a){
        b = cor.test(a$Axis.2, a$val)
        tibble(
            Axis.1 = mean(range(a$Axis.2)), 
            val = mean(range(a$val)), 
            r = b$estimate, 
            p = b$p.value) %>% 
            mutate(
                r = round(r, 2), p = round(p, 3), 
                pp = case_when(p <= 0.001 ~ "***", 
                               p<= 0.01 ~ "**", 
                               p <= 0.05 ~ "*",
                               TRUE ~ ""),
                pp = paste0(r, pp))
    }) %>% 
    map_df(rbind, .id = "param")

plots$cor.axis.2 <- way1 %>% 
    select(-Axis.1) %>% 
    pivot_longer(names_to = "param", values_to = "val", -1:-3) %>% 
    ggplot(aes(val, Axis.2, color = year)) + 
    geom_point() + 
    geom_smooth(method = "lm", formula = 'y ~ x', 
                se = FALSE, ) + 
    geom_label(data = cor.axis.2, mapping = aes(
        val, Axis.1, label = pp),
        color = "black"
    ) + 
    facet_wrap(~param, scales = "free") + 
    theme_bw() + 
    theme(legend.position = "bottom") + 
    labs(x = NULL, y = NULL, subtitle = "Axis 2, 16%")

if(export){
    ggsave("export/corr.axis.1.png", plots$cor.axis.1, 
           height = 9, width = 11, dpi = 300)
    ggsave("export/corr.axis.2.png", plots$cor.axis.2, 
           height = 9, width = 11, dpi = 300)
    lst(axis.1 = cor.axis.1, axis.2 = cor.axis.2) %>% 
        map(~select(.x, Параметр = param, pp)) %>% 
        map_dfr(rbind, .id = "Axis") %>% 
        pivot_wider(names_from = Axis, values_from = pp) %>% 
        writexl::write_xlsx("export/PC, correlation.xlsx")
    
} else { 
    plots$cor.axis.1
    plots$cor.axis.2
}



# corelation --------------------------------------------------------------

# factors <- cbind(km = log10(km), sites[,-1:-4])
# factors[,-1:-2] %>% 
#     cor(method = "spearman") %>% 
#     round(2)

cor.val1 <- cor(sites[,-1:-4], method = "spearman")
rownames(cor.val1) <- c(
    "Logarithmic (log10) distance to KCS", 
    paste0(sites.labs$param.labs, " (", sites.labs$param, ")"))
# cor.pval1 <- expand_grid(v1 = 4:10, v2 = 4:10) %>%
#     split(1:nrow(.)) %>%
#     lapply(function(a){
#         p = cor.test(as_vector(factors[,a$v1]), as_vector(factors[,a$v2]),
#                      method = "spearman")
#         data.frame(p = p$p.value,
#                    v1 = colnames(factors)[a$v1],
#                    v2 = colnames(factors)[a$v2])
#     }) %>%
#     map_df(tibble) %>%
#     mutate(p = p.adjust(p, method = "BY")) %>%
#     pivot_wider(names_from = v2, values_from = p) %>%
#     column_to_rownames("v1") %>%
#     as.matrix()

# colnames(cor.pval1) <- 
#     rownames(cor.pval1) <- 
#     colnames(cor.val1) <- 
#     rownames(cor.val1) <- c(
#         "light level", "wind_20cm", "wind_2m", "temperature_0cm", 
#         "temperature_2m", "humidity_relative", "humidity_absolute"
#     )

# pdf(paste0("export/multicoll_", Sys.Date(), ".pdf"), width = 4, height = 4)
corrplot::corrplot(
    corr = cor.val1, 
    is.corr = F,
    # p.mat = cor.pval1,
    type="upper", 
    order = "original",
    diag = FALSE,
    col =  corrplot::COL2('RdYlBu', 10)[10:1], 
    sig.level = 0.05)
# dev.off() 
