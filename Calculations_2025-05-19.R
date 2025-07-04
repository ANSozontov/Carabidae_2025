export <- FALSE
multip <- 2
# Data load ---------------------------------------------------------------
library(tidyverse)
theme_set(
    theme_bw() +
        theme(
            text = element_text(family = "serif", size = 15),
            legend.position = "bottom"
        )
)

long <- dir("data", pattern = "xlsx") %>% 
    sort(decreasing = TRUE) %>% 
    `[`(1) %>% 
    paste0("data/", .) %>% 
    readxl::read_excel(sheet = "main_data") %>% 
    filter(taxa != "no_insects") %>% 
    select(-duration, -traps) %>% 
    arrange(desc(site), taxa) %>% # -taxa
    mutate(tur = factor(tur), 
        site = fct_inorder(site),
        zone = factor(zone, 
            levels = c("fon", "bufer", "impact", "superimpact")),
        zone = fct_recode(zone, 
            "background" = "fon",
            "industrial barren" = "superimpact"),
        km = str_extract(site, "[:digit:]{1,}"),
        km = as.numeric(km), 
        year = as.factor(year), 
        .after = "site")

# 1 = turs 1 and 2 are united
wide <- long %>% 
    select(-num, -tur) %>% 
    pivot_wider(names_from = taxa, 
                values_from = abu, 
                # quantitative vartiables come from abundance (inds. per 100 trap-days)
                values_fn = sum, values_fill = 0)

nsp_100 <- wide %>% 
    select(-km) %>% 
    unite("id", year, zone, site, plot, sep = "_") %>% 
    column_to_rownames("id") %>% 
    t %>% 
    as.data.frame() %>% 
    sapply(function(x){
        x <- x[x>0]
        x <- as.numeric(as.character(x))
        if(length(x) == 0){0} else 
            if(length(x) == 1){1} else {
                iNEXT::iNEXT(x, size = 100, se = FALSE, q = 0) %>%
                    pluck("iNextEst", "size_based") %>%
                    filter(m == 100) %>%
                    pull(qD)
            }
    })

# 1 = turs 1 and 2 are united
div <- tibble(wide[,1:5], 
               abu = apply(wide[,6:ncol(wide)], 1, sum),
               nsp = apply(wide[,6:ncol(wide)], 1, function(a){length(a[a>0])}),
               nsp100 = nsp_100,
               shan= vegan::diversity(wide[,6:ncol(wide)], 
                                      MARGIN = 1, index = "shannon", base = exp(1))
) %>% 
    mutate(km2 = km^2, kmLog = log(km), .after = km) %>% 
    mutate(abuLog = log10(abu+1), .after = abu)

rm(nsp_100)
res <- list()
plots <- list()
tables <- list()

# Numbers in text ---------------------------------------------------------
tables$t0_text <- div %>% 
    group_by(site) %>% 
    summarise(
        mean_abu = mean(abu), 
        min_abu  = min(abu), 
        max_abu  = max(abu),
        .groups = "drop") %>% 
    arrange(mean_abu)
if(!export){tables$t0_text}
# Models preset ---------------------------------------------------------
models_fit <- function(x){
    x %>% 
        split(1:nrow(.)) %>%
        lapply(function(y){
            df <-  filter(div, year == y$year)
            if(str_detect(y$formula, "Segmented")){
                ###
                segmented::segmented(lm(
                    str_replace_all(y$formula, "Segmented", ""), 
                    data = df), seg.Z = ~km)
                ###
            } else {
                lm(formula = y$formula, data = df)
            }
        })
}

models_pred <- function(fits){
    lapply(fits, function(x){
        c(if("segmented" %in% class(x)){x$psi[2]}else{numeric()}, seq(1, 32, by = 0.5)) %>% 
            sort() %>% 
            unique() %>%
            tibble(
                km = .,
                km2 = .^2, 
                kmLog = log(.)
            ) %>% 
            mutate(predicted = predict(x, .))
    }
)}

manual_aic <- function (a, df_observed = NA, k = 2, simple = TRUE){
    # derived from "AIC.default"    
    lls <- logLik(a)
    # lls is logLik observer
    if(is.na(df_observed)){
        df_observed <- attr(lls, "df")
    }
    if(simple){
        -2 * as.numeric(lls) + k * df_observed
    } else {
        data.frame(
            df = df_observed, 
            AIC = -2 * as.numeric(lls) + k * df_observed)
    }
    
    # old manual AIC
    # manual.aic.old <- function(a){
    #     if(!(a$family$family %in% c("quasipoisson", "poisson"))){return(NA)}
    #     loglik <- sum(dpois(a$model[,1], fitted.values(a), log = TRUE))
    #     phi <- summary(a)$dispersion
    #     -2 * loglik + 2 * summary(a)$df[3] * phi
    # }
}

model_viz_supp <- function(df0){
    df1 <- df0 %>% 
        unite("id", formula, year, sep = " ~ ") %>% 
        pull(id) %>% 
        `names<-`(df0$pred, .) %>% 
        map(~select(.x, km, predicted)) %>% 
        map_dfr(rbind, .id = "model") %>% 
        separate(model, into = c("response", "model", "year"), sep = " ~ ")
    
    df2 <- div %>% 
        select_at(c("km", predicted = df1$response[1], "year")) %>% 
        mutate(model = NA)
    
    df1 %>% 
        ggplot(aes(x = km, y = predicted,  linetype = model)) + #
        geom_point(shape = 21, size = 3, data = df2, color = "black", alpha = 0.5) +
        geom_line() + 
        facet_wrap(~year, scales = "fixed", ncol = 1)  +
        theme(panel.grid = element_blank())
        # guides(linetype = "none")
}

# models count ------------------------------------------------------------
models <- expand_grid(
    resp = c("abu", "abuLog", "nsp", "nsp100", "shan"),
    formula = c("km", "kmSegmented", "km + km2", "kmLog"), 
    year = c(2009, 2014)) %>% 
    mutate(
        formula = paste0(resp, " ~ ", formula)) %>% 
    mutate(fits = models_fit(.), 
        pred = models_pred(fits),
        aic = map_dbl(fits, ~AIC(.x)), 
        r2 = map_dbl(fits, ~summary(.x)$adj.r.squared)
    ) %>% 
    mutate(aic2 = map_dbl(fits, ~manual_aic(.x, df_observed = 3)), 
           aic2 = case_when(str_detect(formula, "Segm") ~ aic2, TRUE ~ aic),
           .after = aic) %>% 
    arrange(resp, year, aic2) %>% 
    split(.$resp)

# Fig. S2. Models all
labs <- list(
    labs(x = NULL, y = "Log10 of individuals per 100 traps-days", 
         subtitle = "log Abundance"), 
    labs(x = NULL, y = "Number of species",
         subtitle = "Number of species"), 
    labs(x = NULL, y = "Number of species per 100 individuals",
         subtitle = "Number of species rarefied to 100 individuals"),
    labs(x = NULL, y = "Shannon index", subtitle = "Diversity")  
)
plots$s2_models.all <- map2(models[-1], labs, ~model_viz_supp(.x)+.y)

plots$s2_models.all <- gridExtra::grid.arrange(
    plots$s2_models.all[[1]],
    plots$s2_models.all[[2]], 
    plots$s2_models.all[[3]], 
    plots$s2_models.all[[4]], 
    ncol = 4
)

# Fig. 2. Models selected
df1 <- models %>% 
    map_dfr(rbind) %>% 
    unite("id", formula, year, sep = " ~ ") %>% 
    pull(id)  %>% 
    `names<-`(pull(map_dfr(models, ~select(.x, pred)), 1), .) %>% 
    map_dfr(rbind, .id = "model") %>% 
    separate(model, into = c("response", "model", "year"), sep = " ~ ") %>% 
    filter(model == 'kmSegmented', response != "abu") %>% 
    select(-km2, -kmLog, -model)

df2 <- div %>% 
    transmute(year = as.character(year), km, abuLog, nsp, nsp100, shan) %>% 
    pivot_longer(names_to = "response", values_to = "predicted", -km:-year) 
plots$f2_models.selected <- df1 %>% 
    ggplot(aes(x = km, y = predicted, fill = year, color = year, shape = year)) +
    geom_point(size = 3, data = df2, alpha = 0.5, color = "black") +
    geom_line() +
    facet_wrap(~response, scales = "free") +
    scale_shape_manual(values = c(21, 24)) + 
    scale_x_continuous(limits = c(0.1, 32)) +
    theme(panel.grid = element_blank()) +
    labs(y = NULL, x = "Distanse, km", shape = "Year: ", color = "Year: ", fill = "Year: ")
if(!export){plots$f2_models.selected}

# + Table S3. Models all
tables$ts3_models.all <- models %>% 
    map(~select(.x, formula, year, aic = aic2)) %>%
    map_dfr(~mutate_if(.x, is.numeric, ~round(.x, 2))) %>% 
    filter(str_detect(formula, "abu ~ ", negate = TRUE)) %>% 
    separate(formula, into = c("predictor", "model"), sep = " ~ ") %>% 
    pivot_wider(names_from = predictor, values_from = aic) %>% 
    arrange(year, model)
if(!export){tables$ts3_models.all }

# Table 2. Models selected comparison
tables$t2_mod.comps <- models %>% 
    `[`(-1) %>% 
    map(~.x %>% 
            filter(str_detect(formula, "Segment")) %>% 
            select(year, fits)
    ) %>% 
    map(~.x %>% 
            pluck("fits") %>% 
            `names<-`(.x$year) %>% 
            lapply(summary) %>% 
            map(~list(.x$coefficients, .x$psi)) %>% 
            map(~map(.x, ~rownames_to_column(as.data.frame(.x), "pred"))) %>% 
            map_dfr(~rbind(
                select(.x[[1]], pred, `Est.` = 2, `St.Err` = 3), 
                transmute(.x[[2]], pred, `Est.`, `St.Err`)
            ) %>% 
                filter(Est. != 0, pred %in% c("km", "psi1.km")), 
            .id = "year"
            ) %>% 
            split(.$pred) %>% 
            map(~arrange(.x, Est.))
    ) %>% 
    map_dfr(~.x %>% 
                map_dfr(rbind) %>% 
                mutate(E = paste0(round(Est., 3), "±", round(St.Err, 3)), 
                       low_CI = Est. - St.Err * multip, 
                       up_CI  = Est. + St.Err * multip, 
                       CI = paste0("(", round(low_CI, 3), "; ", round(up_CI, 3), ")"),
                       .keep = "unused") %>% 
                select(-ends_with("_CI")) %>% 
                unite("E", E, CI, sep = " || ") %>% 
                arrange(year) %>% 
                pivot_wider(names_from = year, values_from = E),
            .id = "response"
    )
if(!export){tables$t2_mod.comps}

# Rarefication ------------------------------------------------------------
library(parallel)
cl <- makeCluster(detectCores()-1)
rar <- long %>% 
    group_by(year, zone, taxa) %>% 
    summarise(abu = sum(abu), .groups = "drop") %>% 
    unite("id", year, zone, sep = "_") %>% 
    split(.$id) %>% 
    map(~sort(.x$abu[.x$abu>0])) %>% 
    parLapply(cl = cl, ., function(a){
        a |> 
            as.character() |>
            as.numeric() |>
            iNEXT::iNEXT(
                size = seq(0, 100, by = 2), 
                ### 999
                nboot = 9, 
                ### 999
                se = TRUE, 
                conf = 0.999) |>
            purrr::pluck("iNextEst", "size_based") |>
            dplyr::select(m, Method, qD, qD.LCL,   qD.UCL) 
    }) %>% 
    map_df(rbind, .id = "id") %>% 
    as_tibble() %>% 
    separate(id, into = c("year", "zone"), sep = "_") %>% 
    mutate(zone = factor(zone, ordered = TRUE, levels = levels(long$zone)))

plots$f3_rarefication <- rar %>% 
    filter(m %in% seq(0, 100, by = 2)) %>% 
    mutate(qD.LCL  = case_when(is.na(qD.LCL) ~ qD, TRUE ~ qD.LCL),
           qD.UCL  = case_when(is.na(qD.UCL) ~ qD, TRUE ~ qD.UCL)) %>% 
    ggplot(aes(x = m, y = qD, color = zone, fill = zone)) + 
    facet_wrap(~year) +
    geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha = 0.2, color = NA) +
    geom_line() + 
    scale_color_manual(values = c("darkgreen", "#ccff00", "#ff9900", "#FF3300")) +
    scale_fill_manual(values = c("darkgreen", "#ccff00", "#ff9900", "#FF3300")) +
    labs(x = "Individuals", y = "Species") + 
    theme(panel.grid = element_blank())

if(!export){plots$f3_rarefication}

# Multidimensional  --------------------------------------------------
dis <- wide %>% 
    select(-site) %>% 
    unite("ID", zone, year, km, plot, sep = "_") %>% 
    column_to_rownames("ID") %>% 
    vegan::vegdist(method = "bray", binary = FALSE)
pc <- ape::pcoa(dis)
eig <- pc$values$Eigenvalues
eig <- round(eig/sum(eig)*100, 1)
pc <- pc$vectors %>% 
    as.data.frame() %>% 
    rownames_to_column("ID") %>% 
    as_tibble() %>% 
    separate(ID, into = c("zone", "year", "site", "plot"), 
             sep = "_") %>% 
    mutate(zone = factor(zone, levels = levels(long$zone)))

L <- expand_grid(levels(long$zone), c(2009, 2014)) %>% 
    apply(1, function(a){paste0(a, collapse = "_")})

distances <- as.matrix(dis)
distances[upper.tri(distances)] <- NA
diag(distances) <- NA
distances <- distances %>%
    as.data.frame() %>%
    rownames_to_column("id1") %>%
    as_tibble() %>%
    pivot_longer(names_to = "id2", values_to = "dis", -id1) %>%
    filter(!is.na(dis)) %>%
    separate(id1, c("zone1", "year1"), sep = "_", extra = "drop") %>%
    separate(id2, c("zone2", "year2"), sep = "_", extra = "drop") %>%
    mutate(
        id1 = factor(paste0(zone1, "_", year1), levels = L, ordered = TRUE),
        id2 = factor(paste0(zone2, "_", year2), levels = L, ordered = TRUE),
        .keep = "unused") %>% 
    split(1:nrow(.)) %>% 
    lapply(function(x){
        if(x$id1 > x$id2) {
            tibble(dis = x$dis, id1 = x$id2, id2 = x$id1)
        } else {
            tibble(dis = x$dis, id1 = x$id1, id2 = x$id2)
        }
    }) %>% 
    map_dfr(rbind) %>% 
    pivot_wider(
        names_from = id2, values_from = dis, 
        values_fill = NA, values_fn = mean)

set.seed(1); res$permanova <- div %>% 
    vegan::adonis2(
    dis ~ zone * year, 
    data = ., 
    permutations = 999, 
    by = "terms")

if(!export){res$permanova} 

plots$f4_pcoa <- pc %>% 
    ggplot(
        aes(x = Axis.1, y = Axis.2, linetype = year,
            fill = zone, color = zone, shape = year)) +
    geom_point(color = "black", size = 2.5) +  
    stat_ellipse() + 
    scale_shape_manual(values = c(21, 22))+
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_fill_manual(values = c("darkgreen", "greenyellow", "orange", "red")) +
    scale_color_manual(values = c("darkgreen", "greenyellow", "orange", "red")) +
    labs(subtitle = "PCoA\nTurs are united", 
         x = paste0("Axis 1 (", eig[1], " %)"), 
         y = paste0("Axis 2 (", eig[2], " %)"), 
         # fill = "Год", color = "Год", shape = "Год"
    ) +
    theme(panel.grid = element_blank())
if(!export){
    plots$f4_pcoa
}

# Export ------------------------------------------------------------
if(export){
# Fig. S2. Models all
ggsave(
    plot = plots$s2_models.all, 
    filename = paste0("export/Fig.S2_models.all_", Sys.Date(), ".pdf"), 
    height = 6, width = 11, dpi = 600)
    
# Fig. 2. Models selected
ggsave(
    plot = plots$f2_models.selected, 
    filename = paste0("export/Fig.2_models.selected_", Sys.Date(), ".pdf"), 
    height = 9, width = 11, dpi = 600)

# + Table S3. Models all
    writexl::write_xlsx(
        tables$ts3_models.all,
        paste0("export/Tab.S2_models.selected_", Sys.Date(), ".xlsx")
    )

# Table 2. Models selected
writexl::write_xlsx(tables$t2_mod.comps,  paste0("export/Tab.2_models.comp_", Sys.Date(), ".xlsx"))
# + Fig. 3. Rarefied number of species
ggsave(
    paste0("export/Fig.3_rarefication_", Sys.Date(), ".pdf"), 
    plots$f3_rarefication, 
    width = 9, height = 5.5, dpi = 600)

# Table 3. Dominant species

# Fig. 4. PCoA ordination
ggsave(paste0("export/Fig.3_ord_", Sys.Date(), ".pdf"), 
       plot = plots$f4_pcoa,
       width = 18, height = 13, units = "cm")

# Permanova in text
res$permanova %>% 
    capture.output() %>% 
    `[`(-c(2, 3, 5, 6, 13, 14)) %>% 
    as.data.frame() %>% 
    writexl::write_xlsx(paste0("export/Tab.4a_permanova_", Sys.Date(), ".xlsx"))

# Table 4. Average distances
if(export){
    distances %>% 
        mutate_if(is.numeric, function(a){round(a, 2)}) %>% 
        writexl::write_xlsx(paste0("export/Tab.4_distances_", Sys.Date(), ".xlsx"))
}

}