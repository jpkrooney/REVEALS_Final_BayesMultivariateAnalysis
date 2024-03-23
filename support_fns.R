
#extract model fixed effect coefficients and write to file
get_model_fixedcoefs <- function(model){
    coefs <- data.frame( fixef(model) )
    coefs$outcome <- substr( rownames(coefs), 4, 7 )
    coefs$variable <- gsub('^\\_|\\_$', '', substr( rownames(coefs), 8 , 1000000L))
    coefs <- coefs %>%
        dplyr::select(outcome, variable, everything())
    coefs <- coefs %>% 
        pivot_wider(names_from = outcome,
                    values_from = c(Estimate, Q2.5, Q97.5, Est.Error))
    return(coefs)
}


get_model_fixedcoefsPCT <- function(model){
    coefs <- data.frame( fixef(model) )
    
    coefs <- coefs %>% 
        mutate(flag = substr( rownames(coefs), 0, 7 ),
               outcome = case_when(flag == "fvcmaxp" ~ "fvcpct",
                                   flag == "svcperc" ~ "svcpct",
                                   flag == "outsnip" ~ "snip",
                                   flag == "outpeak" ~ "peak",
                                   flag == "alsfrsr" ~ "alsfrsrr")) %>% 
        dplyr::select(-flag) %>% 
        data.frame()
    
    coefs$variable <- str_split_fixed(rownames(coefs), pattern = "_", n = 2)[, 2]
    
    coefs <- coefs %>%
        dplyr::select(outcome, variable, everything())
    coefs <- coefs %>% 
        pivot_wider(names_from = outcome,
                    values_from = c(Estimate, Q2.5, Q97.5, Est.Error))
    return(coefs)
}




trim_pop_effects <- function(model, wanted_pars){
    pop_effects <- data.frame( fixef(model) )
    pop_effects$par <- rownames(pop_effects)
    rownames(pop_effects) <- NULL
    pop_effects <- pop_effects %>% select(par, everything())
    pop_effects <- pop_effects %>% filter(par %in% c(wanted_pars))
    return(pop_effects)
}


