Fig3\_MixingModel
================
Tyler Gagne
6/8/2017

**Big MC/MC loop to run and generate plots for each species**

``` r
sppx <- c("SOTE", "WTSH", "BRBO", "BRNO", "LAAL", "BUPE", "WTTR", "WHTE")
# sppx<- 'SOTE'
data_box = NULL

for (x in 1:length(sppx)) {
    mix <- read.csv("mixture_data_may22.csv")
    str(mix)
    mix <- subset(mix, spp == sppx[x])
    write.csv(mix, file = "mixture_data_spp.csv")
    
    # mix data = i.e. consumer
    mix <- load_mix_data(filename = "mixture_data_spp.csv", iso_names = "TL", 
        factors = "spp", fac_random = TRUE, fac_nested = NULL, cont_effects = "year")
    # source data #new omma and caran TL
    source <- load_source_data(filename = "source_data_4grp_may22.csv", source_factors = NULL, 
        conc_dep = FALSE, data_type = "means", mix)
    # discrimination data
    discr <- load_discr_data(filename = "discrimination_data_4grp.csv", mix)
    # one dimensional isospace plot plot_data(filename='isospace_plot',
    # plot_save_pdf=FALSE, plot_save_png=FALSE, mix,source,discr) uninformative
    # prior = alpha.prior = 1 construct informative prior from
    # stomach_proportional_plot_apr13.xlsx
    TP_prior <- c(0.508, 1.212, 1.175, 1.105)
    # plot_prior(alpha.prior = TP_prior, source, plot_save_pdf=FALSE,
    # plot_save_png=FALSE) write jags model
    model_filename <- "MixSIAR_model.txt"
    resid_err <- FALSE
    process_err <- TRUE
    write_JAGS_model(model_filename, resid_err, process_err, mix, source)
    # run model
    jags.1 <- run_model(run = "test", mix, source, discr, model_filename, alpha.prior = TP_prior, 
        resid_err, process_err)
    ########################## 
    R2jags::attach.jags(jags.1)
    n.sources <- source$n.sources
    source_names <- source$source_names
    
    fac.lab <- mix$FAC[[1]]$labels
    label <- mix$cont_effects
    cont <- mix$CE[[1]]  #either CE or CE_orig explore implications?
    ilr.cont <- get(paste("ilr.cont", 1, sep = ""))
    
    get_high <- function(x) {
        return(quantile(x, 0.95))
    }
    get_low <- function(x) {
        return(quantile(x, 0.05))
    }
    
    n.plot = 124  #200 was original consider changing to number of years? which I did
    chain.len = dim(p.global)[1]
    Cont1.plot <- seq(from = round(min(cont), 1), to = round(max(cont), 1), 
        length.out = n.plot)
    ilr.plot <- array(NA, dim = c(n.plot, n.sources - 1, chain.len))
    ilr.median <- array(NA, dim = c(n.plot, n.sources - 1))
    ilr.low <- array(NA, dim = c(n.plot, n.sources - 1))
    ilr.high <- array(NA, dim = c(n.plot, n.sources - 1))
    for (src in 1:n.sources - 1) {
        for (i in 1:n.plot) {
            ilr.plot[i, src, ] <- ilr.global[, src] + ilr.cont[, src] * Cont1.plot[i]  #+ ilr.fac1[,f1,src]
            ilr.low[i, src] <- get_low(ilr.plot[i, src, ])
            ilr.median[i, src] <- median(ilr.plot[i, src, ])  #changed from median to mean 
            ilr.high[i, src] <- get_high(ilr.plot[i, src, ])
        }
    }
    
    # Transform regression lines from ILR-space to p-space
    e <- matrix(rep(0, n.sources * (n.sources - 1)), nrow = n.sources, ncol = (n.sources - 
        1))
    for (i in 1:(n.sources - 1)) {
        e[, i] <- exp(c(rep(sqrt(1/(i * (i + 1))), i), -sqrt(i/(i + 1)), rep(0, 
            n.sources - i - 1)))
        e[, i] <- e[, i]/sum(e[, i])
    }
    cross.med <- array(data = NA, dim = c(n.plot, n.sources, n.sources - 1))  # dummy variable for inverse ILR calculation
    tmp.p.med <- array(data = NA, dim = c(n.plot, n.sources))  # dummy variable for inverse ILR calculation
    p.median <- array(data = NA, dim = c(n.plot, n.sources))
    cross.low <- array(data = NA, dim = c(n.plot, n.sources, n.sources - 1))  # dummy variable for inverse ILR calculation
    tmp.p.low <- array(data = NA, dim = c(n.plot, n.sources))  # dummy variable for inverse ILR calculation
    p.low <- array(data = NA, dim = c(n.plot, n.sources))
    cross.high <- array(data = NA, dim = c(n.plot, n.sources, n.sources - 1))  # dummy variable for inverse ILR calculation
    tmp.p.high <- array(data = NA, dim = c(n.plot, n.sources))  # dummy variable for inverse ILR calculation
    p.high <- array(data = NA, dim = c(n.plot, n.sources))
    for (i in 1:n.plot) {
        for (j in 1:(n.sources - 1)) {
            cross.med[i, , j] <- (e[, j]^ilr.median[i, j])/sum(e[, j]^ilr.median[i, 
                j])
            cross.low[i, , j] <- (e[, j]^ilr.low[i, j])/sum(e[, j]^ilr.low[i, 
                j])
            cross.high[i, , j] <- (e[, j]^ilr.high[i, j])/sum(e[, j]^ilr.high[i, 
                j])
        }
        for (src in 1:n.sources) {
            tmp.p.med[i, src] <- prod(cross.med[i, src, ])
            tmp.p.low[i, src] <- prod(cross.low[i, src, ])
            tmp.p.high[i, src] <- prod(cross.high[i, src, ])
        }
        for (src in 1:n.sources) {
            p.median[i, src] <- tmp.p.med[i, src]/sum(tmp.p.med[i, ])
            p.low[i, src] <- tmp.p.low[i, src]/sum(tmp.p.low[i, ])
            p.high[i, src] <- tmp.p.high[i, src]/sum(tmp.p.high[i, ])
        }
    }
    colnames(p.median) <- source_names
    
    Cont1.plot <- Cont1.plot * 35.93976 + mix$CE_center  # transform Cont1.plot (x-axis) back to the original scale
    df <- data.frame(reshape2::melt(p.median)[, 2:3], rep(Cont1.plot, n.sources), 
        reshape2::melt(p.low)[, 3], reshape2::melt(p.high)[, 3])
    colnames(df) <- c("source", "median", "x", "low", "high")
    
    # Plot of Diet vs. Cont effect Fill##
    
    df$spp <- ifelse(df$high < 100, sppx[x])
    
    data_box = rbind(data_box, df)
}
################### END MCMC LOOP###
```

``` r
data_box$spp <- as.factor(data_box$spp)
data_box$spp <- factor(data_box$spp, levels=c("LAAL", "BUPE", "WTSH", "WTTR", "BRBO","BRNO","WHTE","SOTE","TP"))

#data_box <- subset(data_box, spp == c("LAAL", "BUPE", "WTSH", "WTTR", "BRBO","BRNO","WHTE","SOTE"))

ggplot(data_box, aes(x=x))+geom_line(aes(y=median,color = source),size = 1)+scale_color_brewer(type = "seq",palette = 'Spectral' )+facet_wrap(~spp,ncol =2)+theme_classic()+ 
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.5),
        legend.title=element_blank(),
        strip.text=element_text(hjust=0))+
  xlab("year")+
  ylab("proportion of diet")
```

<img src="Fig3_mixingmodel_files/figure-markdown_github/unnamed-chunk-3-1.png" width="672" />