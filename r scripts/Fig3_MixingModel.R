library(data.table)
library(MixSIAR)
library(dplyr)
library(plyr)
library(ggplot2)
library(R2jags)
####generate LOESS timelines of TP by species and total
setwd('/Users/tgagne/Documents/Spring 2017/Seabird TP/Data/CSV')
data<-read.csv('seabird_CSSIA_NPacific_1890_2016_rawdata.csv', header = T)
#generate dataset from which to draw
data <- droplevels( data[-which(data$spp == "RFBO"), ] )
dataList = split(data, data$uc_id)
set.seed(1234)
##################
###DO NOT RERUN###
##################
ID_1 <-dataList[[1]]
GetTP <- function(anID){
  b = 2.4167
  TEF = 5.6333
  ala.est<- rnorm(1000, mean = as.numeric(subset(anID, value == "ave", select = "ala")) , sd = as.numeric(subset(anID, value == "sd", select = "ala")))
  glu.est<- rnorm(1000, mean = as.numeric(subset(anID, value == "ave", select = "glue")) , sd = as.numeric(subset(anID, value == "sd", select = "glue")))
  ile.est<- rnorm(1000, mean = as.numeric(subset(anID, value == "ave", select = "ile")) , sd = as.numeric(subset(anID, value == "sd", select = "ile")))
  leu.est<- rnorm(1000, mean = as.numeric(subset(anID, value == "ave", select = "leu")) , sd = as.numeric(subset(anID, value == "sd", select = "leu")))
  pro.est<- rnorm(1000, mean = as.numeric(subset(anID, value == "ave", select = "pro")) , sd = as.numeric(subset(anID, value == "sd", select = "pro")))
  val.est<- rnorm(1000, mean = as.numeric(subset(anID, value == "ave", select = "val")) , sd = as.numeric(subset(anID, value == "sd", select = "val")))
  phe.est<- rnorm(1000, mean = as.numeric(subset(anID, value == "ave", select = "phe")) , sd = as.numeric(subset(anID, value == "sd", select = "phe")))
  TPforID <- (((((ala.est + glu.est + ile.est + leu.est + pro.est + val.est)/6)) - phe.est - b)/TEF ) + 1 
  return(TPforID)
}
GetTP(ID_1)
TPforID_all = lapply(dataList, GetTP ) # apply function to list
glimpse(TPforID_all)
TPforID_all = melt(TPforID_all)
TPforID_all$L1 <- as.factor(TPforID_all$L1)
setnames(TPforID_all, old=c("L1","value"), new=c("uc_id", "tp"))
str(TPforID_all)
data<-subset(data, value == "ave")
str(data)
total<-merge(data,TPforID_all, by="uc_id")
str(total)
total <- subset(bigger_data, select = c(uc_id,
                                            spp,
                                            year,
                                            tp))
                                    
####Building and exeute the loop
sppx<-c("SOTE","WTSH","BRBO","BRNO","LAAL","BUPE","WTTR","WHTE","TP")
d = NULL
for(x in 1:length(sppx)){
  if(x < 9) {
    species_data <- subset(total, spp == sppx[x])
  } else{
    species_data <- total
  }
tp_predict = as.data.frame(matrix(ncol=1000,nrow=124))
for (i in 1:1000){
  new_df <- ddply(species_data,.(year),function(x) x[sample(nrow(x),1),])
  tp_est <- loess(tp ~ year, new_df,span = 1)
  tp_predict[i] <- as.data.frame(predict(tp_est,data.frame(year = seq(1891, 2014,1))))
}
tp_predict$year <- seq(1891, 2014, by = 1)
test_data_long <- melt(tp_predict, id="year")
sum_data<- do.call(data.frame,aggregate(value ~ year, data = test_data_long, function(x) c(quantile(x,.841),median(x) )))
sum_data$sd <- sum_data$value.84.1. - sum_data$value.V2
sum_data$value.84.1. <- NULL
setnames(sum_data, "value.V2", "TL")
sum_data$spp <- ifelse(sum_data$year < 3000, sppx[x])
d = rbind(d,sum_data)

}
##################
#####End##########
##################

setwd('/Users/tgagne/Documents/Spring 2017/Mixing Model Run')
#write.csv(d,file = "mixture_data_may22.csv")
################################
##Building absolute loss table##
################################
absolute<-read.csv('/Users/tgagne/Documents/Spring 2017/Mixing Model Run/mixture_data_may22.csv')
ggplot(absolute,aes(x=year,y=TL,color=spp))+
  geom_line(size=0.5)+theme_classic()
df<-subset(absolute,year == 1891 | year == 1950 | year == 2014)
df$sd<-NULL
df$X <- NULL
ggplot(df,aes(x=year))+geom_line(aes(y=TL,color=spp))+geom_point(aes(y=TL,color=spp),shape = 21,fill="white")+theme_classic()+ylab("trophic position")
###############################################################
##Big MC/MC loop to run and generate plots for each species##
###############################################################
sppx<-c("SOTE","WTSH","BRBO","BRNO","LAAL","BUPE","WTTR","WHTE","TP")
data_box = NULL
for(x in 1:length(sppx)){
mix<-read.csv('/Users/tgagne/Documents/Spring 2017/Mixing Model Run/mixture_data_may22.csv')
str(mix)
mix<-subset(mix, spp == sppx[x])
write.csv(mix,file = "mixture_data_spp.csv")
#mix data = i.e. consumer
mix<- load_mix_data(filename = "mixture_data_spp.csv",iso_names = "TL",factors = "spp",
                          fac_random = TRUE,fac_nested = NULL,cont_effects = "year")
#source data #new omma and caran TL
source <- load_source_data( filename = "source_data_4grp_may22.csv",
                            source_factors = NULL,conc_dep = FALSE,data_type = "means",mix)
#discrimination data
discr <- load_discr_data(filename = "discrimination_data_4grp.csv",mix)
#one dimensional isospace plot
#plot_data(filename="isospace_plot",
 #     plot_save_pdf=FALSE,
  #    plot_save_png=FALSE,
   #   mix,source,discr)
#uninformative prior = alpha.prior = 1
#construct informative prior from stomach_proportional_plot_apr13.xlsx
TP_prior <- c(0.508,1.212,1.175,1.105)
#plot_prior(alpha.prior = TP_prior, source,
 #         plot_save_pdf=FALSE,
  #       plot_save_png=FALSE)
#write jags model
model_filename <- "MixSIAR_model.txt"
resid_err <- FALSE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err,process_err, mix, source)
#run model
jags.1 <- run_model(run="fast",mix,source,discr,model_filename,
                    alpha.prior = TP_prior,resid_err,process_err)
##########################
##########################
R2jags::attach.jags(jags.1)
n.sources <- source$n.sources
source_names <- source$source_names

      fac.lab <- mix$FAC[[1]]$labels
      label <- mix$cont_effects
      cont <- mix$CE[[1]] #either CE or CE_orig explore implications?
      ilr.cont <- get(paste("ilr.cont",1,sep=""))
      
      get_high <- function(x){return(quantile(x,.95))}
      get_low <- function(x){return(quantile(x,.05))}
      
      n.plot = 124 #200 was original consider changing to number of years? which I did
      chain.len = dim(p.global)[1]
      Cont1.plot <- seq(from=round(min(cont),1), to=round(max(cont),1), length.out=n.plot)
      ilr.plot <- array(NA,dim=c(n.plot, n.sources-1, chain.len))
      ilr.median <- array(NA,dim=c(n.plot, n.sources-1))
      ilr.low <- array(NA,dim=c(n.plot, n.sources-1))
      ilr.high <- array(NA,dim=c(n.plot, n.sources-1))
      for(src in 1:n.sources-1){
        for(i in 1:n.plot){
          ilr.plot[i,src,] <- ilr.global[,src] + ilr.cont[,src]*Cont1.plot[i] #+ ilr.fac1[,f1,src]
          ilr.low[i,src] <- get_low(ilr.plot[i,src,])
          ilr.median[i,src] <- median(ilr.plot[i,src,])  #changed from median to mean 
          ilr.high[i,src] <- get_high(ilr.plot[i,src,])
        }
      }
      
      # Transform regression lines from ILR-space to p-space
      e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
      for(i in 1:(n.sources-1)){
        e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
        e[,i] <- e[,i]/sum(e[,i])
      }
      cross.med <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
      tmp.p.med <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
      p.median <- array(data=NA,dim=c(n.plot, n.sources))
      cross.low <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
      tmp.p.low <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
      p.low <- array(data=NA,dim=c(n.plot, n.sources))
      cross.high <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
      tmp.p.high <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
      p.high <- array(data=NA,dim=c(n.plot, n.sources))
      for(i in 1:n.plot){
        for(j in 1:(n.sources-1)){
          cross.med[i,,j] <- (e[,j]^ilr.median[i,j])/sum(e[,j]^ilr.median[i,j]);
          cross.low[i,,j] <- (e[,j]^ilr.low[i,j])/sum(e[,j]^ilr.low[i,j]);
          cross.high[i,,j] <- (e[,j]^ilr.high[i,j])/sum(e[,j]^ilr.high[i,j]);
        }
        for(src in 1:n.sources){
          tmp.p.med[i,src] <- prod(cross.med[i,src,]);
          tmp.p.low[i,src] <- prod(cross.low[i,src,]);
          tmp.p.high[i,src] <- prod(cross.high[i,src,]);
        }
        for(src in 1:n.sources){
          p.median[i,src] <- tmp.p.med[i,src]/sum(tmp.p.med[i,]);
          p.low[i,src] <- tmp.p.low[i,src]/sum(tmp.p.low[i,]);
          p.high[i,src] <- tmp.p.high[i,src]/sum(tmp.p.high[i,]);
        }
      }
      colnames(p.median) <- source_names
    
      Cont1.plot <-  Cont1.plot*35.93976 + mix$CE_center # transform Cont1.plot (x-axis) back to the original scale
      df <- data.frame(reshape2::melt(p.median)[,2:3], rep(Cont1.plot,n.sources), reshape2::melt(p.low)[,3], reshape2::melt(p.high)[,3])
      colnames(df) <- c("source","median","x","low","high")
      
      # Plot of Diet vs. Cont effect
      #dev.new()
      print(ggplot2::ggplot(data=df,ggplot2::aes(x=x,y=median)) +
              ggplot2::geom_line(ggplot2::aes(x=x, y=median,group=source,colour=source),size=1.5) +
              #ggplot2::geom_ribbon(ggplot2::aes(ymin=low, ymax=high, group=source, fill=source), alpha=0.35) +
              ggplot2::labs(title = fac.lab) +
              ggplot2::ylab("Diet Proportion") +
              ggplot2::xlab(label) +
              ggplot2::scale_y_continuous(expand = c(0, 0), limits=c(0,1)) +
              ggplot2::theme_bw() +
              ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), 
                             panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                             axis.line = ggplot2::element_line(colour = "black"), axis.title=ggplot2::element_text(size=16), 
                             axis.text=ggplot2::element_text(size=14), legend.text=ggplot2::element_text(size=14), legend.position=c(.02,1), 
                             legend.justification=c(0,1), legend.title=ggplot2::element_blank()))

      df$spp <- ifelse(df$high < 100, sppx[x])
      
      data_box = rbind(data_box,df)
}
###################
###END MCMC LOOP###
###################
#write.csv(data_box,file = "MCMCmix_runMay22.csv")
data_box<-read.csv("MCMCmix_runMay21.csv",header = T)

data_box$spp <- as.factor(data_box$spp)
data_box$SE <- (abs(data_box$high-data_box$low)/4)
data_box$spp <- factor(data_box$spp, levels=c("LAAL", "BUPE", "WTSH", "WTTR", "BRBO","BRNO","WHTE","SOTE","TP"))
#data_box$source <- factor(data_box$source, levels=c("carangidae", "ommastrephidae", "exocoetidae", "mullidae", "myctophidae","sternoptychidae"))
#Facetted plot of distribution of diet
ggplot(data_box, aes(x=x))+geom_line(aes(y=median,color = source),size = 1)+scale_color_brewer(type = "seq",palette = 'Spectral' )+facet_wrap(~spp,ncol =2)+theme_classic()

##################################
#Attempt to build a mean ensemble#
##################################
ensemble <- data_box
ensemble<-ensemble[!ensemble$spp %in% "TP",]
ensemble<-droplevels(ensemble)
colnames(ensemble)[4] <- "year"
ensemble$spp <- NULL
ensemble$X <- NULL
aggdata <- aggregate(ensemble$median, by=list(ensemble$year,ensemble$source),FUN=mean)
ggplot(aggdata,aes(x=Group.1))+geom_line(aes(y=x,color=Group.2),size=1)+theme_classic()+scale_color_brewer(type = "seq",palette = 'Spectral' )+scale_y_continuous(limits = c(0,1))
##########################
##Trying to add confidence bands
##########################

df<-data_box
#df<-subset(df, spp == c("LAAL", "BUPE", "WTSH", "WTTR", "BRBO","BRNO","WHTE","SOTE") )
#df$spp <- droplevels(df$spp)
#rough range rule estimated SE error bars
df$y_lower <- ifelse(df$median-2*df$SE>0,df$median-2*df$SE,0)
df$y_upper <- ifelse(df$median+2*df$SE<1,df$median+2*df$SE,1)



ggplot(df,aes(x=x))+
  geom_line(aes(y=median,color = source))+
  geom_ribbon(aes(ymin=y_lower,ymax=y_upper,fill=source),alpha = .3)+
  scale_color_brewer(type = "seq",palette = 'Spectral' )+
  facet_wrap(~spp,ncol=2)+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.5),
        legend.key.size = unit(.75, "cm"),
        legend.title=element_blank(),
        strip.text=element_text(hjust=0))+
  ylab("trophic position")+
  scale_x_continuous(limits=c(1891,2014))

####################################################
##########################Breaking in to one iso plot
####################################################
plot_data_one_iso_test <- function(mix,source,discr,filename,plot_save_pdf,plot_save_png){

  if(!exists("x_label")) x_label <- mix$iso_names
  
  y_data <- 0.5
  y <- rep(y_data,mix$N)
  df <- data.frame(x = mix$data_iso, y = y)
  spacing <- 0.1
  
  if(!is.na(source$by_factor)){} else {  # source$by_factor = NA
    source_linetype <- 1:source$n.sources    # each source gets a different linetype
    source_color <- factor(rep("black",source$n.sources))  # this doesn't work...solution was to make separate ggplot calls for by_factor and not_by_factor
    index <- 1:source$n.sources              # "index" gets the row in S_MU of the first instance of each source (since not by factor, only one instance of each source)
    discr_mu_plot <- discr$mu
    discr_sig2_plot <- discr$sig2
    y_sources <- seq(y_data+0.2,(source$n.sources*spacing)-spacing+y_data+0.2,by=spacing)
    source$S_factor_levels <- 0.5   # just for correctly spacing the source labels (y in source.labels, line 61)
    MU_plot <- source$S_MU + discr_mu_plot    # add discrimination mean to the source mean values
    SIG_plot <- sqrt(source$S_SIG^2 + discr_sig2_plot)  # add discrimination SD to the source SD values
    }
  MU_plot <- as.vector(MU_plot)
  SIG_plot <- as.vector(SIG_plot)
  df_sources <- data.frame(x=MU_plot,
                           xmin = MU_plot - SIG_plot,
                           xmax = MU_plot + SIG_plot,
                           y = y_sources,
                           linetype = source_linetype,
                           scolour = source_color)
  source.labels <- data.frame(
    x = MU_plot[index],    # label sources just left
    #y = MU_plot[index,2] + rep(1.5,n.sources),    # and up from their means
    y = y_sources[index] + spacing*source$S_factor_levels,
    label = source$source_names)
  .e <- environment()
  dev.new()
   # end n.effects==2
  if(mix$n.effects==1){
    if(!is.na(source$by_factor)){ } else { # sources not by factor (make the sources black)
      g <- ggplot2::ggplot(data = df,ggplot2::aes(x = x,y = y),environment=.e) +
        ggplot2::geom_point(ggplot2::aes(colour = factor(mix$FAC[[1]]$values)), position=position_jitter(width=.2,height=.1), show.legend=T) +  # Factor.1
        #ggplot2::scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)),  # Factor.1
        #                               labels = mix$FAC[[1]]$labels) +     # factor1_names
        ggplot2::geom_point(data=df_sources,
                            ggplot2::aes(x = x, y = y),
                            size=2,
                            show.legend=F) +
        ggplot2::geom_errorbarh(data=df_sources,
                                ggplot2::aes(xmin=xmin,xmax=xmax),
                                size=1,
                                height=0,
                                linetype=source_linetype,
                                show.legend=F) +
        ggplot2::geom_text(data=source.labels, ggplot2::aes(x=x,y=y,label=label), show.legend=F) +
        ggplot2::scale_y_continuous(breaks = NULL) +
        ggplot2::ylab("") +
        ggplot2::xlab(x_label) +
        ggplot2::theme_classic() +
        ggplot2::scale_color_brewer(breaks = levels(factor(mix$FAC[[1]]$values)),labels = mix$FAC[[1]]$labels,type = "seq",palette = 'RdYlBu' ) + 
        ggplot2::theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=ggplot2::element_blank())
      print(g)
    } } }



mix<- load_mix_data(filename = "mixture_data_may22.csv" ,iso_names = "TL",factors = "spp",
                    fac_random = TRUE,fac_nested = NULL,cont_effects = "year")
plot_data_one_iso_test(filename="isospace_plot",mix,source,discr)
