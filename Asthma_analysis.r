library(dplyr)
library(ggplot2)
library(data.table)
library(survival)
library(foreign)
library(survminer)
asthma <- read.table("..\\asthma-prevention-trial.txt",header=TRUE,sep = "")
old_names <- c("id.w"  ,  "trt.w"  , "start.w", "stop.w"  ,"st.w"   , "nn"   ,   "fevent" )
new_names <- c( "id","treatment","start","stop","event","no.events","first event")
setnames(asthma, old = old_names, new = new_names)
freq <- table(asthma$treatment,asthma$id)
freq
asthma <- asthma %>% group_by(id) %>% mutate (asthma_count = row_number())
asthma$treatmentn[asthma$treatment == 0] = 1
asthma$treatmentn[asthma$treatment == 1] = 0
asthma$treatmentc[asthma$treatmentn == 0] = "placebo"
asthma$treatmentc[asthma$treatmentn == 1] = "medication"
asthma$eventc[asthma$event == 0] = "censored"
asthma$eventc[asthma$event == 1] = "event"

asthma_censored<- asthma %>% filter(eventc=="censored")
asthma2_censored <- asthma_censored %>% filter(id %in% c("3","8","9", '11',"13",'14',"16","31","33","43")) 
asthma2<- asthma%>% filter(id %in% c("3","8","9", '11',"13",'14',"16","31","33","43")) 
plot_t1<-asthma2 %>% ggplot(aes(c(0,max(stop)), as.factor(id))) + geom_segment(aes(x = 0,xend =stop, y=as.factor(id), yend = as.factor(id))) + geom_point(aes(x = stop, y= as.factor(id),colour=factor(eventc))) +
  facet_grid(.~treatmentc)
plot_t1
plot_t2<- plot_t1 +  geom_point(aes(x = stop, y = as.factor(id),colour = factor(eventc)), data = asthma2_censored)  +
  labs(x="follow up time in days",y="patient id", col ="Event")+ 
  ggtitle("recurrent event history") + 
  theme_bw()
plot_t2

#Subject Trial

#Total Time
sub1 <- asthma2 %>% filter(id %in% c("13"))
plot_sub1<-sub1 %>% ggplot(aes(c(0,max(stop)), as.factor(asthma_count))) + geom_segment(aes(x = 0,xend =stop, y=as.factor(asthma_count), yend = as.factor(asthma_count))) + geom_point(aes(x = stop, y= as.factor(asthma_count),colour=factor(eventc)))
plot_sub1
plot_sub1<- plot_sub1+
  labs(x="follow-up (in days)",y="no. of asthma attack occurrence", col ="")+ 
  ggtitle("Total time") + 
  theme_bw()
plot_sub1


#Counting Process
sub1_2 <- asthma2 %>% filter(id %in% c("13"))
plot_sub1<-sub1_2 %>% ggplot(aes(c(0,max(stop)), as.factor(asthma_count))) + geom_segment(aes(x = start,xend =stop, y=as.factor(asthma_count), yend = as.factor(asthma_count))) + geom_point(aes(x = stop, y= as.factor(asthma_count),colour=factor(eventc)))
plot_sub1
plot_sub1<- plot_sub1+
  labs(x="follow-up time (in days)",y="no. of asthma attack occurrence", col ="")+ 
  ggtitle("Counting process") + 
  theme_bw()
plot_sub1

#Gap Time 

sub_gap <- asthma2 %>% filter(id %in% c("13")) %>% mutate(gap_time = stop - start)

gap_plot<-sub_gap %>% ggplot(aes(c(0,max(gap_time)), as.factor(asthma_count))) + geom_segment(aes(x = 0,xend =gap_time, y=as.factor(asthma_count), yend = as.factor(asthma_count))) + geom_point(aes(x = gap_time, y= as.factor(asthma_count),colour=factor(eventc)))
gap_plot
gap_plot<- gap_plot+
  labs(x="follow-up time (in days)",y="no. of asthma attack occurrence", col ="")+ 
  ggtitle("Gap time") + 
  theme_bw()
gap_plot

# counting process for all the subjects
plot_counting <- asthma2 %>% ggplot(aes(c(0,max(stop)), as.factor(id))) + geom_segment(aes(x = start,xend =stop, y=as.factor(id), yend = as.factor(id))) + geom_point(aes(x = stop, y= as.factor(id),colour=factor(eventc))) 
plot_counting

plot_counting_2<- plot_counting +  geom_point(aes(x = stop, y = as.factor(id),colour = factor(eventc)), data = asthma2_censored) + facet_wrap(~treatmentc) +
  labs(x="follow up period",y="'patient id'", col ="Event")+ 
  ggtitle("recurrent event history") + 
  theme_bw()
plot_counting_2


#AG model

AG_Model <- coxph (Surv (start,stop,event)~treatmentn + cluster (id), data = asthma)
summary (AG_Model)

AG_model_plot <- survfit(Surv(start, stop, event)~treatmentc + cluster (id), data=asthma, id=id, )

AG_plot <- ggsurvplot(AG_model_plot,
           conf.int = 0.95,
           risk.table.col = "treatmentc", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           fun = "cumhaz")
AG_plot


#-Total time Model:

PWP_TT <- coxph (Surv (start, stop,event) ~ treatmentn + cluster (id) + strata(asthma_count), data = asthma)
summary (PWP_TT)

TT_model_plot <- survfit(Surv(start, stop,event) ~ treatmentn + cluster (id) + strata(asthma_count), data = asthma )
TT_model_plot
TT_plot <- ggsurvplot(TT_model_plot,
                      conf.int = FALSE,
                      risk.table.col = "treatmentc", # Change risk table color by groups
                      ggtheme = theme_bw(), 
                      fun = "cumhaz")
TT_plot
print(TT_plot)

#PWP gap time model
PWP_GT <-coxph (Surv(stop-start,event) ~ treatmentn + cluster(id) + strata(asthma_count),data = asthma)
summary (PWP_GT)

GT_model_plot <- survfit(Surv(stop-start,event) ~ treatmentn + cluster(id) + strata(asthma_count),data = asthma)

GT_plot <- ggsurvplot(GT_model_plot,
                      conf.int = 0.95,
                      risk.table.col = "treatmentc", # Change risk table color by groups
                      ggtheme = theme_bw(), # Change ggplot2 theme
                      palette = c("#E7B800", "#2E9FDF"),
                      fun = "cumhaz")

GT_plot 

#Marginal Model

Marginal <- coxph (Surv (start,stop,event) ~ treatmentn + cluster(id) + strata(asthma_count),data = asthma)
summary(Marginal)

Marginal_model_plot <- survfit(Surv (start,stop,event) ~ treatment + cluster(id) + strata(asthma_count),data = asthma)

marginal_plot <- ggsurvplot(Marginal_model_plot,
                      conf.int = 0.95,
                      risk.table.col = "treatmentc", # Change risk table color by groups
                      ggtheme = theme_bw(), # Change ggplot2 theme
                      palette = c("#E7B800", "#2E9FDF"),
                      fun = "cumhaz")
marginal_plot


#Frailty model

Frailty <-coxph (Surv (start,stop,event) ~ treatmentn + frailty (id, dist = "gamma"),data = asthma)
summary (Frailty)

Frailty_model_plot <- survfit(Surv (start,stop,event) ~ treatment + frailty (id, dist = "gamma"),data = asthma)
Frailty_model_plot
Frailty_plot <- ggsurvplot(Frailty_model_plot,
                            conf.int = 0.95,
                            risk.table.col = "treatmentc", # Change risk table color by groups
                            ggtheme = theme_bw(),
                            fun = "cumhaz")

Frailty_plot
