## loading up necessary packages
library(tidyverse)
library(ggplot2)
library(cowplot)
library(lme4)
library(lmerTest)
library(gridExtra)
library(logistf)
library(bridgesampling)
library(brms)
library(here)
library(cmdstanr)
library(bayesplot)

setwd("G:/My Drive/Everything/Belgeler/RDirectory")

#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#cmdstanr::check_cmdstan_toolchain(fix = TRUE)
#install_cmdstan()


## my little function that takes logits and returns probabilities
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1+odds)
  return(prob)
}


# Processing... -----------------------------------------------------------

# my analyses of the comprehension data
## read in data
setwd("G:/My Drive/Everything/Belgeler/RDirectory")
b_compdat <- read.csv("G:/My Drive/Everything/Belgeler/RDirectory/comprehension_data_processedv3_pnoandagefixed.csv")

## processing

## the first 5, the 79th (break) and the last trial template in each persons experiment are non-experimental trials
b_non_explist <- c(unique(b_compdat$trial_template)[1:5],unique(b_compdat$trial_template)[length(unique(b_compdat$trial_template))],unique(b_compdat$trial_template)[79])
## all the others in-between are experimental (144) or catch trials (11)
b_explist <- c(unique(b_compdat$trial_template)[6:78],unique(b_compdat$trial_template)[80:(length(unique(b_compdat$trial_template))-1)])
## defining experimental and non-experimental trials
b_experimental_trials <- filter(b_compdat,trial_template %in% b_explist)
b_nonexperimental_trials <- filter(b_compdat,trial_template %in% b_non_explist)

## processing to give word order, case, sentence type, scenario, and quarter information
b_infolist <- strsplit(b_experimental_trials$stimuli_presented,"_")

b_case_list <- c()
b_wo_list <- c()
b_quart_list <- c()
b_scenario_list <- c()
b_sentence_type_list <- c()

for(i in 1:length(b_infolist)) {
  temp_wo_list <- substr(b_infolist[[i]][1], 1, 3)
  temp_case_list <- substr(b_infolist[[i]][1], 4, 4)
  temp_quart_list <- b_infolist[[i]][5]
  temp_scenario_list  =  paste(b_infolist[[i]][2], b_infolist[[i]][3], b_infolist[[i]][4],sep="")
  temp_sentence_type_list <- substr(b_infolist[[i]][1], 1, 4)
  b_scenario_list <- append(b_scenario_list,temp_scenario_list)
  b_case_list <- append(b_case_list,temp_case_list)
  b_wo_list <- append(b_wo_list,temp_wo_list)
  b_quart_list <- append(b_quart_list,temp_quart_list)
  b_sentence_type_list <- append(b_sentence_type_list,temp_sentence_type_list)
}

b_experimental_trials$scenario <- b_scenario_list
b_experimental_trials$case <- b_case_list
b_experimental_trials$wo <- b_wo_list
b_experimental_trials$quarter <- b_quart_list
b_experimental_trials$sentence_type <- b_sentence_type_list

b_catch_indexes <- which(substr(b_experimental_trials$trial_template,1,3)=="cat")
b_experimental_trials[b_catch_indexes,26] <- "catch"
b_experimental_trials[-b_catch_indexes,26] <- "experimental" 
colnames(b_experimental_trials)[26]<- "trial_type"

b_experimental_trials$scenario_no <- c(NA)
# Morpheme-level slice information
obj <- read.csv("G:/My Drive/Everything/Belgeler/RDirectory/elifdat.csv")
obj <- arrange(obj, obj$sentence)
objpr<-obj[-which(duplicated(obj$sentence)),] ## since every sentence has one signal associated, the 9 thirds are not needed, we only keep whatever was first from a given sentence, just to have a link to the signal. The gram role of the word that remain just depends on the sentence type and alphabetic, tells nothing theoretically.

b_experimental_trials$scenario_no <- NA
## Adding signals and scenario number as columns
for (i in 1:length(b_experimental_trials$scenario)) {
  index = which(b_experimental_trials$scenario[i] == objpr$propcontent) ## objective signal comes from production data which had different headers than subjective signal headers where it comes from FindingFive (i think)
  b_experimental_trials$scenario_no[i] <- objpr[index,18] # 18th column is scenario number
}

b_experimental_trials$scenario_no <- as.numeric(b_experimental_trials$scenario_no)

b_experimental_trials[,28] <- NA
colnames(b_experimental_trials)[28]<- "question_type"

b_experimental_trials[which(b_experimental_trials$scenario_no %% 2 == 0),28] <- "active"
b_experimental_trials[which(b_experimental_trials$scenario_no %% 2 == 1),28] <- "passive"


b_proc_experimental_trials <- b_experimental_trials

## Removing participant 13 because non-native
## and 17 because at-chance performance on all conditions
## for making things cleaner
## and not selecting catch trials
b_proc_experimental_trials <- filter(b_experimental_trials, trial_type == "experimental")
b_proc_experimental_trials <- filter(b_proc_experimental_trials, participant_id != 13)
b_proc_experimental_trials <- filter(b_proc_experimental_trials, participant_id != 17)

b_a<-length(b_proc_experimental_trials$expt_id)

b_proc_experimental_trials <- filter(b_proc_experimental_trials, response_rt < 10000)
b_proc_experimental_trials <- filter(b_proc_experimental_trials, response_rt > 200)
b_b<-length(b_proc_experimental_trials$expt_id)
print(c("Lost Data is %",((b_a)-(b_b))/((b_a))))

b_proc_experimental_trials$pd <- NA
b_proc_experimental_trials[which(b_proc_experimental_trials$sentence_type=="sovm"),29] <- .1
b_proc_experimental_trials[which(b_proc_experimental_trials$sentence_type=="sovu"),29] <- .1
b_proc_experimental_trials[which(b_proc_experimental_trials$sentence_type=="ovsm"),29] <- 0
b_proc_experimental_trials[which(b_proc_experimental_trials$sentence_type=="ovsu"),29] <- 1.746


b_dat <- b_proc_experimental_trials

b_model.dat <- filter(b_dat)

## determining contrast coding schemes
b_model.dat$case <- as.factor(b_model.dat$case)
contrasts(b_model.dat$case) <- contr.sum(2)

b_model.dat$wo <- as.factor(b_model.dat$wo)
contrasts(b_model.dat$wo) <- contr.sum(2)

b_model.dat$question_type <- as.factor(b_model.dat$question_type)
contrasts(b_model.dat$question_type) <- contr.sum(2)

b_model.dat$sentence_type <- as.factor(b_model.dat$sentence_type)
contrasts(b_model.dat$sentence_type) <- contr.treatment(4)




# DDM ---------------------------------------------------------------------


# Formula of the model
b_model.dat

ddm_model_formula <- brms::bf(
  response_rt | dec(response_correct) ~ 0 + wo,
  bs ~ 0 + wo, 
  ndt ~ 0 + wo, 
  bias ~ 0 + wo
)

#T ake a look at priors
get_prior(ddm_model_formula, 
          family = wiener(
            link_bs = "identity", 
            link_ndt = "identity", 
            link_bias = "identity"), 
          data = b_model.dat
)

# Set weakly informative priors
prior <- c(
  prior("normal(0, 1)", class = "b"), # work for categorical factor, add more for more fixed factors in the model
  prior("normal(0, 5)", class = "b", dpar = "bs"),
  prior("normal(0.2, 1)", class = "b", dpar = "ndt"),
  prior("normal(0.5, 1)", class = "b", dpar = "bias")
)

# Print stan code to check specification
make_stancode(
  ddm_model_formula, 
  family = wiener(
    link_bs = "identity", 
    link_ndt = "identity", 
    link_bias = "identity"),
  data = b_model.dat, 
  prior = prior)

# Generate temp data for possible init values
tmp_dat <- make_standata(
  formula = ddm_model_formula, 
  family = wiener(
    link_bs = "identity", 
    link_ndt = "identity", 
    link_bias = "identity"), 
  data = b_model.dat, 
  prior = prior)

str(tmp_dat, 1, give.attr = FALSE)

# Create init function to be used when fitting model
initfun <- function() {
  list(
    b = rnorm(tmp_dat$K),
    b_bs = runif(tmp_dat$K_bs, 0, 5), #1, 2
    b_ndt = runif(tmp_dat$K_ndt, 0, 5), 
    b_bias = rnorm(tmp_dat$K_bias, .5, 0.1) 
  )
}

# Test model on single participant
test <- brm(
  formula = ddm_model_formula,
  data = filter(b_model.dat, participant_id == 23),
  prior = prior,
  sample_prior = TRUE,
  init = initfun,
  family = wiener(
    link_bs = "identity", 
    link_ndt = "identity",
    link_bias = "identity"),
  backend = "cmdstanr", 
  chains = 4, warmup = 10000, iter = 15000, thin = 10, 
  cores = 4, 
  control = list(adapt_delta = 0.95, max_treedepth = 1000000000), 
  file = "G:/My Drive/Everything/Belgeler/RDirectory/comp1/ddms"
)
