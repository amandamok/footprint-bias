rm(list=ls())

library(here)
library(patchwork)

# load codon properties [[ origin ?? ]]
codon_properties <- read.delim(file.path(here(), "reference_data", "yeast_codon_properties.txt"),
                               stringsAsFactors=F)

# add regression coefficients
plot_data <- lapply(c("green", "lareau", "weinberg"),
                    function(expt) {
                      expt_dir <- paste0(expt, "_20cds20_f5_2_f3_3")
                      # load model fit
                      model_obj <- paste0(expt, "_fit_150")
                      load(file.path(here(), "expts", expt_dir, paste0(model_obj, ".Rda")))
                      # pull coefficients from model
                      model_coef <- data.frame(summary(get(model_obj))$coefficients)
                      model_coef <- model_coef[match(paste0("A", codon_properties$codon),
                                                     rownames(model_coef)),]
                      # compute fold-change and upper/lower bounds
                      model_coef$FC <- with(model_coef, exp(Estimate))
                      model_coef$FC_lower <- with(model_coef, exp(Estimate - Std..Error))
                      model_coef$FC_upper <- with(model_coef, exp(Estimate + Std..Error))
                      # return codon_properties with added coefficients
                      return(cbind(codon_properties, model_coef, expt=expt))
                    })
plot_data <- do.call(rbind, plot_data)
levels(plot_data$expt) <- c("Green", "Lareau", "Weinberg")
plot_data <- plot_data[!(plot_data$codon %in% c("TAA", "TAG", "TGA", "AAA")),]

corr_text <- data.frame(expt = levels(plot_data$expt),
                        corr = sapply(levels(plot_data$expt),
                                      function(x) {
                                        paste("cor = ",
                                              round(with(subset(plot_data, expt==x),
                                                         cor(tAI, FC)), 2))
                                      }),
                        pvalue = sapply(levels(plot_data$expt),
                                        function(x) {
                                          paste("p =",
                                                signif(with(subset(plot_data, expt==x),
                                                            signif(cor.test(FC, tAI)$p.value, 2))))
                                        }))

full_plot <- ggplot(plot_data, aes(x=tAI, y=FC)) + geom_point() + facet_grid(~expt) +
  geom_errorbar(aes(ymin=FC_lower, ymax=FC_upper)) + geom_smooth(method="lm", formula="y~x") +
  geom_text(data=corr_text, mapping=aes(x=0.85, y=0.95*max(plot_data$FC_upper), label=corr)) +
  geom_text(data=corr_text, mapping=aes(x=0.85, y=0.9*max(plot_data$FC_upper), label=pvalue)) +
  theme_classic() + theme(panel.background=element_rect(fill=NA, color="black")) +
  xlab("codon adaptiveness") + ylab(expression("exp("*beta*")")) +
  ggtitle("Regression model with library preparation terms (d5, d3, f5, f3)")

# regression model without library prep terms -----------------------------

alt_model <- formula(count ~ transcript + A + P + E)

alt_plot_data <- lapply(c("green", "lareau", "weinberg"),
                        function(expt) {
                          expt_dir <- file.path(here(), "expts", paste0(expt, "_20cds20_f5_2_f3_3"))
                          # load alternative model fit
                          model_obj <- paste0(expt, "_alt_fit_150")
                          if(!file.exists(file.path(expt_dir, paste0(model_obj, ".Rda")))) {
                            # load training data
                            load(file.path(here(), "expts", paste0(expt, "_20cds20_f5_2_f3_3"),
                                           paste0(expt, "_training_data.Rda")))
                            # compute fit for alternative model
                            assign(model_obj,
                                   MASS::glm.nb(alt_model, data=get(paste0(expt, "_training")), model=F))
                            # save fit
                            save(list=model_obj, file=file.path(expt_dir, paste0(model_obj, ".Rda")))
                          } else {
                            load(file.path(expt_dir, paste0(model_obj, ".Rda")))
                          }
                          # pull coefficients from model
                          model_coef <- data.frame(summary(get(model_obj))$coefficients)
                          model_coef <- model_coef[match(paste0("A", codon_properties$codon),
                                                         rownames(model_coef)),]
                          # compute fold-change and upper/lower bounds
                          model_coef$FC <- with(model_coef, exp(Estimate))
                          model_coef$FC_lower <- with(model_coef, exp(Estimate - Std..Error))
                          model_coef$FC_upper <- with(model_coef, exp(Estimate + Std..Error))
                          # return codon_properties with added coefficients
                          return(cbind(codon_properties, model_coef, expt=expt))
                        })
alt_plot_data <- do.call(rbind, alt_plot_data)
levels(alt_plot_data$expt) <- c("Green", "Lareau", "Weinberg")
alt_plot_data <- alt_plot_data[!(alt_plot_data$codon %in% c("TAA", "TAG", "TGA", "AAA")),]

alt_corr_text <- data.frame(expt = levels(alt_plot_data$expt),
                            corr = sapply(levels(alt_plot_data$expt),
                                          function(x) {
                                            paste("cor = ",
                                                  round(with(subset(alt_plot_data, expt==x),
                                                             cor(tAI, FC)), 2))
                                          }),
                            pvalue = sapply(levels(alt_plot_data$expt),
                                            function(x) {
                                              paste("p =",
                                                    signif(with(subset(alt_plot_data, expt==x),
                                                                signif(cor.test(FC, tAI)$p.value, 2))))
                                            }))

alt_plot <- ggplot(alt_plot_data, aes(x=tAI, y=FC)) + geom_point() + facet_grid(~expt) +
  geom_errorbar(aes(ymin=FC_lower, ymax=FC_upper)) + geom_smooth(method="lm", formula="y~x") +
  geom_text(data=alt_corr_text, mapping=aes(x=0.85, y=0.95*max(alt_plot_data$FC_upper), label=corr)) +
  geom_text(data=alt_corr_text, mapping=aes(x=0.85, y=0.9*max(alt_plot_data$FC_upper), label=pvalue)) +
  theme_classic() + theme(panel.background=element_rect(fill=NA, color="black")) +
  xlab("codon adaptiveness") + ylab(expression("exp("*beta*")")) +
  ggtitle("Regression model without library preparation terms")

full_plot / alt_plot
