library(here)

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

ggplot(plot_data, aes(x=tAI, y=FC)) + geom_point() + facet_grid(~expt) +
  geom_errorbar(aes(ymin=FC_lower, ymax=FC_upper)) + geom_smooth(method="lm", formula="y~x") +
  geom_text(data=corr_text, mapping=aes(x=0.85, y=0.95*max(plot_data$FC_upper), label=corr)) +
  geom_text(data=corr_text, mapping=aes(x=0.85, y=0.9*max(plot_data$FC_upper), label=pvalue)) +
  theme_classic() + theme(panel.background=element_rect(fill=NA, color="black")) +
  xlab("codon adaptiveness") + ylab(expression("exp("*beta*")"))
