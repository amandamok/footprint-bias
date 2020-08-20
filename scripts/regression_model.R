## requires libraries: MASS

##### ASSUMES correct_bias() IN GLOBAL ENVIRONMENT
##### NEED TO source("~/footprint-bias/scripts/correct_bias.R")

model_frame_by_size <- function(regression_data, min_prop = 0.9,
                                model = formula(count ~ transcript + A + P + E + f5 + f3)) {
  # compute regression per frame/size pair; return corrected counts for subsets
  ## regression_data: data.frame; output from init_data() and count_footprints()
  ## min_prop: numeric; minimum proportion of footprints to model
  ## model: formula for MASS::glm.nb()
  # 1. establish subsets to compute regression over
  subsets <- aggregate(count ~ d5 + d3, data=regression_data, FUN=sum)
  subsets <- subsets[order(subsets$count, decreasing=T),]
  subsets$prop <- sapply(seq(nrow(subsets)), function(x) sum(subsets$count[1:x]))/sum(subsets$count)
  num_subsets <- min(which(subsets$prop >= min_prop))
  # 2. subset data
  subset_data <- lapply(seq(num_subsets),
                        function(x) {
                          subset(regression_data, d5 == subsets$d5[x] & d3 == subsets$d3[x])
                        })
  # 3. compute regression
  subset_models <- lapply(subset_data, function(x) { MASS::glm.nb(model, data=x, model=F) })
  # 4. pull regression coefficients
  subset_coefs <- sapply(subset_models, function(x) coef(x))
  colnames(subset_coefs) <- sapply(seq(num_subsets),
                                   function(x) {
                                     paste("d5", subsets$d5[x], "d3", subsets$d3[x], sep="_")
                                   })
  # 5. correct biases
  subset_data <- do.call(rbind, 
                         lapply(seq(num_subsets),
                                function(x) {
                                  correct_bias(subset_data[[x]], subset_models[[x]])
                                }))
  return(list(corrected_data=subset_data, regression_coefs=subset_coefs))
}
