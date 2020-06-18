simulate_simple_mediation <- function(n,p_otu, p_metabolite, mediations = 1){
  not_mediation_otu <- p_otu - 2*mediations
  not_mediation_metabolite <- p_metabolite - mediations
  otutable <- t(replicate(not_mediation_otu, runif(n)))
  metabolitetable <- t(replicate(not_mediation_metabolite, runif(n)))
  for (m in 1:mediations){
    temp = simple_mediation_direct(n)
    otutable <- rbind(otutable, temp$OTU)
    metabolitetable <- rbind(metabolitetable, temp$Metabolite)
  }
  final <-list(OTU = otutable, Metabolite = metabolitetable,
               OtuAnn = data.frame(Species = paste0("OTU", 1:p_otu), stringsAsFactors = F),
               MetAnn = data.frame(Species = paste0("M", 1:p_metabolite), stringsAsFactors = F))
  final$MetAnn$Species[not_mediation_metabolite:length(final$MetAnn$Species)] <- paste(final$MetAnn$Species[not_mediation_metabolite:length(final$MetAnn$Species)], "Mediator")
  return(final)
  }

simple_mediation <- function(n, beta11 = -2, beta21 = 2, b32 = -1.5){
  producer = runif(n)
  mm = beta21*producer + rnorm(n)
  target = beta11*producer + b32*mm + rnorm(n)
  return(list(OTU = rbind(producer, target), Metabolite = mm))
}

simple_mediation_direct <- function(n, beta11 = -2, beta21 = 2, b32 = -1.5){
  producer = runif(n)
  mm = beta21*producer + rnorm(n)
  target =  b32*mm + rnorm(n)
  return(list(OTU = rbind(producer, target), Metabolite = mm))
}