library(sybil)
library(sybilcycleFreeFlux)
library(reshape2)
library(ggplot2)

sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1

#recon <- readRDS("~/uni/dat/mod/r/recon2v04.RDS")
recon <- readRDS("~/uni/dat/mod/r/recon2_2_corrected.RDS")
ex <- findExchReact(recon)

do_not_exclude <- c("biomass_reaction",
                    "R_group_phosphotase_1", # otherwise no growth possible, determined see below
                    "FAOXC4C2m",# beta oxidation
                    "PPCOACm", #reaction for enzyme 6.4.1.3
                    "MMMm",#reaction for enzyme 5.4.99.2
                    "CSm",#reaction for enzyme 2.3.3.1
                    "ICDHxm",#reaction for Isocitrate dehydrogenase (NAD+)
                    "r0384"#reaction for enzyme KEGG:1.2.4.2
)


# sigmoidal colon non inflammed
rea_cd_ni_sig <- scan("~/uni/dat/mod/tissue_models_201611/CD_ni_SIG-biomass.lst", sep=" ", what="")
rea_uc_ni_sig <- scan("~/uni/dat/mod/tissue_models_201611/UC_ni_SIG-biomass.lst", sep=" ", what="")
rea_dc_ni_sig <- scan("~/uni/dat/mod/tissue_models_201611/DC_ni_SIG-biomass.lst", sep=" ", what="")
subset <- union(rea_cd_ni_sig, union(rea_uc_ni_sig, rea_dc_ni_sig))
unused <- recon@react_id[which(!recon@react_id %in% subset & !recon@react_id %in% ex@react_id)]
tissue <- rmReact(recon, react=setdiff(unused, do_not_exclude))
mod_id <- "recon_ni_sig"; tissue@mod_desc <- mod_id; tissue@mod_id <- mod_id



system.time(
  fluxes <- (ACHR(tissue,react_num(tissue)*4,5000))$Points)
rownames(fluxes) <- tissue@react_id

mlt <- melt(fluxes)
ggplot(mlt, aes(Var1, value)) + geom_boxplot()



#
# check recon
#

recon <- readRDS("~/uni/dat/mod/r/recon2v04.RDS") # minimal medium is not working with recon 2.04
# only for recon 2.04
minimal <- c("EX_glc_D(e)", 
             "EX_lnlc(e)", "EX_lnlnca(e)",
             "EX_his_L(e)", "EX_ile_L(e)", "EX_leu_L(e)", "EX_lys_L(e)", "EX_met_L(e)", "EX_phe_L(e)", "EX_thr_L(e)", "EX_trp_L(e)", "EX_val_L(e)",
             "EX_pnto_R(e)", "EX_ribflv(e)",
             "EX_nh4(e)", "EX_o2(e)", "EX_pi(e)", "EX_so4(e)", "EX_ca2(e)", "EX_cl(e)", "EX_fe2(e)", "EX_k(e)", "EX_na1(e)", "EX_h(e)", "EX_h2o(e)")
ex <- findExchReact(recon)
for(e in ex@react_id){
  lb <- ex@lowbnd[match(e, ex@react_id)]
  if(e %in% minimal) recon <- changeBounds(recon, react = e, lb = ifelse(lb<0,lb,-1), ub = 1000)
  else recon <- changeBounds(recon, react = e, lb = 0, ub = 1000)
}

recon <- readRDS("~/uni/dat/mod/r/recon2_2_corrected.RDS")
#recon <- tissue

uptake  <- c("EX_but(e)", "EX_glc(e)", "EX_h2o(e)", "EX_na1(e)", "EX_gln_L(e)")
product <- c("EX_co2(e)", "EX_lac_D(e)", "EX_lac_L(e)" ,"EX_cl(e)", "EX_bhb(e)", "EX_acac(e)", "EX_acetone(e)", "EX_hco3(e)", "EX_k(e)")
# "EX_bhb(e)", "EX_acac(e)", "EX_acetone(e)" => ketone bodies
for(s in c(uptake, product)){
  if(s %in% uptake) recon <- changeBounds(recon, react = s, lb = -1, ub = 1000)
  else recon <- changeBounds(recon, react = s, lb = 0, ub = 1000)
}

recon <- changeBounds(recon, react = "EX_glc(e)", lb = -0, ub = 1000)
recon <- changeBounds(recon, react = "EX_but(e)", lb = -0, ub = 1000)
recon <- changeBounds(recon, react = "EX_bhb(e)", lb = -0, ub = 1000) # hydroxybutyrate can also be substrate
recon <- changeBounds(recon, react = "EX_gln_L(e)",lb= -0, ub = 1000) # hydroxybutyrate can also be substrate

recon <- changeBounds(recon, react = "EX_hx(e)", lb = 0, ub = 1000) # hexanoate ~fe2/3 p+ shuffle bug
recon <- changeBounds(recon, react = "EX_fe2(e)", lb = -1, ub = 1000) # limit fe/3 p+ shuffle bug
#recon <- changeBounds(recon, react = "EX_o2(e)", lb = -1000, ub = 1000) # limit o2

recon <- changeBounds(recon, react = "EX_lnlc(e)", lb = 0, ub = 1000) # fatty acids not necessay
recon <- changeBounds(recon, react = "EX_lnlnca(e)", lb = 0, ub = 1000) # fatty acids not necessay
 
ex <- findExchReact(recon)
#ex[which(ex@react_id=="EX_but(e)"),]
sol <- optimizeProb(recon)
sol
exch <- getFluxDist(sol, ex)
flux <- getNetFlux(exch)
exch[flux@uptake]
exch[flux@product]

rel_ex <- c(uptake, product)
fva         <- fluxVar(recon, percentage=90, react = rel_ex) 
ex_max      <- maxSol(fva, "lp_obj")
ex_min      <- minSol(fva, "lp_obj")
pos_uptake <- ex_min[match(uptake,rel_ex)]; names(pos_uptake) <- uptake; pos_uptake
pos_product<- ex_max[match(product,rel_ex)];names(pos_product)<- product;pos_product


all <- getFluxDist(sol)
react_id(recon)[which(abs(all)>900)]

rel_ex <- names(exch[flux@uptake])
fva         <- fluxVar(recon, percentage=10, react = rel_ex) 
ex_max      <- maxSol(fva, "lp_obj")
ex_min      <- minSol(fva, "lp_obj")

exch[which(exch < -100)]
paste(names(exch[which(exch < -100)]), sep = ",")
exch[which(exch < -100 & ex_max < 0)]

pro <- c("EX_k(e)","EX_cl(e)", "EX_hco3(e)", "EX_lac_L(e)", "EX_co2(e)")
exch[match(pro, ex@react_id)]
sub <- c("EX_but(e)","EX_glc(e)", "EX_gln_L(e)", "EX_na1(e)", "EX_h2o(e)")
exch[match(sub, ex@react_id)]

ex[match(c(pro,sub), ex@react_id),]



#
# recon22 model update
#

# fix recon 2.2 (exchange reactions other way around defined)
ex <- findExchReact(recon)
ex2 <- ex[grep("^EX_", ex@react_id),]
for(i in 1:length(ex2@met_id)){
  old <- S(recon)[ex2@met_pos[i], ex2@react_pos[i]]
  S(recon)[ex2@met_pos[i], ex2@react_pos[i]] <- old * -1
}
lb <- recon@lowbnd[ex2@react_pos]; ub <- recon@uppbnd[ex2@react_pos]
recon@lowbnd[ex2@react_pos] <- -ub; recon@uppbnd[ex2@react_pos] <- -lb
# saveRDS(recon, "r/recon2_2_corrected.RDS")

# add u compartment
ex <- findExchReact(recon)
ex <- ex[grep("^EX_", ex@react_id),]
for(i in seq_along(ex@react_id)){
  rea <- gsub("\\(e\\)","\\(u\\)",ex@react_id[i])
  met <- ex@met_id[i]
  ub  <- ex@uppbnd[i]
  lb  <- ex@lowbnd[i]
  recon <- addReact(recon, id=rea, met = met, Scoef = -1, reversible = T, ub = ub, lb=lb)
}
# saveRDS(recon, "~/uni/dat/mod/r/recon2_2_updated.RDS")

recon <- readRDS("~/uni/dat/mod/r/recon2_2_updated.RDS")

ex <- findExchReact(recon)
ex <- ex[grep("^EX_.*\\(e\\)$", ex@react_id),]
lowbnd(recon)[ex@react_pos] <- 0
optimizeProb(recon)

c("EX_o2(e)", EX_fe2(u),,EX_h(u),EX_his_L(u),EX_ile_L(u),EX_leu_L(u),EX_lys_L(u),EX_met_L(u),EX_phe_L(u),
EX_pi(u),EX_thr_L(u),EX_trp_L(u),EX_val_L(u))



sol  <- optimizeProb(recon) # retOptSol=F
exch <- getFluxDist(sol, findExchReact(recon))
flux <- getNetFlux(exch)
exch[flux@uptake]
exch[flux@product]


printReaction(recon, react="EX_glc(u)")
