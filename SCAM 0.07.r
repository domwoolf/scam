#############################################
#       Initialisation
#############################################
data_path = "~/Documents/NIFA/data"
time_step = 'month'
clm.fit = FALSE
source('SCAM packages.r')
sourceCpp('SCAM1.01.cpp')
source('SCAM Functions 1.02.r')

exp.const = fread(paste0(data_path, "/experiment_constants.csv"))
exp.const[, max_tsmd := -(20+1.3*clay - 0.01*clay^2) * depth/23]
all.data = get_data(data_path)
all.data = initialise.monthly.dayprimc.data(all.data)
training.set = grep('hoos', unique(experiments), value = T)
training.set = unique(experiments)
training.set = sample(unique(experiments),length(experiments)%/%2)
training.set = training.set[!grepl('pendleton', training.set)]

#############################################
#
#      Solve for optimal parameter values
#
#############################################
scam_objective <- function(p, dt, experiments) {
  use_atsmd = TRUE
  monthly = TRUE
  for (experiment in experiments) {
    max_tsmd = exp.const[exp==experiment, max_tsmd]
    clay = exp.const[exp==experiment, clay]
    res = scam(dt[exp==experiment], 
               mic_fact = 1.0, 
               mic_offset = p['mic_offset'],
               kdissolution = p['kdissolution'],
               kdepoly = p['kdepoly'],
               kdeath_and_exudates = p['kdeath_and_exudates'],
               kdesorb = p['kdesorb'],
               ksorb = p['ksorb'],
               kmicrobial_uptake = p['kmicrobial_uptake'],
               cue_0 = p['cue_0'],
               mcue = p['mcue'], 
               mclay = p['mclay'], 
               clay, max_tsmd, use_atsmd, monthly)
    dt[exp==experiment, soc := res$soc]
  }
  sum((dt$soc - dt$measured.soc)^2, na.rm = TRUE)
}

train.data = all.data[exp %in% training.set]
p = fit.par[!(names(fit.par) %in% c(''))]
fit = optim(p, scam_objective, dt = train.data, experiments=training.set,
            method = "L-BFGS-B", control=list(maxit=500, factr=1e9, trace=1, REPORT=1),
            lower = fit.lower[names(p)],  upper = fit.upper[names(p)])
fit.par[names(fit$par)] <- fit$par
cat(paste0('fit.par = c(', paste0(names(fit.par), ' = ', round(fit.par,3), collapse = ', '), ')'))
#############################################
#
#      Run Model
#
#############################################
for (exp.index in 1:length(experiments)) {
  exp.data <- all.data[exp==experiments[exp.index]]
  primc <- as.data.frame(scam(exp.data, 1.0,
      fit.par[1], fit.par[2], fit.par[3],fit.par[4],fit.par[5], fit.par[6], fit.par[7], fit.par[8], fit.par[9], fit.par[10],
      exp.const$clay[exp.index], exp.const$max_tsmd[exp.index], 
      use_atsmd = TRUE, monthly = TRUE))
  all.data[exp==experiments[exp.index], (names(primc)) := primc]
}

all.data[, facet_label := gsub('expt_','',exp)]
all.data[, facet_label := gsub('_',' ',facet_label)]

ggplot(all.data, aes(x=month/12)) +
  geom_line(aes(y=soc), color="blue") +
  geom_point(aes(y=measured.soc), color="red", data=all.data[!is.na(measured.soc)]) +
  scale_y_continuous(name = expression(SOC~'(Mg ha'^{-1}*')'), limits=c(0,90)) +
  scale_x_continuous(name = 'years') +
  facet_wrap(~facet_label, ncol=6) +
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"),
    panel.border = element_rect(colour = 'grey50', fill=NA, size=0.5),
    strip.text = element_text(size=18),
    axis.text = element_text(size=14),
    axis.title =element_text(size=18))

modelled.soc <- all.data[!is.na(measured.soc) , soc]
measured.soc <- all.data[!is.na(measured.soc) , measured.soc]
# modelled.soc <- all.data[!is.na(measured.soc) & !(exp %in% training.set), soc]
# measured.soc <- all.data[!is.na(measured.soc) & !(exp %in% training.set), measured.soc]
mod.fit <- lm(modelled.soc~measured.soc)

prim.fit.data = data.table(
  modelled = modelled.soc,
  measured = measured.soc,
  model = 'SCAM')
all.fit.data = prim.fit.data

lm_eqn = function(df){
  m = lm(modelled ~ measured, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
    list(a = format(coef(m)[1], digits = 3),
      b = format(coef(m)[2], digits = 3),
      r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}
eq <- ddply(all.fit.data,.(model),lm_eqn)

ggplot(all.fit.data, aes(measured, modelled)) +
  geom_point() +
  geom_abline(slope = 1) +
  geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
  geom_text(data=eq, aes(30, 80, label=V1), parse = TRUE, inherit.aes=FALSE,size=5) +
  scale_x_continuous(name = expression(Measured~SOC~'(Mg ha'^{-1}*')'), limits = c(0,90)) +
  scale_y_continuous(name = expression(Modelled~SOC~'(Mg ha'^{-1}*')'), limits = c(0,90)) +
  facet_grid(. ~ model) +
  theme_bw() +
  theme(axis.text = element_text(size=14),
    axis.title = element_text(size=18),
    strip.text = element_text(size=18))
