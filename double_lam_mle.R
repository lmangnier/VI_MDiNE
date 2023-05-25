species_crohns = mdine::crohns$otu.counts
covar_crohns = cbind(mdine::crohns$covars[,2], mdine::crohns$covars[,3],mdine::crohns$covars[,4])
status = mdine::crohns$covars[,1]

ref.counts = species_crohns[,ncol(species_crohns)]
nonref.counts = species_crohns[,-ncol(species_crohns)]

r.counts = cbind(species_crohns[,-ncol(species_crohns)], species_crohns[,ncol(species_crohns)])

fit = nnet::multinom(r.counts~covar_crohns, trace=F)
summary(fit)
coefs = as.matrix(coef(fit), ncol=4)


residus = log((nonref.counts+1)/(ref.counts+1)) - tcrossprod(cbind(1,covar_crohns), coefs)

cov0 = cov(residus[status=="no",])
cov1 = cov(residus[status=="CD",])

prec0 = MASS::ginv(cov0)
prec1 = MASS::ginv(cov1)

lt0 = sum(abs(prec0[lower.tri(prec0)]))
lt1 = sum(abs(prec1[lower.tri(prec1)]))

lambda0 = (ncol(nonref.counts)^2 + ncol(nonref.counts)) / (lt0 + 0.5*sum(diag(prec0)))
lambda1 = (ncol(nonref.counts)^2 + ncol(nonref.counts)) / (lt1 + 0.5*sum(diag(prec1)))

lam_mle = (ncol(nonref.counts)^2 + ncol(nonref.counts))/((lt0+lt1)+sum(diag(prec0))+sum(diag(prec1))/2)
