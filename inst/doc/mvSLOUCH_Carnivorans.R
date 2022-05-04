## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
## Code for installation of packages inside vignette taken from PCMBase vignette by Venelin Mitov
if(!requireNamespace("ggplot2")) {
  message("Building the vignette requires ggplot2 R-package. Trying to install.")
  status.ggplot2 <- try({
    install.packages("ggplot2")
  }, silent = TRUE)
  if(class(status.ggplot2 == "try-error")) {
    stop(
      "The ggplot2 installation did not succeed. The vignette cannot be built.")
  }
}

## -----------------------------------------------------------------------------
library(ggplot2)
library(ape)
library(mvSLOUCH)

## -----------------------------------------------------------------------------
load("./Carnivora_mvSLOUCH_objects.RData")

## -----------------------------------------------------------------------------
RNGversion(min(as.character(getRversion()),"3.6.1"))
set.seed(5, kind = "Mersenne-Twister", normal.kind = "Inversion")

## ----eval=TRUE, echo=TRUE, results="hide", message=FALSE----------------------
b_correct_dryad_download<-FALSE
temp <- tempfile()
tryCatch({
download.file("datadryad.org/api/v2/datasets/doi%253A10.5061%252Fdryad.77tm4/download",temp)
b_correct_dryad_download<-TRUE
},error=function(e){cat("Could not download data file from Dryad! No analysis can be done! Vignette not built!")},
warning=function(w){cat("Problem with downloading data file from Dryad! No analysis can be done! Vignette not built!")})

if (b_correct_dryad_download){
    b_correct_dryad_download<-FALSE
    tryCatch({
	dfcarnivores_postcranial <- read.table(unz(temp, "Carnivore Postcranial Morphology Data Samuels et al. 2013.txt"),header=TRUE,sep="\t",stringsAsFactors =FALSE)
	b_correct_dryad_download<-TRUE
    },error=function(e){cat("Corrupt data file from Dryad! No analysis can be done! Vignette not built!")},
    warning=function(w){cat("Problem with accessing data file from Dryad! No analysis can be done! Vignette not built!")})
}


## -----------------------------------------------------------------------------
unlink(temp)
Urocyon_cinereoargenteus_duplicated<-which(dfcarnivores_postcranial[,2]=="Urocyon cinereoargenteus")[2]
dfcarnivores_postcranial<-dfcarnivores_postcranial[-Urocyon_cinereoargenteus_duplicated,]
v_species_to_remove<-c("Bdeogale crassicauda","Lycalopex sp.","Daphoenus vetus","Barbourofelis loveorum","Archaeocyon leptodus","Canis armbrusteri","Canis dirus","Canis latrans orcutti","Canis lupus (Pleistocene)","Desmocyon thomsoni","Hesperocyon gregarius","Mesocyon coryphaeus","Paraenhydrocyon josephi","Phlaocyon leucosteus","Vulpes macrotis (Pleistocene)","Vulpes vulpes (Pleistocene)","Homotherium ischyrus","Homotherium serum","Leopardus wiedii (Pleistocene)","Lynx rufus (Pleistocene)","Miracinonyx inexpectatus","Panthera atrox","Puma concolor (Pleistocene)","Puma lacustris","Smilodon fatalis","Miacis gracilis","Mephitis mephitis (Pleistocene)","Spilogale gracilis (Pleistocene)","Spilogale putorius (Pleistocene)","Gulo gulo (Pleistocene)","Martes americana (Pleistocene)","Martes nobilis (Pleistocene)","Mustela nigripes (Pleistocene)","Neovison frenata (Pleistocene)","Neovison vison (Pleistocene)","Satherium piscinarium","Taxidea taxus (Pleistocene)","Dinictis felina","Dinictis major","Hoplophoneus primaevus","Nimravus brachyops","Agriotherium africanum","Arctodus simus","Ursus arctos (Pleistocene)")
v_indices_to_remove<-which(sapply(dfcarnivores_postcranial[,2],function(x,v_species_to_remove){is.element(x,v_species_to_remove)},v_species_to_remove=v_species_to_remove,simplify=TRUE))
dfcarnivores_postcranial<-dfcarnivores_postcranial[-v_indices_to_remove,]
m_names_change<-rbind(c("Alopex lagopus","Vulpes lagopus"),c("Lycalopex gymnocerus","Lycalopex gymnocercus"),c("Caracal serval","Leptailurus serval"),c("Felis silvestris libyca","Felis silvestris"),c("Atilax palundinosus","Atilax paludinosus"),c("Gallerella pulverulenta","Galerella pulverulenta"),c("Gallerella sanguinea","Galerella sanguinea"),c("Conepatus mesoleucus","Conepatus leuconotus"),c("Mephitis macroura vittata","Mephitis macroura"),c("Aonyx cinereus","Aonyx cinerea"),c("Paradoxurus hemaphroditus","Paradoxurus hermaphroditus"))
for (i in 1:nrow(m_names_change)){
    dfcarnivores_postcranial[which(dfcarnivores_postcranial[,2]==m_names_change[i,1]),2]<-m_names_change[i,2]
}
dfcarnivores_postcranial[,2]<-gsub(" ", "_", dfcarnivores_postcranial[,2])
row.names(dfcarnivores_postcranial)<-dfcarnivores_postcranial[,2]
dat<-dfcarnivores_postcranial[,c("Ecology","HuL","HuPCL","RaL")]

## -----------------------------------------------------------------------------
head(dat)

## -----------------------------------------------------------------------------
b_correct_tree_download<-FALSE
tryCatch({
phyltree_mammals<-ape::read.tree("http://www.biodiversitycenter.org/files/ttol/8.TTOL_mammals_unsmoothed.nwk")
b_correct_tree_download<-TRUE
},error=function(e){cat("Could not download tree file! No analysis can be done! Vignette not built!")},
warning=function(w){cat("Problem with downloading tree file! No analysis can be done! Vignette not built!")}
)

## -----------------------------------------------------------------------------
phyltree_mammals$tip.label[which(phyltree_mammals$tip.label=="Uncia_uncia")]<-"Panthera_uncia"
phyltree_mammals$tip.label[which(phyltree_mammals$tip.label=="Parahyaena_brunnea")]<-"Hyaena_brunnea"
phyltree_mammals$tip.label[which(phyltree_mammals$tip.label=="Bdeogale_crassicauda")]<-"Bdeogale_jacksoni"
tips_todrop<-setdiff(phyltree_mammals$tip.label,rownames(dat))
PrunedTree<-ape::drop.tip(phyltree_mammals,tips_todrop)
Tree<-ape::di2multi(PrunedTree)
dat<-dat[Tree$tip.label,]

## -----------------------------------------------------------------------------
mvSLOUCH::phyltree_paths(Tree)$tree_height

## -----------------------------------------------------------------------------
tree_height<-mvSLOUCH::phyltree_paths(Tree)$tree_height
ScaledTree<-Tree
ScaledTree$edge.length<-ScaledTree$edge.length/tree_height
mvSLOUCH::phyltree_paths(ScaledTree)$tree_height

## -----------------------------------------------------------------------------
isTRUE(all.equal(ScaledTree$tip.label,rownames(dat)))

## -----------------------------------------------------------------------------
row.names(dat)==ScaledTree$tip.label

## -----------------------------------------------------------------------------
regimes<-dat$Ecology[order(match(row.names(dat), ScaledTree$tip.label))]

## -----------------------------------------------------------------------------
regimesFitch<-mvSLOUCH::fitch.mvsl(ScaledTree,regimes)

## -----------------------------------------------------------------------------
regimesFitch<-mvSLOUCH::fitch.mvsl(ScaledTree,regimes,deltran=TRUE)

## -----------------------------------------------------------------------------
reg.col<-regimesFitch$branch_regimes
reg.col[reg.col=="generalist"]<-"purple"
reg.col[reg.col=="arboreal"]<-"green"
reg.col[reg.col=="cursorial"]<-"red"
reg.col[reg.col=="scansorial"]<-"orange"
reg.col[reg.col=="semiaquatic"]<-"blue"
reg.col[reg.col=="semifossorial"]<-"brown"

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  plot(ScaledTree, cex = 1,  edge.color = reg.col, edge.width=3.5, type="fan", font=4)

## ----eval=TRUE, echo=FALSE, out.width = "100%", fig.pos="h"-------------------
knitr::include_graphics("./ScaledTree_fan.png", auto_pdf=TRUE)

## -----------------------------------------------------------------------------
mvData<-data.matrix(dat[,c("HuPCL","RaL","HuL")])
mvData<-log(mvData)

## -----------------------------------------------------------------------------
mvStree<-mvSLOUCH::phyltree_paths(ScaledTree)

## -----------------------------------------------------------------------------
BMestim<-mvSLOUCH::BrownianMotionModel(mvStree,mvData)

## -----------------------------------------------------------------------------
BMestim$ParamsInModel$Sxx
BMestim$ParamsInModel$vX0

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OU1BM<-mvSLOUCH::mvslouchModel(mvStree, mvData, kY = 2, predictors = c(3))

## -----------------------------------------------------------------------------
OU1BM$FinalFound$ParamsInModel$A 

## -----------------------------------------------------------------------------
OU1BM$FinalFound$ParamSummary$phyl.halflife$halflives

## -----------------------------------------------------------------------------
OU1BM$FinalFound$ParamSummary$optimal.regression
OU1BM$FinalFound$ParamSummary$evolutionary.regression

## -----------------------------------------------------------------------------
OU1BM$FinalFound$ParamsInModel$mPsi

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OUBMestim <- mvSLOUCH::mvslouchModel(mvStree, mvData, kY = 2, predictors = c(3), regimes = regimesFitch$branch_regimes, root.regime = "generalist")

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OUBMestim$FinalFound
#  OUBMestim$MaxLikFound

## -----------------------------------------------------------------------------
OUBMestim$FinalFound$ParamsInModel$mPsi

## -----------------------------------------------------------------------------
OUBMestim$FinalFound$ParamSummary$phyl.halflife$halflives

## -----------------------------------------------------------------------------
OUBMestim$FinalFound$ParamSummary$optimal.regression
OUBMestim$FinalFound$ParamSummary$evolutionary.regression

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OU1OU <- mvSLOUCH::ouchModel(mvStree, mvData, predictors = c(3))

## -----------------------------------------------------------------------------
OU1OU$FinalFound$ParamsInModel$A

## -----------------------------------------------------------------------------
OU1OU$FinalFound$ParamSummary$phyl.halflife$halflives

## -----------------------------------------------------------------------------
OU1OU$FinalFound$ParamSummary$RSS$R2_non_phylogenetic_conditional_on_predictors

## -----------------------------------------------------------------------------
OU1BM$FinalFound$ParamSummary$RSS$R2_non_phylogenetic_conditional_on_predictors
OUBMestim$FinalFound$ParamSummary$RSS$R2_non_phylogenetic_conditional_on_predictors

## -----------------------------------------------------------------------------
OUBMestim$FinalFound$ParamsInModel$A

## -----------------------------------------------------------------------------
OUBMestim$FinalFound$ParamSummary$phyl.halflife$halflives

## -----------------------------------------------------------------------------
OU1OU$FinalFound$ParamSummary$evolutionary.regression

## -----------------------------------------------------------------------------
OU1OU$FinalFound$ParamsInModel$mPsi

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OUOUestim <- mvSLOUCH::ouchModel(mvStree, mvData, regimes = regimesFitch$branch_regimes, root.regime = "generalist", predictors = c(3))

## -----------------------------------------------------------------------------
OUOUestim$FinalFound$ParamsInModel$mPsi

## -----------------------------------------------------------------------------
OUOUestim$FinalFound$ParamSummary$phyl.halflife$halflives

## -----------------------------------------------------------------------------
OUOUestim$FinalFound$ParamSummary$RSS$R2_non_phylogenetic_conditional_on_predictors

## -----------------------------------------------------------------------------
OUOUestim$FinalFound$ParamSummary$evolutionary.regression

## -----------------------------------------------------------------------------
cbind(BMestim$ParamSummary$dof,
      OU1BM$FinalFound$ParamSummary$dof,
      OU1OU$FinalFound$ParamSummary$dof,
      OUBMestim$FinalFound$ParamSummary$dof,
      OUOUestim$FinalFound$ParamSummary$dof)

## -----------------------------------------------------------------------------
OUOUestim$FinalFound$ParamsInModel$Syy

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OUBMestim.mod <- mvSLOUCH::mvslouchModel(mvStree, mvData, kY = 2, predictors = c(3), regimes = regimesFitch$branch_regimes, root.regime = "generalist", Syytype = "LowerTri", Atype = "Diagonal")

## -----------------------------------------------------------------------------
OUBMestim.mod$MaxLikFound

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OU1 <- mvSLOUCH::estimate.evolutionary.model(mvStree, mvData, repeats = 5, model.setups = "basic", predictors = c(3), kY = 2, doPrint = TRUE)

## -----------------------------------------------------------------------------
capture.output(OU1,file = "OU1.txt")

## -----------------------------------------------------------------------------
OU1$BestModel$model

## -----------------------------------------------------------------------------
OU1$BestModel$BestModel$ParamSummary$phyl.halflife$halflives

## -----------------------------------------------------------------------------
OU1$BestModel$key.properties$R2_non_phylogenetic_conditional_on_predictors

## -----------------------------------------------------------------------------
OU1$BestModel$BestModel$ParamsInModel$A

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OU1$BestModel

## -----------------------------------------------------------------------------
OU1$BestModel$i

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OU1$testedModels[[11]]

## -----------------------------------------------------------------------------
OU1$testedModels[[21]]

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OU1$model.setups

## -----------------------------------------------------------------------------
OU1$BestModel$BestModel$ParamSummary$aic.c

## -----------------------------------------------------------------------------
OU1$BestModel$BestModel$ParamSummary$aic
OU1$BestModel$BestModel$ParamSummary$sic
OU1$BestModel$BestModel$ParamSummary$bic

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OUf <- mvSLOUCH::estimate.evolutionary.model(mvStree, mvData, regimes = regimesFitch$branch_regimes, root.regime = "generalist", repeats = 5, model.setups = "basic", predictors = c(3), kY = 2, doPrint = TRUE)

## -----------------------------------------------------------------------------
capture.output(OUf, file = "OUf.txt")

## -----------------------------------------------------------------------------
OUf$BestModel$model

## -----------------------------------------------------------------------------
OUf$BestModel$BestModel$ParamSummary$phyl.halflife$halflives
OUf$BestModel$key.properties$R2_non_phylogenetic_conditional_on_predictors

## -----------------------------------------------------------------------------
OUOUestim$FinalFound$ParamSummary$aic.c
OUf$BestModel$BestModel$ParamSummary$aic.c

## -----------------------------------------------------------------------------
OUOUestim$FinalFound$ParamSummary$aic.c-OUf$BestModel$BestModel$ParamSummary$aic.c
OU1$BestModel$BestModel$ParamSummary$aic.c-OUf$BestModel$BestModel$ParamSummary$aic.c

## -----------------------------------------------------------------------------
climb.reg <- regimesFitch$branch_regimes
climb.reg[climb.reg=="arboreal"] <- "climber"
climb.reg[climb.reg=="scansorial"] <- "climber"

## -----------------------------------------------------------------------------
climb.col <- climb.reg
climb.col[climb.col=="generalist"] <- "purple"
climb.col[climb.col=="climber"] <- "green"
climb.col[climb.col=="cursorial"] <- "red"
climb.col[climb.col=="semiaquatic"] <- "blue"
climb.col[climb.col=="semifossorial"] <- "brown"

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  plot(ScaledTree, cex = 1,  edge.color = climb.col, edge.width=3.5, type="fan", font=4)

## ----eval=TRUE, echo=FALSE, out.width = "100%", fig.pos="h"-------------------
knitr::include_graphics("./ScaledTree2_fan.png", auto_pdf=TRUE)

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OUc <- mvSLOUCH::estimate.evolutionary.model(mvStree, mvData, regimes = climb.reg, root.regime = "generalist", repeats = 5, model.setups = "basic", predictors = c(3), kY = 2, doPrint = TRUE)

## -----------------------------------------------------------------------------
OUc$BestModel$model

## -----------------------------------------------------------------------------
OUc$BestModel$BestModel$ParamSummary$aic.c

## -----------------------------------------------------------------------------
strok.reg<-regimesFitch$branch_regimes
strok.reg[strok.reg=="semiaquatic"]<-"stroker"
strok.reg[strok.reg=="semifossorial"]<-"stroker"

## -----------------------------------------------------------------------------
strok.col<-strok.reg
strok.col[strok.col=="generalist"]<-"purple"
strok.col[strok.col=="stroker"]<-"blue"
strok.col[strok.col=="cursorial"]<-"red"
strok.col[strok.col=="arboreal"]<-"green"
strok.col[strok.col=="scansorial"]<-"orange"

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  plot(ScaledTree, cex = 1,  edge.color = strok.col, edge.width=3.5, type="fan", font=4)

## ----eval=TRUE, echo=FALSE, out.width = "100%", fig.pos="h"-------------------
knitr::include_graphics("./ScaledTree3_fan.png", auto_pdf=TRUE)

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OUs <- mvSLOUCH::estimate.evolutionary.model(mvStree, mvData, regimes = strok.reg, root.regime = "generalist", repeats = 5, model.setups = "basic", predictors = c(3), kY = 2, doPrint = TRUE)

## -----------------------------------------------------------------------------
OUs$BestModel$model

## -----------------------------------------------------------------------------
OUs$BestModel$BestModel$ParamSummary$aic.c

## -----------------------------------------------------------------------------
OUf$BestModel$BestModel$ParamSummary$dof
OUs$BestModel$BestModel$ParamSummary$dof

## -----------------------------------------------------------------------------
red.reg<-regimesFitch$branch_regimes
red.reg[red.reg=="arboreal"]<-"climber"
red.reg[red.reg=="scansorial"]<-"climber"
red.reg[red.reg=="semiaquatic"]<-"stroker"
red.reg[red.reg=="semifossorial"]<-"stroker"

## -----------------------------------------------------------------------------
red.col<-red.reg
red.col[red.col=="generalist"]<-"purple"
red.col[red.col=="climber"]<-"green"
red.col[red.col=="cursorial"]<-"red"
red.col[red.col=="stroker"]<-"blue"

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  plot(ScaledTree, cex = 1,  edge.color = red.col, edge.width=3.5, type="fan", font=4)

## ----eval=TRUE, echo=FALSE, out.width = "100%", fig.pos="h"-------------------
knitr::include_graphics("./ScaledTree4_fan.png", auto_pdf=TRUE)

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OUr <- mvSLOUCH::estimate.evolutionary.model(mvStree, mvData, regimes = red.reg, root.regime = "generalist", repeats = 5, model.setups = "basic", predictors = c(3), kY = 2, doPrint = TRUE)

## -----------------------------------------------------------------------------
OUr$BestModel$model

## -----------------------------------------------------------------------------
OUr$BestModel$BestModel$ParamSummary$aic.c

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OUOUstart<-mvSLOUCH::ouchModel(mvStree, mvData, predictors = c(3), Atype = "Diagonal", diagA = NULL)

## -----------------------------------------------------------------------------
OUOUstart$FinalFound$ParamsInModel$A
OUOUstart$FinalFound$ParamsInModel$Syy

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OptOUs1<-mvSLOUCH::ouchModel(mvStree, mvData, regimes = strok.reg, root.regime = "generalist", predictors = c(3), Atype = OUs$BestModel$model$Atype, Syytype = OUs$BestModel$model$Syytype, diagA = OUs$BestModel$model$diagA, start_point_for_optim=list(A = OUOUstart$FinalFound$ParamsInModel$A, Syy = OUOUstart$FinalFound$ParamsInModel$Syy))

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OptOUf1<-mvSLOUCH::ouchModel(mvStree, mvData, regimes = regimesFitch$branch_regimes, root.regime = "generalist", predictors = c(3), Atype = OUf$BestModel$model$Atype, Syytype = OUf$BestModel$model$Syytype, diagA = OUf$BestModel$model$diagA, start_point_for_optim=list(A = OUOUstart$FinalFound$ParamsInModel$A, Syy = OUOUstart$FinalFound$ParamsInModel$Syy))

## -----------------------------------------------------------------------------
OptOUs1$FinalFound$LogLik
OptOUf1$FinalFound$LogLik

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OptOUs2<-mvSLOUCH::ouchModel(mvStree, mvData, regimes = strok.reg, root.regime = "generalist", predictors = c(3), Atype = OUs$BestModel$model$Atype, Syytype = OUs$BestModel$model$Syytype, diagA = OUs$BestModel$model$diagA, start_point_for_optim=list(A = OptOUs1$FinalFound$ParamsInModel$A, Syy = OptOUs1$FinalFound$ParamsInModel$Syy))
#  OptOUf2<-mvSLOUCH::ouchModel(mvStree, mvData, regimes = regimesFitch$branch_regimes, root.regime = "generalist", predictors = c(3), Atype = OUf$BestModel$model$Atype, Syytype = OUf$BestModel$model$Syytype, diagA = OUf$BestModel$model$diagA, start_point_for_optim=list(A = OptOUf1$FinalFound$ParamsInModel$A, Syy = OptOUf1$FinalFound$ParamsInModel$Syy))

## -----------------------------------------------------------------------------
OptOUs2$FinalFound$LogLik
OptOUf2$FinalFound$LogLik

## -----------------------------------------------------------------------------
OptOUs2$FinalFound$LogLik
OptOUf2$FinalFound$LogLik

## -----------------------------------------------------------------------------
OUs$BestModel$BestModel$LogLik
OUf$BestModel$BestModel$LogLik

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  OUOUreStart<-mvSLOUCH::ouchModel(mvStree, mvData, predictors = c(3), Atype = "Diagonal", diagA = NULL)

## -----------------------------------------------------------------------------
OUOUreStart$FinalFound$ParamsInModel$A
OUOUreStart$FinalFound$ParamsInModel$Syy

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  FinalOUs1<-mvSLOUCH::ouchModel(mvStree, mvData, regimes = strok.reg, root.regime = "generalist", predictors = c(3), Atype = OUs$BestModel$model$Atype, Syytype = OUs$BestModel$model$Syytype, diagA = OUs$BestModel$model$diagA, start_point_for_optim=list(A = OUOUreStart$FinalFound$ParamsInModel$A, Syy = OUOUreStart$FinalFound$ParamsInModel$Syy))
#  FinalOUf1<-mvSLOUCH::ouchModel(mvStree, mvData, regimes = regimesFitch$branch_regimes, root.regime = "generalist", predictors = c(3), Atype = OUf$BestModel$model$Atype, Syytype = OUf$BestModel$model$Syytype, diagA = OUf$BestModel$model$diagA, start_point_for_optim=list(A = OUOUreStart$FinalFound$ParamsInModel$A, Syy = OUOUreStart$FinalFound$ParamsInModel$Syy))

## -----------------------------------------------------------------------------
FinalOUs1$FinalFound$LogLik
FinalOUf1$FinalFound$LogLik

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  FinalOUs2<-mvSLOUCH::ouchModel(mvStree, mvData, regimes = strok.reg, root.regime = "generalist", predictors = c(3), Atype = OUs$BestModel$model$Atype, Syytype =OUs$BestModel$model$Syytype, diagA = OUs$BestModel$model$diagA, start_point_for_optim=list(A = FinalOUs1$FinalFound$ParamsInModel$A, Syy = FinalOUs1$FinalFound$ParamsInModel$Syy))
#  FinalOUf2<-mvSLOUCH::ouchModel(mvStree, mvData, regimes = regimesFitch$branch_regimes, root.regime = "generalist", predictors = c(3), Atype = OUf$BestModel$model$Atype, Syytype = OUf$BestModel$model$Syytype, diagA = OUf$BestModel$model$diagA, start_point_for_optim=list(A = FinalOUf1$FinalFound$ParamsInModel$A, Syy = FinalOUf1$FinalFound$ParamsInModel$Syy))

## -----------------------------------------------------------------------------
FinalOUs2$FinalFound$LogLik
FinalOUf2$FinalFound$LogLik

## -----------------------------------------------------------------------------
FinalOUf2$FinalFound$ParamSummary$aic.c - FinalOUs2$FinalFound$ParamSummary$aic.c

## -----------------------------------------------------------------------------
FinalOUs2$FinalFound$ParamSummary$RSS$R2

## -----------------------------------------------------------------------------
FinalOUs2$FinalFound$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval

## -----------------------------------------------------------------------------
DFs<-data.frame(
 Ecology=factor(colnames(FinalOUs2$FinalFound$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Estimated.Point)),
  RaL=FinalOUs2$FinalFound$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Estimated.Point["RaL",],
  upper=FinalOUs2$FinalFound$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Upper.end["RaL",],
  lower=FinalOUs2$FinalFound$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Lower.end["RaL",]
)

## -----------------------------------------------------------------------------
ggplot2::ggplot(DFs, ggplot2::aes(Ecology, RaL))+
  ggplot2::geom_point(size=4, colour=c("green","red","purple","orange","blue")) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=lower,ymax=upper),width=0.1,lwd=1.5, colour=c("green","red","purple","orange","blue"))+
  ggplot2::xlab("Locomotor habits")+
  ggplot2::ylab("RaL(log)") + ggplot2::coord_flip()

## -----------------------------------------------------------------------------
FinalOUs2$FinalFound$ParamSummary$evolutionary.regression
FinalOUs2$FinalFound$ParamSummary$corr.matrix

## -----------------------------------------------------------------------------
FinalOUs2$FinalFound$ParamSummary$trait.regression
FinalOUs2$FinalFound$ParamSummary$conditional.corr.matrix

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  BT<-mvSLOUCH::parametric.bootstrap(estimated.model = FinalOUs2, phyltree = mvStree,
#  values.to.bootstrap = c("evolutionary.regression", "corr.matrix", "trait.regression", "conditional.corr.matrix"),
#  regimes = strok.reg, root.regime = "generalist", predictors = c(3), numboot = 1000,
#  Atype = OUs$BestModel$model$Atype,
#  Syytype = OUs$BestModel$model$Syytype,
#  diagA = OUs$BestModel$model$diagA)

## -----------------------------------------------------------------------------
BTU.EvoReg <- FinalOUs2$FinalFound$ParamSummary$evolutionary.regression
BTU.EvoReg[] <- 0L
BTL.EvoReg <- BTU.EvoReg
for(i in 1:nrow(FinalOUs2$FinalFound$ParamSummary$evolutionary.regression)){
  BT.EvoReg<-quantile(sapply(BT$bootstrapped.parameters$evolutionary.regression,function(x) 
    x[i]), c(0.025, 0.975))
  BTL.EvoReg[i,] <- BT.EvoReg[1]
  BTU.EvoReg[i,] <- BT.EvoReg[2]
}

## -----------------------------------------------------------------------------
BTL.EvoReg
BTU.EvoReg

## -----------------------------------------------------------------------------
BTU.CorrMat <- rep(NA,length(as.vector(FinalOUs2$FinalFound$ParamSummary$corr.matrix)))
BTL.CorrMat<-BTU.CorrMat
for(i in 1:length(as.vector(FinalOUs2$FinalFound$ParamSummary$corr.matrix))){
  BT.CorrMat<-quantile(sapply(BT$bootstrapped.parameters$corr.matrix,function(x) x[i]),c(0.025,0.975))
  BTL.CorrMat[i] <- BT.CorrMat[1]
  BTU.CorrMat[i] <- BT.CorrMat[2]
}
BTL.CorrMat <- matrix(BTL.CorrMat, nrow =
  nrow(FinalOUs2$FinalFound$ParamSummary$corr.matrix))
BTU.CorrMat <- matrix(BTU.CorrMat, nrow =
  nrow(FinalOUs2$FinalFound$ParamSummary$corr.matrix))
dimnames(BTL.CorrMat) <- dimnames(BTU.CorrMat)<-
  list(row.names(FinalOUs2$FinalFound$ParamSummary$corr.matrix),
      colnames(FinalOUs2$FinalFound$ParamSummary$corr.matrix))

## -----------------------------------------------------------------------------
BTL.CorrMat
BTU.CorrMat

## -----------------------------------------------------------------------------
NA.TrtReg<-lapply(1:length(BT$bootstrapped.parameters$trait.regression), function(x) 
  rep(NA,length(unlist(FinalOUs2$FinalFound$ParamSummary$trait.regression))))
BTU.TrtReg <- rep(NA, length(unlist(FinalOUs2$FinalFound$ParamSummary$trait.regression)))
BTL.TrtReg <- BTU.TrtReg
for(i in 1:length(unlist(FinalOUs2$FinalFound$ParamSummary$trait.regression))){
  BT.TrtReg<-quantile(sapply(relist(unlist(BT$bootstrapped.parameters$trait.regression),
                                     NA.TrtReg), function(x) x[i]), c(0.025, 0.975))
  BTL.TrtReg[i] <- BT.TrtReg[1]
  BTU.TrtReg[i] <- BT.TrtReg[2]
}
BTL.TrtReg <- relist(BTL.TrtReg, FinalOUs2$FinalFound$ParamSummary$trait.regression)
BTU.TrtReg <- relist(BTU.TrtReg, FinalOUs2$FinalFound$ParamSummary$trait.regression)

## -----------------------------------------------------------------------------
BTL.TrtReg
BTU.TrtReg

## -----------------------------------------------------------------------------
BTU.CondCorr <-
  rep(NA,length(as.vector(FinalOUs2$FinalFound$ParamSummary$conditional.corr.matrix)))
BTL.CondCorr<-BTU.CondCorr 
for(i in 1:length(as.vector(FinalOUs2$FinalFound$ParamSummary$conditional.corr.matrix))){
  BT.CondCorr<-quantile(sapply(BT$bootstrapped.parameters$conditional.corr.matrix,function(x) x[i]),c(0.025,0.975))
  BTL.CondCorr[i] <- BT.CondCorr[1]
  BTU.CondCorr[i] <- BT.CondCorr[2]
}
BTL.CondCorr <- matrix(BTL.CondCorr, nrow = 
  nrow(FinalOUs2$FinalFound$ParamSummary$conditional.corr.matrix))
BTU.CondCorr<-matrix(BTU.CondCorr,nrow =
  nrow(FinalOUs2$FinalFound$ParamSummary$conditional.corr.matrix))
dimnames(BTL.CondCorr)<-dimnames(BTU.CondCorr)<-list(row.names(FinalOUs2$FinalFound$ParamSummary$conditional.corr.matrix), 
  colnames(FinalOUs2$FinalFound$ParamSummary$conditional.corr.matrix))

## -----------------------------------------------------------------------------
BTL.CondCorr
BTU.CondCorr

## ----conditional_print_treeerror, eval=!b_correct_tree_download, echo=FALSE----
#  cat("Error: Could not download tree file! No analysis can be done! Vignette not built!")

## ----conditional_print_Dryaderror, eval=!b_correct_dryad_download, echo=FALSE----
#  cat("Error: Could not download data file from Dryad! No analysis can be done! Vignette not built!")
#  unlink(temp)

