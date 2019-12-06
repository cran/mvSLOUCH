## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


.OUtoBM_pcmbase_model_box<-function(pcmbase_model_box,StS,vX0,kX,vVars=NULL){
## length of vVars HAS to equal lengths of StS and vX0
## if we have vVars then the whole dimension changes!!
## vVars has to be a vector of numbers, NOT variable names

    pcmbase_model_box[which(names(pcmbase_model_box)=="H")]<-NULL 
    pcmbase_model_box[which(names(pcmbase_model_box)=="Theta")]<-NULL

    if (!is.null(vVars)){
	kX<-length(vX0)
	if ((length(vVars)!=kX) || (!all(dim(StS)==kX))){.my_stop(".OUtoBM_pcmbase_model_box: wrong number of BM variables passed.",TRUE)}
	
	## here we only set the dimension
	pcmbase_model_box$Sigma_x<-pcmbase_model_box$Sigma_x[vVars,vVars,,drop=FALSE]
	fullSigmae_x<-pcmbase_model_box$Sigmae_x
	pcmbase_model_box$Sigmae_x<-pcmbase_model_box$Sigmae_x[vVars,vVars,,drop=FALSE]
	for (cpcmbase_reg in names(pcmbase_model_box$Sigma_x[1,1,])){    	    
	    Sigmae_x<-fullSigmae_x[,,which(names(pcmbase_model_box$Sigmae_x[1,1,])==cpcmbase_reg)]
	    Sigmae<-Sigmae_x%*%t(Sigmae_x)
	    Sigmae<-Sigmae[vVars,vVars,drop=FALSE]
	    Sigmae_x<-.changeSigmatoSyy(Sigmae,"UpperTri",NULL,NULL,TRUE)
	    pcmbase_model_box$Sigmae_x[,,which(names(pcmbase_model_box$Sigmae_x[1,1,])==cpcmbase_reg)]<-Sigmae_x
	}
    }
    class(pcmbase_model_box)<-'BM'
    for (cpcmbase_reg in names(pcmbase_model_box$Sigma_x[1,1,])){
        ## Brownian motion with drift can also be considered here
            pcmbase_model_box$Sigma_x[,,which(names(pcmbase_model_box$Sigma_x[1,1,])==cpcmbase_reg)]<-.changeSigmatoSyy(StS,"UpperTri",NULL,NULL,FALSE)
            ## no need to touch the measurement error Sigmae object as this is the same in the models
    }
    pcmbase_model_box$X0<-vX0 ## at this stage vX0 is a vector
    pcmbase_model_box    
}

.set_1var_pcmbase_model_box<-function(pcmbase_model_box,varid,sdvalue=NULL){
## So far this is for OU (BM, OUOU, OUBM) type models
## if drift will be added this function needs to be expanded

    pcmbase_model_box_1var<-pcmbase_model_box
    if (is.element("Theta",names(pcmbase_model_box))){
        pcmbase_model_box_1var$Theta[,] <-pcmbase_model_box$Theta[varid,,drop=FALSE]
    }
    if (is.element("H",names(pcmbase_model_box))){
        pcmbase_model_box_1var$H <-pcmbase_model_box$H[varid,varid,,drop=FALSE]
    }
    if (is.element("Sigmae_x",names(pcmbase_model_box))){
        pcmbase_model_box_1var$Sigmae_x <-pcmbase_model_box$Sigmae_x[varid,varid,,drop=FALSE]
    }

    if (is.element("Sigma_x",names(pcmbase_model_box))){
	pcmbase_model_box_1var$Sigma_x <-pcmbase_model_box$Sigma_x[varid,varid,,drop=FALSE]
    	if (!is.null(sdvalue)){pcmbase_model_box_1var$Sigma_x[1,1,] <-sdvalue }
    }
    if (is.element("X0",names(pcmbase_model_box))){
        pcmbase_model_box_1var$X0<-pcmbase_model_box$X0[varid,drop=FALSE]
    }

    pcmbase_model_box_1var
}

.is0_Merror<-function(pcmbase_model_box){
    sum_merror<-0
    if (is.element("Sigmae_x",names(pcmbase_model_box))){
        for (cpcmbase_reg in names(pcmbase_model_box$Sigmae_x[1,1,])){
    	    sum_merror<-sum_merror+sum(abs(c(pcmbase_model_box$Sigmae_x[,,which(names(pcmbase_model_box$Sigmae_x[1,1,])==cpcmbase_reg),drop=FALSE])))
        }
    }
    isTRUE(all.equal(sum_merror,0))
}

.bm.estim<-function(mX,phyltree,pcmbase_model_box,regimes_types_orig,bRSScalc=TRUE,vVars=NULL,minLogLik=-Inf){
## for PCMbase mX can has to be a matrix 
## CANNOT be a data frame as when ouch::brown() was called
## PCMbase internally takes care of NA observations

    LogLik<- -Inf
    RSS <- Inf
    
    kX<-ncol(mX)
    n<-nrow(mX)
    
    bdoneX0<-FALSE
    if ((length(which(is.na(c(mX))))==0)&&(.is0_Merror(pcmbase_model_box))){
	pcmbase_model_box<-.OUtoBM_pcmbase_model_box(pcmbase_model_box,diag(1,kX,kX),rep(0,kX),kX,vVars)
	pcmbase_model_box_1var<-.set_1var_pcmbase_model_box(pcmbase_model_box,1,1)
	vX0<-matrix(NA,nrow=kX,ncol=1)
	## estimate root state separately for each dimension
	mRegressCovar<-matrix(0,kX,kX)
	for (i in 1:kX){	    
	    mD<-.design_matrix_construction(evolmodel="bm",n=phyltree$Ntips,kYX=1)[,-1,drop=FALSE]
	    lX0glsest<-.pcmbaseDphylOU_GLS(mX[,i,drop=FALSE],mD=mD,phyltree,pcmbase_model_box_1var,glsmodel=NA)
	    vX0[i,1]<-lX0glsest$vGLSest[1,1]
	    mRegressCovar<-lX0glsest$minDV1D		
	}

	mRes<-mX-matrix(vX0[,1],nrow=n,ncol=kX,byrow=TRUE)
	pcmbase_model_box_1var<-.set_1var_pcmbase_model_box(pcmbase_model_box,1,1)
	pcmbase_model_box_1var_mean0<-.set_mean0_pcmbase_model_box(pcmbase_model_box_1var,glsmodel=NA)
        StS<-matrix(0,kX,kX)
        for (i in 1:kX){
            vDi<-mRes[,i]
            mvD<-matrix(vDi,ncol=1)
            StS[i,i]<-.pcmbaseDphylGaussian_RSS(mvD,phyltree,pcmbase_model_box_1var_mean0,glsmodel=NA)            
        }
        if (kX>1){## double for loop cannot be done if kX is 1
            for (i in 1:(kX-1)){
                for(j in (i+1):kX){ ## upper triangle
                    vDi<-mRes[,i]
                    vDj<-mRes[,j]
                    mvDij<-matrix(vDi-vDj,ncol=1,byrow=TRUE)
                    rss_calc_ij<-.pcmbaseDphylGaussian_RSS(mvDij,phyltree,pcmbase_model_box_1var_mean0,glsmodel=NA)         
                    StS[i,j]<-(StS[i,i]+StS[j,j]-rss_calc_ij)/2
                    StS[j,i]<-StS[i,j]
                }
            }
        }
        StS<-StS/n
	pcmbase_model_box$X0<-vX0[,1]
	bdoneX0<-TRUE
	for (cpcmbase_reg in names(pcmbase_model_box$Sigma_x[1,1,])){
    		    ## Brownian motion with drift can also be considered here
    	    pcmbase_model_box$Sigma_x[,,which(names(pcmbase_model_box$Sigma_x[1,1,])==cpcmbase_reg)]<-.changeSigmatoSyy(StS,"UpperTri",NULL,NULL,FALSE)
	}
	LogLik<- minLogLik
	tryCatch({
	    LogLik<-.callPCMBase_mvlik(mX,phyltree, pcmbase_model_box,b_islog=TRUE,minLogLik=minLogLik)
	},error=function(e){.my_message(paste("Error in BM optim when calling PCMBse::mvlik",e),FALSE);.my_message("\n",FALSE)})
	optPar<-.sym.unpar(StS)
    }
    else{## there is measurement error or missing values
	## initial conditions are motivated by the 1D study case in
	## Bartoszek, K. and Sagitov S. (2015) "A consistent estimator of the evolutionary rate" Journal of Theoretical Biology 371:69-78.

	vX0<-apply(mX,2,mean,na.rm=TRUE)
	StS<-cov(mX,use="pairwise.complete.obs")
	## resulting matrix might not be symmetric-positive-definite due to presence of NA values
	## we do not delete whole observations but try to estimate pair-wise entries
	## resulting in potential loss of sym-pos-def, see ?cov
	vNAStS<-which(is.na(StS))
    
	if (length(vNAStS)>0){
	    StS[vNAStS]<-runif(length(vNAStS))/10
	}
	if ((!matrixcalc::is.symmetric.matrix(StS)) || (!matrixcalc::is.positive.definite(StS))){
	    StS<-as.matrix(Matrix::nearPD(StS)$mat)
	}
	mRegressCovar<-StS/n
    
	if(is.element("tree_height",names(phyltree))){StS<-StS/phyltree$tree_height}else{StS<-StS/(log(n))}
	Sxx<-t(.my_chol(StS))
	##Sxx<-t(chol(StS))

	## Brownian motion with drift can also be considered here
	## model was originally setup as OU now 'dropped' to BM
	pcmbase_model_box<-.OUtoBM_pcmbase_model_box(pcmbase_model_box,StS,vX0,kX,vVars)
    
	optPar<-NA
	tryCatch({
	    c_optim_method<-"Nelder-Mead"
	    if (length(.sym.unpar(StS))==1){c_optim_method<-"BFGS"}
	    optSxx<-optim(par=.sym.unpar(StS),fn=function(parStS,mX,phyltree,pcmbase_model_box){
		    LogLik<- 1000000
		    for (cpcmbase_reg in names(pcmbase_model_box$Sigma_x[1,1,])){
    			## Brownian motion with drift can also be considered here
        		pcmbase_model_box$Sigma_x[,,which(names(pcmbase_model_box$Sigma_x[1,1,])==cpcmbase_reg)]<-.changeSigmatoSyy(StS,"UpperTri",NULL,NULL,FALSE)
		    }
	    	
	    	    ## X0 has to be estimated anyway here, does not matter what its initial condition was		
		    mDesignBM<-.design_matrix_construction(evolmodel="bm",n=phyltree$Ntips,kYX=kX)[,-1,drop=FALSE]
		
	    	    ## in BM case no need to calculate an intercept of the GLS, at least for the moment
		    #mIntercept<-.calculate_intercept(evolmodel="bm",NA,NA)
		    
		    lX0glsest<-.pcmbaseDphylOU_GLS(mX,mDesignBM,phyltree,pcmbase_model_box)
		    ## Brownian motion with drift can also be considered here
    		
    		    pcmbase_model_box$X0<-lX0glsest$vGLSest[,1] ## the output of .pcmbaseDphylOU_GLS()$vGLSest is a matrix
    		    if (any(sapply(pcmbase_model_box$X0,is.na,simplify=TRUE))){pcmbase_model_box$X0<-apply(mX,2,mean,na.rm=TRUE)}
	    
		    tryCatch({
			LogLik<-(-1)*.callPCMBase_mvlik(mX,phyltree, pcmbase_model_box,b_islog=TRUE,minLogLik=minLogLik)
	    	    },error=function(e){.my_message(paste("Error in BM optim when calling PCMBse::mvlik",e),FALSE);.my_message("\n",FALSE)})
		    if (is.nan(LogLik)||is.na(LogLik)||(LogLik> 1000000)||is.infinite(LogLik)){LogLik<- 1000000}
		    LogLik
		},method=c_optim_method,mX=mX,phyltree=phyltree,pcmbase_model_box=pcmbase_model_box);
		optPar<-optSxx$par
		LogLik<- (-1)*optSxx$value
	    },error=function(e){.my_warning(paste("Error in BM optim ",e,"\n"),TRUE,TRUE);}	
	)
    }
    if (!is.na(optPar[1])){
	    StS<-.sym.par(optPar,kX)
	    if ((!matrixcalc::is.symmetric.matrix(StS)) || (!matrixcalc::is.positive.definite(StS))){
		StS<-as.matrix(Matrix::nearPD(StS)$mat)
	    }

	    Sxx<-t(.my_chol(StS))
	    ##Sxx<-t(chol(StS))
	    StS[which(abs(StS)<1e-15)]<-0
    	    Sxx[which(abs(Sxx)<1e-15)]<-0
	    
	    for (cpcmbase_reg in names(pcmbase_model_box$Sigma_x[1,1,])){
    		## Brownian motion with drift can also be considered here
        	pcmbase_model_box$Sigma_x[,,which(names(pcmbase_model_box$Sigma_x[1,1,])==cpcmbase_reg)]<-.changeSigmatoSyy(StS,"UpperTri",NULL,NULL,FALSE)
        	## no need to touch the measurement error Sigmae object as is the same as in other PCMBase models
	    }	    
	    if (!bdoneX0){
		pcmbase_model_box$X0<-vX0 ## here vX0 is constructed as a vector
		mDesignBM<-.design_matrix_construction(evolmodel="bm",n=phyltree$Ntips,kYX=kX)[,-1,drop=FALSE]
		lX0glsest<-.pcmbaseDphylOU_GLS(mX,mDesignBM,phyltree,pcmbase_model_box)
		vX0<-lX0glsest$vGLSest ## here vX0 is returned as a matrix
	        if (any(sapply(vX0[,1],is.na,simplify=TRUE))){vX0[,1]<-apply(mX,2,mean,na.rm=TRUE)}
	        mRegressCovar<-lX0glsest$minDV1D		
		
		pcmbase_model_box$X0<-vX0[,1]
	    }
	    RSS<-NA
	    if (bRSScalc){
		RSS<-list(RSS=NA,R2=NA)
		vIntercp<-apply(mX,2,mean,na.rm=TRUE)
                mInterceptCentredData<-mX-matrix(c(vIntercp),nrow=n,ncol=length(vIntercp),byrow=TRUE)
                RSS_null_model<-sum((mInterceptCentredData)^2,na.rm=TRUE)
		RSS$RSS<- .pcmbaseDphylGaussian_RSS(mX,phyltree,pcmbase_model_box)		
		RSS$R2<-1-(RSS$RSS)/RSS_null_model
		if (!is.element("RSS",names(RSS))||is.na(RSS$RSS)||(RSS$RSS<0)){
            	    RSS$RSS_comment<-"RSS is negative, consider rerunning estimation or a different (also non-phylogenetic) model of evolution"
                }
                if (!is.element("R2",names(RSS))||is.na(RSS$R2)||(RSS$R2<0)){
            	    RSS$R2_comment<-"R2 is negative, consider rerunning estimation or a different (also non-phylogenetic) model of evolution"
                }
	    }
    }    
    list(vX0=matrix(vX0,ncol=1,nrow=length(vX0)),StS=StS,Sxx=Sxx,LogLik=LogLik,RSS=RSS,regressCovar=mRegressCovar)
}

