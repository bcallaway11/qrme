## does the heavy lifting for computing nonlinear models with measurement error
compute.nlme <- function(data, Yformla, Tformla, tau, tvals, copType="gumbel", ndraws=250,
                         reportTmat=TRUE, reportSP=TRUE, reportUM=TRUE, reportPov=TRUE,
                         povline=log(20000), reportQ=c(.1,.5,.9),
                         Ynmix=1, Tnmix=1, tol=1, drawsd=1, messages=FALSE) {
    
    yname <- lhs.vars(Yformla)
    tname <- lhs.vars(Tformla)
    Qyx <- qrme(Yformla, data=data, tau=tau, nmix=Ynmix,
                tol=tol, drawsd=drawsd, messages=messages, se=FALSE) ## don't bootstrap these
    Qtx  <- qrme(Tformla, data=data, tau=tau, nmix=Tnmix,
                 tol=tol, drawsd=drawsd, messages=messages, se=FALSE)

    if (reportTmat) {
        meTmat <- tmat(Qyx$Ystar, Qtx$Ystar)
        nomeTmat <- tmat(data[,yname], data[,tname])
    }

    if (reportSP) {
        mePs <- cor(Qyx$Ystar, Qtx$Ystar, method="spearman")
        nomePs  <- cor(data[,yname], data[,tname], method="spearman")
    }

    if (reportUM) {
        meUm <- upMob(Qyx$Ystar, Qtx$Ystar)
        nomeUm <- upMob(data[,yname], data[,tname])
    }

    
    ## drop simulated draws as these take up a lot of memory and are no longer needed
    Qyx$Ystar <- NULL
    Qtx$Ystar <- NULL

    ## now get joint distribution
    meres <- qr2me(yname, tname, Yformla, tau=tau, data=data,
                   tvals=tvals,
                   method="Nelder-Mead", copula=copType, Qyx=Qyx, Qtx=Qtx, ndraws=ndraws, retFytxlist=FALSE,
                   messages=FALSE)


    ## get results without measurement error
    rqyx <- rq(Yformla, tau=tau, data=data)
    rqtx <- rq(Tformla, tau=tau, data=data)
    class(rqyx) <- c("merr", "rqs")
    class(rqtx) <- c("merr", "rqs")
    rqyx$pi <- rqtx$pi <- 1
    rqyx$mu <- rqtx$mu <- 0
    rqyx$sig <- rqtx$sig <- 0

    nomeres <- qr2me(yname, tname, Yformla, tau=tau, data=data,
                     tvals=tvals, ndraws=ndraws,
                     method="Nelder-Mead", copula=copType, Qyx=rqyx, Qtx=rqtx, retFytxlist=FALSE,
                     messages=FALSE)

    browser()

    ## results just using quantile regression
    qrformla <- toformula(yname, c(tname, rhs.vars(Yformla)))
    qrytx <- rq(qrformla, tau=tau, data=data)
    qrytx$Fyt <- lapply(1:length(tvals), function(i) {
        newdta <- data
        newdta[,tname] <- tvals[i]
        Fytx <- predict(qrytx, newdata=newdta, type="Fhat", stepfun=TRUE)
        Fytx <- rearrange(Fytx)
        combineDfs(seq(min(data[,yname]), max(data[,yname]), length.out=500), Fytx)
    })
    


    if (reportPov) {
        mePovrate <- sapply(1:(length(meres$tvals)), function(i) meres$Fyt[[i]](povline))
        nomePovrate <- sapply(1:(length(nomeres$tvals)), function(i) nomeres$Fyt[[i]](povline))
        qrPovrate <- sapply(1:(length(nomeres$tvals)), function(i) qrytx$Fyt[[i]](povline))
    }

    meresQ <- sapply(meres$Fyt, function(Fy) quantile(Fy, type=1, probs=reportQ))
    nomeresQ <- sapply(nomeres$Fyt, function(Fy) quantile(Fy, type=1, probs=reportQ))
    qrytxQ <- sapply(qrytx$Fyt, function(Fy) quantile(Fy, type=1, probs=reportQ))

    meCopParam <- meres$cop.param
    nomeCopParam <- nomeres$cop.param
    
    out <- list(Yformla=Yformla, Tformla=Tformla, tau=tau, tvals=tvals, copType=copType, meCopParam=meCopParam, nomeCopParam=nomeCopParam, Ynmix=Ynmix, Tnmix=Tnmix, reportQ=reportQ,
                meQyx=Qyx, meQtx=Qtx, meresQ=meresQ, nomeresQ=nomeresQ, qrytxQ=qrytxQ, meTmat=meTmat, nomeTmat=nomeTmat,
                mePs=mePs, nomePs=nomePs, meUm=meUm, nomeUm=nomeUm, rqyx=rqyx, rqtx=rqtx, #nomeres=nomeres,
                qrytx=qrytx, mePovrate=mePovrate, nomePovrate=nomePovrate, qrPovrate=qrPovrate)

    out

}



## function to compute nonlinear models with two sided measurement error
nlme <- function(data, Yformla, Tformla, tau, tvals, copType="gumbel", ndraws=250,
                 reportTmat=TRUE, reportSP=TRUE, reportUM=TRUE, reportPov=TRUE, povline=log(20000),
                 Ynmix=1, Tnmix=1, tol=1, drawsd=1, messages=FALSE,
                 se=FALSE, biters=100, cl=1) {

    res <- compute.nlme(data=data,
                        Yformla=Yformla,
                        Tformla=Tformla,
                        tau=tau,
                        tvals=tvals,
                        copType=copType,
                        ndraws=ndraws,
                        reportTmat=reportTmat,
                        reportSP=reportSP,
                        reportUM=reportUM,
                        reportPov=reportPov,
                        povline=povline,
                        Ynmix=Ynmix,
                        Tnmix=Tnmix,
                        tol=tol,
                        drawsd=drawsd,
                        messages=messages)
                        

    if (se) {

        eachIter <- pbapply::pblapply(1:biters, function(b) {
            n <- nrow(data)
            brows <- sample(1:n, size=n, replace=TRUE)
            bdata <- data[brows,]
            tryCatch({
                out <- compute.nlme(data=bdata,
                                    Yformla=Yformla,
                                    Tformla=Tformla,
                                    tau=tau,
                                    tvals=tvals,
                                    copType=copType,
                                    ndraws=ndraws,
                                    reportTmat=reportTmat,
                                    reportSP=reportSP,
                                    reportUM=reportUM,
                                    reportPov=reportPov,
                                    povline=povline,
                                    Ynmix=Ynmix,
                                    Tnmix=Tnmix,
                                    tol=tol,
                                    drawsd=drawsd,
                                    messages=messages)
                out
            }, error=function(cond) {
                return(-99) ## use this as code for error on that bootstrap iteration
                ##return(cond)
            })
        }, cl=cl)

        #browser()
        ## drop list elements where bootstrap failed
        eachIter <- eachIter[!sapply(eachIter, is.null)]

        ## first step estimators
        ## outcome
        res$meQyx$sig.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$meQyx$sig)), 2, sd)
        res$meQyx$mu.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$meQyx$mu)), 2, sd)
        res$meQyx$pi.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$meQyx$pi)), 2, sd)
        res$meQyx$bet.se <- apply(simplify2array(lapply(eachIter, function(e) e$meQyx$bet)), 1:2, sd)
        ## treatment
        res$meQtx$sig.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$meQtx$sig)), 2, sd)
        res$meQtx$mu.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$meQtx$mu)), 2, sd)
        res$meQtx$pi.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$meQtx$pi)), 2, sd)
        res$meQtx$bet.se <- apply(simplify2array(lapply(eachIter, function(e) e$meQtx$bet)), 1:2, sd)
        ## qr
        res$qrytx$bet.se <- apply(simplify2array(lapply(eachIter, function(e) coef(e$qrytx))), 1:2, sd)
        ## nome
        res$rqyx$bet.se <- apply(simplify2array(lapply(eachIter, function(e) t(coef(e$rqyx)))), 1:2, sd)
        res$rqtx$bet.se <- apply(simplify2array(lapply(eachIter, function(e) t(coef(e$rqtx)))), 1:2, sd)
        
        ## transition matrix
        res$meTmat.se <- apply(simplify2array(lapply(eachIter, function(e) e$meTmat)), 1:2, sd)
        res$nomeTmat.se <- apply(simplify2array(lapply(eachIter, function(e) e$nomeTmat)), 1:2, sd)

        ## spearmans rho
        res$mePs.se <- sd(sapply(eachIter, function(e) e$mePs))
        res$nomePs.se <- sd(sapply(eachIter, function(e) e$nomePs))

        ## upward mobility
        res$meUm.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$meUm)), 2, sd)
        res$nomeUm.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$nomeUm)), 2, sd)

        ## poverty rate
        res$mePovrate.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$mePovrate)), 2, sd)
        res$nomePovrate.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$nomePovrate)), 2, sd)
        res$qrPovrate.se  <- apply(do.call(rbind, lapply(eachIter, function(e) e$qrPovrate)), 2, sd)

        ## quantiles
        res$meresQ.se <- apply(simplify2array(lapply(eachIter, function(e) e$meresQ)), 1:2, sd)
        res$nomeresQ.se <- apply(simplify2array(lapply(eachIter, function(e) e$nomeresQ)), 1:2, sd)
        res$qrytxQ.se <- apply(simplify2array(lapply(eachIter, function(e) e$qrytxQ)), 1:2, sd)

        res$biters <- length(eachIter)
        
    }

    res
}
