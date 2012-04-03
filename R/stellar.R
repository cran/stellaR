testComposition <- function(Z, Y, ML, AFE) {
# test for the existence of a given combination of
# chemical composition (Z, Y, AFE) and mixing length

    z <- c( (1:9)*1e-4, (1:9)*1e-3, 1e-2 )
    y <- cbind( c(rep(0.249,4), rep(0.250,6), seq(0.252,0.268, by=0.002) ), rep(0.250,19), rep(0.270,19), rep(0.33,19), rep(0.38,19), rep(0.42,19) )
    ml <- c(1.7, 1.8, 1.9)
    afe <- 0:1
    
    testAFE <- AFE %in% afe
    if(!testAFE) 
        warning("[alpha/Fe] value not in the database")
    
    testML <- ML %in% ml
    if(!testML) 
        warning("mixing-length value not in the database")
    
    testZ <- Z %in% z
    if(!testZ) 
        warning("z value not in the database")
    
    if( testZ ) {
        testY <- Y %in% y[z == Z,]
    } else {
        testY <- FALSE
    }
    if(!testY) 
        warning("y value not in the database")
    
    return( testAFE & testML & testZ & testY )
}

showComposition <- function() {
# Show the possible combinations of
# chemical composition (Z, Y, AFE) and mixing length

    z <- c( (1:9)*1e-4, (1:9)*1e-3, 1e-2 )
    y <- cbind( c(rep(0.249,4), rep(0.250,6), seq(0.252,0.268, by=0.002) ), rep(0.250,19), rep(0.270,19), rep(0.33,19), rep(0.38,19), rep(0.42,19) )
    colnames(y) <- rep("y", 6)
    ml <- c(1.7, 1.8, 1.9)
    afe <- 0:1
    
    cat("Mixing-length values:\n")
    cat("\t", paste(ml,collapse=", "), "\n\n")
    
    cat("alpha-enhancement values:\n")
    cat("\t", paste(afe,collapse=", "), " (i.e. [alpha/Fe] = 0.0 [alpha/Fe] = 0.3)", "\n\n")
    
    cat("Chemical compositions:\n")
    df <- as.data.frame(cbind(z,y))
    print(df, print.gap=2, row.names=rep("    ", dim(df)[1]))
}

getZahb <- function(z, y, ml, afe, baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {
# retrieve ZAHB for given parameters
  
    specificURL <- "zahb/ZAHB_Z"
    
    if(substr(baseURL, nchar(baseURL), nchar(baseURL)) != "/")
        baseURL <- paste(baseURL, "/", sep="")
    
    if( !testComposition(z, y, ml, afe))
        stop("required data not present in the database")
    
    Z <- format(z, nsmall=5, scientific=FALSE)	   
    Y <- format(y, nsmall=4)	   
    ML <- format(ml, nsmall=2)	   
    AFE <- ifelse(afe, "_AS09a3", "_AS09a0")
    
    URL <- paste(baseURL, specificURL, Z, "_He", Y, "_ML", ML, AFE, ".DAT", sep="")
  
    DATA <- read.table(URL, skip=6)
    colnames(DATA) <- c("mass", "logTe", "logL")
    L <- list(z=z, y=y, ml=ml, alpha.enh=ifelse(afe,0.3,0), data=DATA)
    class(L) <- c("zahb", "stellar")
    return(L)
}

getTrk <- function(m, z, y, ml, afe, baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {
# retrieve track (from PMS to He flash) for given parameters

    specificURL <- "trk/TRK_Z"
    
    if(substr(baseURL, nchar(baseURL), nchar(baseURL)) != "/")
        baseURL <- paste(baseURL, "/", sep="")
    
    if( !testComposition(z, y, ml, afe))
        stop("required data not present in the database")
    
    M <- format(m, nsmall=2)
    Z <- format(z, nsmall=5, scientific=FALSE)	   
    Y <- format(y, nsmall=4)	   
    ML <- format(ml, nsmall=2)	   
    AFE <- ifelse(afe, "_AS09a3", "_AS09a0")
    
    dirURL <- paste(baseURL, specificURL, Z, "_He", Y, "_ML", ML, AFE, "/", sep="")
    
    URL <- paste(dirURL, "OUT_M", M, "_Z", Z, "_He", Y, "_ML", ML, AFE, ".DAT", sep="")
    
    DATA <- read.table(URL, skip=5)
    colnames(DATA) <- c("mod", "time", "logL" ,"logTe", "mass", "Hc", "logTc", "logRHOc", "MHEc", "Lpp", "LCNO", "L3a", "Lg", "radius", "logg")
    L <- list(mass=m, z=z, y=y, ml=ml, alpha.enh=ifelse(afe,0.3,0), data=DATA)
    class(L) <- c("trk", "stellar")
    return(L)
}

getIso <- function(age, z, y, ml, afe, baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {
# retrieve iso (from PMS to He flash) for given parameters
  
    specificURL <- "iso/ISO_Z"
    
    if(substr(baseURL, nchar(baseURL), nchar(baseURL)) != "/")
        baseURL <- paste(baseURL, "/", sep="")
    
    if( !testComposition(z, y, ml, afe))
        stop("required data not present in the database")
    
    AGE <- age*1000
    zero <- ifelse(AGE < 10000, "0", "")
    Z <- format(z, nsmall=5, scientific=FALSE)	   
    Y <- format(y, nsmall=4)	   
    ML <- format(ml, nsmall=2)	   
    AFE <- ifelse(afe, "_AS09a3", "_AS09a0")
    
    dirURL <- paste(baseURL, specificURL, Z, "_He", Y, "_ML", ML, AFE, "/", sep="")
    
    URL <- paste(dirURL, "AGE", zero, AGE, "_Z", Z, "_He", Y, "_ML", ML, AFE, ".DAT", sep="")
    
    DATA <- read.table(URL, skip=6)
    colnames(DATA) <- c("logL" ,"logTe", "mass", "radius", "logg")
    L <- list(age=age, z=z, y=y, ml=ml, alpha.enh=ifelse(afe,0.3,0), data=DATA)
    class(L) <- c("iso", "stellar")
    return(L)
}


getHb <- function(m, z, y, ml, afe, baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {
# retrieve track (from ZAHB to thermal pulses) for given parameters
  
    specificURL <- "hb/TRK_Z"
    
    if(substr(baseURL, nchar(baseURL), nchar(baseURL)) != "/")
        baseURL <- paste(baseURL, "/", sep="")
    
    if( !testComposition(z, y, ml, afe))
        stop("required data not present in the database")
    
    M <- format(m, nsmall=2)
    Z <- format(z, nsmall=5, scientific=FALSE)	   
    Y <- format(y, nsmall=4)	   
    ML <- format(ml, nsmall=2)	   
    AFE <- ifelse(afe, "_AS09a3", "_AS09a0")
    
    dirURL <- paste(baseURL, specificURL, Z, "_He", Y, "_ML", ML, AFE, "_HB/", sep="")
  
  # search the mass of RGB progenitor...
    data(masshb)
    T <- c(m, z, y, ml)
    idx <- apply(masshb[,1:4], 1, function(x) all(as.numeric(x) == as.numeric(T)))
    masshb.ext <- masshb[idx,]
    sel <- masshb.ext[, 5] == substr(AFE, 2, nchar(AFE))
    massRGB <- format(masshb.ext[sel, 6], nsmall=4)
    
    URL <- paste(dirURL, "OUT_M", M, "_Z", Z, "_He", Y, "_ML", ML, AFE, "_ZAHB", massRGB, ".DAT", sep="")
    
    DATA <- read.table(URL, skip=5)
    colnames(DATA) <- c("mod", "time", "logL" ,"logTe", "mass", "Hc", "logTc", "logRHOc", "MHEc", "Lpp", "LCNO", "L3a", "Lg", "radius", "logg")
    L <- list(mass=m, massRGB=m, z=z, y=y, ml=ml, alpha.enh=ifelse(afe,0.3,0), data=DATA)
    class(L) <- c("hb", "stellar")
    return(L)
}


getHbgrid <- function(z, y, ml, afe, baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {
# retrieve track (from ZAHB to thermal pulses) for given parameters
  
    specificURL <- "hb/TRK_Z"
    
    if(substr(baseURL, nchar(baseURL), nchar(baseURL)) != "/")
        baseURL <- paste(baseURL, "/", sep="")
    
    if( !testComposition(z, y, ml, afe))
        stop("required data not present in the database")
    
    Z <- format(z, nsmall=5, scientific=FALSE)	   
    Y <- format(y, nsmall=4)	   
    ML <- format(ml, nsmall=2)	   
    AFE <- ifelse(afe, "_AS09a3", "_AS09a0")
    
    dirURL <- paste(baseURL, specificURL, Z, "_He", Y, "_ML", ML, AFE, "_HB/grid/", sep="")
  
  # search the mass of RGB progenitor...
    data(masshbgrid)
    T <- c(z, y, ml)
    idx <- apply(masshbgrid[,2:4], 1, function(x) all(as.numeric(x) == as.numeric(T)))
    masshb.ext <- masshbgrid[idx,]
    sel <- masshb.ext[, 5] == substr(AFE, 2, nchar(AFE))
    massRGB <- format(masshb.ext[sel, 6], nsmall=4)
    M <- format(masshb.ext[1, 1], nsmall=2)
    
    n.trk <- length(massRGB)
    L <- list()
    for(i in 1:n.trk) {
        URL <- paste(dirURL, "OUT_M", M, "_Z", Z, "_He", Y, "_ML", ML, AFE, "_ZAHB", massRGB[i], ".DAT", sep="")
        
        DATA <- read.table(URL, skip=5)
        colnames(DATA) <- c("mod", "time", "logL" ,"logTe", "mass", "Hc", "logTc", "logRHOc", "MHEc", "Lpp", "LCNO", "L3a", "Lg", "radius", "logg")
        L[[i]] <- list(mass=round(as.numeric(massRGB[i]),2), massRGB=M, z=z, y=y, ml=ml, alpha.enh=ifelse(afe,0.3,0), data=DATA)
        class(L[[i]]) <- c("hb", "stellar")
    }
    class(L) <- c("hbset", "stellar")
    return(L)
}


###############################


getTrkSet <- function(m, z, y, ml, afe, baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {

    grid <- expand.grid(m, z, y, ml, afe)
    n <- dim(grid)[1]
    
    if(n == 1) {
        if( !testComposition(grid[1,2], grid[1,3], grid[1,4], grid[1,5]))
            stop("required data not present in the database")
        return(getTrk(grid[1,1], grid[1,2], grid[1,3], grid[1,4], grid[1,5], baseURL))
    }
    
    trk <- list()
    for(i in 1:n) {
        if( !testComposition(grid[i,2], grid[i,3], grid[i,4], grid[i,5]))
            stop("required data not present in the database")
        trk[[i]] <- getTrk(grid[i,1], grid[i,2], grid[i,3], grid[i,4], grid[i,5], baseURL)
    }
    class(trk) <- c("trkset", "stellar")
    return(trk)
}

getIsoSet <- function(age, z, y, ml, afe, baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {

    grid <- expand.grid(age, z, y, ml, afe)
    n <- dim(grid)[1]
    
    if(n == 1) {
        if( !testComposition(grid[1,2], grid[1,3], grid[1,4], grid[1,5]))
            stop("required data not present in the database")
        return(getIso(grid[1,1], grid[1,2], grid[1,3], grid[1,4], grid[1,5], baseURL))
    }
    
    iso <- list()
    for(i in 1:n) {
        if( !testComposition(grid[i,2], grid[i,3], grid[i,4], grid[i,5]))
            stop("required data not present in the database")
        iso[[i]] <- getIso(grid[i,1], grid[i,2], grid[i,3], grid[i,4], grid[i,5], baseURL)
    }
    class(iso) <- c("isoset", "stellar")
    return(iso)
}
