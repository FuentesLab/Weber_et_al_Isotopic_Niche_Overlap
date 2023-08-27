


### Analysis associated with Weber et al. publication: “Isotopic niche overlap among foraging marine turtle species in the Gulf of Mexico”

### Codes modified from Lysy, Stasko, & Swanson:  https://cran.r-project.org/package=nicheROVER and vignettes from https://cran.r-project.org/web/packages/nicheROVER/vignettes/ecol-vignette.html

#load required packages
library(tidyverse)
library(nicheROVER)



####### initial NICHE PLOT CODE to re-implement niche.plot() function to include the cex argument, to increase the font size of the later plots ###
niche.plot <- function(niche.par, niche.data, alpha = .95,
                       species.names, iso.names, lims,
                       col, ndens = 512, pfrac = 0, xlab, cex = 1) {
  niso <- ncol(niche.par[[1]]$mu)
  nspec <- length(niche.par)
  npts <- 100 # number of points for each ellipse
  nell <- sapply(niche.par, function(x) nrow(x$mu)) # number of ellipses per species
  if(missing(species.names)) species.names <- names(niche.par)
  if(missing(iso.names)) iso.names <- colnames(niche.par[[1]]$mu)
  # create all the ellipses to get the plot limits right.
  ell <- vector("list", nspec)
  names(ell) <- names(species.names)
  D <- combn(niso, 2)
  for(ii in 1:nspec) {
    ell.tmp <- array(NA, c(nell[ii], ncol(D), npts+1, 2))
    for(jj in 1:nell[ii]) {
      for(kk in 1:ncol(D)) {
        ell.coord <- ellipse(niche.par[[ii]]$mu[jj, D[,kk]],
                             V = niche.par[[ii]]$Sigma[D[,kk], D[,kk], jj],
                             alpha = alpha, n = npts)
        ell.tmp[jj,kk,,] <- ell.coord
      }
    }
    ell[[ii]] <- ell.tmp
  }
  # plot limits.
  if(!missing(lims)) {
    # check that user didn't specify `col` positionally by mistake
    if(is.character(lims)) {
      stop("`lims` argument passed as a character vector.\nDid you mean to specify `col` instead?")
    }
  } else {
    # data
    lims <- array(sapply(niche.data, function(x) apply(x, 2, range)),
                  dim = c(2, niso, nspec))
    lims <- apply(lims, 2, range)
    # ellipse
    DD <- t(sapply(1:niso, function(ii) which(D == ii, arr.ind = TRUE)[1,]))
    elims <- array(sapply(ell, function(x) {
      tmp <- apply(x, c(4,2), range)
      rbind(as.matrix(tmp[1,,])[DD], as.matrix(tmp[2,,])[DD])
    }), dim = c(2, niso, nspec))
    elims <- apply(elims, 2, range)
    lims <- apply(rbind(lims, elims), 2, range)
  }
  # plots
  opar <- par(mfcol = c(niso,niso), mar = rep(.5, 4), oma = rep(4,4), cex = cex)
  for(ci in 1:niso) {
    for(ri in 1:niso) {
      # initialize plot
      plot.new()
      plot.window(lims[,ci], lims[,ri])
      if (ci == ri) {
        # diagonals: density plots
        xdens <- matrix(NA, ndens, nspec)
        ydens <- xdens
        for(ii in 1:nspec) {
          den <- density(niche.data[[ii]][,ci], n = ndens)
          xdens[,ii] <- den$x
          ydens[,ii] <- den$y
        }
        for(ii in 1:nspec) {
          ly <- par("usr")[1:2]
          ly[2] <- ly[1] + pfrac*(ly[2]-ly[1])
          ly[3] <- (ly[2]-ly[1])/nspec
          segments(x0 = niche.data[[ii]][,ci],
                   y0 = ly[1]+(ii-1)*ly[3], y1 = ly[1]+ii*ly[3], col = col[ii])
          ly <- ly[2] + ydens[,ii]/max(ydens)*(lims[2,ci]-ly[2])
          lines(xdens[,ii], ly, col = col[ii])
        }
      }
      if (ri > ci) {
        # lower triangle: point plots
        for(ii in 1:nspec) {
          points(niche.data[[ii]][,c(ci,ri)], col = col[ii], pch = 16)
        }
      }
      if (ri < ci) {
        # upper triangle: ellipses
        for(ii in 1:nspec) {
          for(jj in 1:nell[ii]) {
            lines(ell[[ii]][jj,which(D[1,] == ri & D[2,] == ci),,2:1], col = col[ii])
          }
        }
      }
      box()
      if(ci == niso) axis(side = 4, cex = cex) else axis(side = 4, labels = FALSE)
      if(ri == niso) axis(side = 1, cex = cex) else axis(side = 1, labels = FALSE)
      if(ri == 1) mtext(text = iso.names[ci], side = 3, line = 1, cex = cex)
      if(ci == 1) mtext(text = iso.names[ri], side = 2, line = 1, cex = cex)
    }
  }
  if(!missing(xlab)) {
    mtext(text = xlab, side = 1, outer = TRUE, line = 2.2, cex = 1.25 * cex)
  }
  legend(x = "topleft", legend = species.names, fill = col, bty = "n", cex = 0.9 * cex)
  par(opar) # reset par
}
### initial OVERLAP PLOT CODE to re-implement overlap.plot() function to include the cex argument, to increase the font size of the later plots ###

overlap.plot <- function(over.stat, nbreaks = 50, equal.axis = FALSE, species.names, col,
                         mean.cred = TRUE, mean.cred.col = "green", xlab, cex = 1) {
  if(length(dim(over.stat)) != 3 || dim(over.stat)[1] != dim(over.stat)[2])
    stop("Incorrect specification of over.stat.")
  nspec <- dim(over.stat)[1]
  if(missing(species.names)) species.names <- dimnames(over.stat)[[1]]
  # histograms
  over.hist <- apply(over.stat, 1:2,
                     function(x) {
                       if(any(is.na(x))) return(NULL)
                       else {
                         tmp <- hist(x*100, breaks = nbreaks, plot = FALSE)
                         tmp$density <- tmp$density/max(tmp$density)
                         tmp$counts <- tmp$counts/max(tmp$counts)
                       }
                       tmp
                     })
  # x-axis limits
  xlim <- lapply(over.hist,
                 function(h) {
                   if(is.null(h)) tmp <- matrix(NA, 2, 2)
                   else {
                     tmp <- cbind(range(h$breaks), range(h$density))
                   }
                   tmp
                 })
  xlim <- matrix(xlim, nspec, nspec)
  if(equal.axis) {
    for(ii in 1:nspec) {
      xlim[,ii] <- rep(list(cbind(range(sapply(xlim[,ii],
                                               function(ll) ll[,1]), na.rm = TRUE),
                                  range(sapply(xlim[,ii],
                                               function(ll) ll[,2]), na.rm = TRUE))),
                       nspec)
    }
  }
  # mean and credible intervals
  if(mean.cred) {
    over.mean <- apply(over.stat, 1:2, mean)*100
    over.cred <- apply(over.stat, 1:2, quantile, prob = c(.025, .975), na.rm = TRUE)*100
    over.cred <- array(over.cred, dim = c(2, nspec, nspec))
  }
  # plot
  par(mfcol = c(nspec,nspec), mar = c(1.5,rep(.5, 3)), oma = rep(4,4), cex = cex)
  for(ci in 1:nspec) {
    for(ri in 1:nspec) {
      # initialize plot
      plot.new()
      if (ri != ci) {
        # off diagonals: overlap histograms
        plot.window(xlim[ri,ci][[1]][,1], xlim[ri,ci][[1]][,2])
        if(equal.axis) {
          # recalculate histograms with new limits
          tmp <- hist(over.stat[ri,ci,]*100,
                      breaks = seq(xlim[ri,ci][[1]][1,1],
                                   xlim[ri,ci][[1]][2,1], len = nbreaks+1),
                      plot = FALSE)
          tmp$density <- tmp$density/max(tmp$density)
          over.hist[[ri,ci]] <- tmp
        }
        plot(over.hist[[ri,ci]], add = TRUE, freq = FALSE, col = col[ci],
             border = "white")
        if(mean.cred) {
          # mean and 95% credible intervals
          abline(v = c(over.mean[ri,ci], over.cred[,ri,ci]),
                 col = mean.cred.col, lty = c(1,2,2), lwd = 2)
        }
      } else {
        # diagonals: empty plots
        plot.window(xlim = c(0,1), ylim = c(0,1))
      }
      if(ri == 1 && ci == 1) {
        text(x = .5, y = .5,
             labels = expression("Niche Overlap: "*bgroup("(",
                                                          atop("Row", "Column"), ")")),
             adj = c(.5, .5), cex = 1)
      }
      box()
      if(ci != ri) axis(side = 1, cex = cex)
      if(ci > 1) axis(side = 2, cex = cex, labels = FALSE)
      if(ci < nspec) axis(side = 4, labels = FALSE)
      if(ri == 1) mtext(text = species.names[ci], side = 3, line = 1, cex = cex, col = col[ci])
      if(ci == 1) mtext(text = species.names[ri], side = 2, line = 1, cex = cex)
      if(mean.cred && ri == nspec && ci == nspec) {
        legend(x = "center", legend = c("Mean", "95% Credible Interval"),
               lty = c(1, 2), lwd = 2, col = mean.cred.col)
      }
    }
  }
  if(!missing(xlab)) {
    mtext(text = xlab, side = 1, line = 1.5, cex = 1.25 * cex, outer = TRUE)
  }
}

#load in data, make it into a dataframe
skin <- read.csv("SIA_Values.csv") 
turtle<- data.frame(skin)
print(turtle) #check to make sure its displaying correctly
summary(turtle)

#Calculate the means for each isotope and species
aggregate(turtle[2:4], turtle[1], mean)

skin %>% group_by(Species) %>%
  summarize(count = n(),
            mC = mean(d13C),
            sdC = sd(d13C), 
            rC = range(d13C),
            mN = mean(d15N), 
            sdN = sd(d15N),
            rN = range(d15N),
            mS = mean(d34S),
            sdS = sd(d34S), 
            rS = range(d34S))

set.seed(123)

#####Generate the posterior distributions of μ (mean) and Σ (variance) for isotope values for each species with the default prior
# generate parameter draws from the "default" posteriors of each species
nsamples <- 1e4
system.time({
  turtle.par <- tapply(1:nrow(turtle), turtle$Species,
                     function(ii) niw.post(nsamples = nsamples, X = turtle[ii,2:4]))
})

# various parameter plots
clrs <- c("red", "green", "blue") # colors for each Species

# mu1 (del15N), mu2 (del13C), and Sigma12
#par(mar = c(4, 4, .5, .1)+.1, mfrow = c(1,3))
niche.par.plot(turtle.par, col = clrs, plot.index = 1)
niche.par.plot(turtle.par, col = clrs, plot.index = 2)
niche.par.plot(turtle.par, col = clrs, plot.index = 1:2)
legend("topleft", legend = names(turtle.par), fill = clrs)

# all mu (del15N, del13C, del34S)
niche.par.plot(turtle.par, col = clrs, plot.mu = TRUE, plot.Sigma = FALSE) #posterior distribution of mu per species per isotope
legend("topleft", legend = names(turtle.par), fill = clrs)

# all mu and Sigma
par(mar = c(4.2, 4.2, 2, 1)+.1)
niche.par.plot(turtle.par, col = clrs, plot.mu = TRUE, plot.Sigma = TRUE) #posterior distribution of mu and sigma
legend("topright", legend = names(turtle.par), fill = clrs)

# 4. Create 2-d projections of a subset of niche regions
# See ?niche.plot for more details on parameter values.

# The resulting figure generates niche plots, density distributions, and raw data for each pairwise combination of isotope data for all three turtle Species (i.e., bivariate projections of 3-dimensional isotope data).

# 2-d projections of 10 niche regions

clrs <- c("#E31A1C", "#33A02C", "#0000FF") 

nsamples <- 10
turtle.par <- tapply(1:nrow(turtle), turtle$Species,
                   function(ii) niw.post(nsamples = nsamples, X = turtle[ii,2:4]))

# format data for plotting function
turtle.data <- tapply(1:nrow(turtle), turtle$Species, function(ii) X = turtle[ii,2:4])

niche.plot(niche.par = turtle.par, niche.data = turtle.data, pfrac = .1, alpha = 0.95, 
           iso.names = expression(delta^{13}*C, delta^{15}*N, delta^{34}*S),
           col = clrs, xlab = expression("Isotope Ratio (\u2030)"), cex = 0.98)



### Bottom-left = raw data points graphed
### Middle-diagonal = density distributions
### Top-right = niche plots


# 5. Calculate and display niche overlap estimates
# We use the function overlap() to calculate overlap metric estimates from a specified number of Monte Carlo draws (nsamples) from the turtle.par parameter list. It is necessary to specify the α-level. In this case, we have decided to calculate the overlap metric at two niche regions sizes for comparison: alpha=0.95 and alpha=0.99, or 95% and 99%.
# 
# Then, we calculate the mean overlap metric between each Species. Remember that the overlap metric is directional, such that it represents the probability that an individual from Species A will be found in the niche of Species B.

# niche overlap plots for 95% niche region sizes
nsamples <- 10000
turtle.par <- tapply(1:nrow(turtle), turtle$Species,
                   function(ii) niw.post(nsamples = nsamples, X = turtle[ii,2:4]))

# Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher accuracy.(doing 1000, 1e3 now bc increasing will take awhile/will need to let code run for awhile)
# the variable over.stat can be supplied directly to the overlap.plot function

over.stat <- overlap(turtle.par, nreps = nsamples, nprob = 1e4, alpha = c(.95, 0.99))

#The mean overlap metrics calculated across iteratations for both niche 
#region sizes (alpha = .95 and alpha = .99) can be calculated and displayed in an array.
over.mean <- apply(over.stat, c(1:2,4), mean)*100
round(over.mean, 2)

#the credible intervals overlap metrics
over.cred <- apply(over.stat*100, c(1:2, 4), quantile, prob = c(.025, .975), na.rm = TRUE)
round(over.cred[,,,1]) # display alpha = .95 niche region

# In the returned plot, Species A is along the rows and Species B is along columns. The plots represent the posterior probability that an individual from the Species indicated by the row will be found within the niche of the Species indicated by the column header. Before you plot, you must decide upon your α-level, and make sure the variable over.stat reflects this choice of α.

# Overlap plot. Before you run this, make sure that you have chosen your 
#alpha level.


clrs <- c("#E31A1C", "#33A02C", "blue") # colors for each Species
over.stat <- overlap(turtle.par, nreps = nsamples, nprob = 1e4, alpha = .95)

jpeg("RFig4prac.jpeg", width = 10, height = 8, units = 'in', res = 600)
overlap.plot(over.stat, col = clrs, mean.cred.col = "black", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%", cex = 1)
dev.off()
###axes have different ranges of values



# 6. Calculate and display niche size estimates
# See ?niche.size for exactly how niche size is defined as a function of the parameters μ and Σ. In a Bayesian context, we calculate the posterior distribution of niche size by Species. This done by calculating the niche size for every posterior sample of μ and Σ.

# posterior distribution of (mu, Sigma) for each Species
nsamples <- 10000
turtle.par <- tapply(1:nrow(turtle), turtle$Species,
                   function(ii) niw.post(nsamples = nsamples, X = turtle[ii,2:4]))

# posterior distribution of niche size by Species
turtle.size <- sapply(turtle.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})

# point estimate and standard error
rbind(est = colMeans(turtle.size), #point estimate of the mean of the size (or volume)
      se = apply(turtle.size, 2, sd)) #standard deviation of each column


# boxplots
dev.off()

clrs <- c("#E31A1C", "#33A02C", "#0000FF") # colors for each Species
boxplot(turtle.size, xlab = "Species", ylab = "Hypervolume Niche Size", col = clrs, pch = 16, cex = .5)



#how many times is the niche size of the loggerhead greater than green for each row (based on T/F)
sum(
  turtle.size[,"Cc"]>turtle.size[,"Cm"]
)/nsamples

sum(
  turtle.size[,"Cc"]>turtle.size[,"Lk"]
)/nsamples


sum(
  turtle.size[,"Cm"]>turtle.size[,"Cc"]
)/nsamples

sum(
  turtle.size[,"Cm"]>turtle.size[,"Lk"]
)/nsamples



sum(
  turtle.size[,"Lk"]>turtle.size[,"Cc"]
)/nsamples

sum(
  turtle.size[,"Lk"]>turtle.size[,"Cm"]
)/nsamples





