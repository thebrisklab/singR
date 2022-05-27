

################################################
#' New dataset for simulation of FMRI
#'
#' \code{newSimFMRI} Create a dataset containing 4 components correspond to 2D images, roughly similar to a brain activation map.
#'
#' @param snr numeric, which shows the ratio of non-Gaussian (both joint and individual) to Gaussian components in dataset.
#' @param noisyICA whether to use noisyICA method to simulate gaussian noise.
#' @param nTR the number of time points, which is not changeable at here.
#' @param nImages the number of images for each component generated in the dataset.
#' @param phi default value is 0.5.
#' @param dim.data vector of length 2, the dimension of the simulation data.
#' @param var.inactive the background variance of non-active voxels or edges within a non-Gaussian component. Setting this equal to zero results in true sparsity in a non-Gaussian components.
#'
#' @export
#' @import neuRosim
#' @import steadyICA
newSimFMRI = function(snr = 1, noisyICA=FALSE, nTR=50, nImages=1, phi=0.5, dim.data=c(33,33), var.inactive=0.0001) {
  ##ASSUME 1,000 samples

  m = nImages
  #Latent components are fixed for each simulation:
  x1 = c(rep(3,5),4:7)
  y1 = c(3:7,rep(3,4))
  s1.coords = cbind(x1,y1)
  s1 = specifyregion(dim = dim.data, coord = s1.coords, form = "manual")
  s1[s1!=0] = seq(0.5,1,length=length(x1))

  x2 = c(8,8,8,9,10,9,10,10,10,9,8)
  y2 = c(15,14,13,13,13,15,15,16,17,17,17)
  s2.coords = cbind(c(x2,x2+7),c(y2,y2))
  s2 = specifyregion(dim=dim.data, coord = s2.coords, form = 'manual')
  s2[s2!=0] = seq(0.5,1,length=2*length(x2))

  x3 = c(13,14,15,15,15,14,13,15,15,14,13)
  y3 = c(19,19,19,20,21,21,21,22,23,23,23)
  s3.coords = cbind(c(x3,x3+7,x3+14),c(y3,y3,y3))
  s3 = specifyregion(dim=dim.data, coord = s3.coords, form = 'manual')
  s3[s3!=0] = seq(0.5,1,length=3*length(x3))

  ## add the fourth component
  x4=c(27:30,rep(30,5))
  y4=c(rep(3,4),4:8)
  s4.coords = cbind(x4,y4)
  s4 = specifyregion(dim=dim.data,coord = s4.coords, form = 'manual')
  s4[s4!=0] = seq(0.5,1,length=length(x4))





  sim.S = cbind(as.vector(s1),as.vector(s2),as.vector(s3),as.vector(s4))

  if(m>1) {
    t.sim.S = sim.S
    for(i in 1:(m-1)) t.sim.S = rbind(t.sim.S,sim.S)
    sim.S = t.sim.S
    rm(t.sim.S)
  }

  ## Add small amount of Gaussian noise to inactive voxels
  nInactive = sum(sim.S == 0)
  baseline = rnorm(nInactive,mean=0,sd=sqrt(var.inactive))
  sim.S[sim.S==0] = baseline

  ##For noise, simulate Gaussian random field. Unique for each simulation:
  if(noisyICA)  nscan = nTR else nscan = nTR-3
  sim.GRF = NULL
  for(k in 1:m) {
    t.sim.GRF <- spatialnoise(dim = dim.data, sigma=1, nscan = nscan, method = "gaussRF", FWHM = 6)
    dim(t.sim.GRF) <- c(prod(dim.data),nscan)
    sim.GRF = rbind(sim.GRF,t.sim.GRF)
  }

  ##Mixmat:
  #create timecourses for latent components:
  totaltime <- nTR
  nOnsets = 5+1
  onsets <- seq(from=1, to=totaltime, length=nOnsets)
  dur <- totaltime/10
  #s <- stimfunction(totaltime = totaltime, onsets = onsets, durations = dur, accuracy = 1)
  row1 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(1,3)]), durations = list(dur), effectsize = 1, TR = 1, conv = "gamma")
  row2 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(2,5)]), durations = list(dur), effectsize = 1, TR=1, conv='gamma')
  #NOTE: Time courses can not be identical.
  row3 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(2,4)]), durations=list(dur), effectsize=1, TR=1, conv='gamma')
  row4 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(3,6)]), durations=list(dur), effectsize = 1, TR=1,conv='gamma' )


  sim.Ms = matrix(c(row1,row2,row3,row4),nrow=4,byrow=TRUE)
  sim.Xs = sim.S%*%sim.Ms

  if(noisyICA)  {
    sim.Mn = NULL
    sim.Xn = sim.GRF
    for(t in 2:nTR) sim.Xn[,t] = phi*sim.Xn[,t-1]+sim.Xn[,t]
  }  else {
    sim.Mn = matrix(rnorm(nscan*nTR,0,1),nrow=nscan,ncol=nTR)
    for(t in 2:nTR) sim.Mn[,t] = phi*sim.Mn[,t-1] + sim.Mn[,t]
    sim.Xn = sim.GRF%*%sim.Mn
  }
  #sim.Xs = sim.Xs/sqrt(mean(sim.Xs^2))
  #sim.Xn = sim.Xn/sqrt(mean(sim.Xn^2))
  sim.Xs = sim.Xs/sd(as.vector(sim.Xs)) #standardize so we can control SNR
  sim.Xn = sim.Xn/sd(as.vector(sim.Xn))
  sim.Xs = sqrt(snr)*sim.Xs
  sim.X = sim.Xs + sim.Xn
  sim.X.whitened = whitener(X=sim.X)

  if(noisyICA) {
    return(list(S = sim.S, Ms = sim.Ms, X=sim.X, Mn = sim.Mn, N = sim.Xn, scaled.S = scale(sim.S),scaled.X = scale(sim.X), whitened.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
  } else {
    return(list(S = sim.S, Ms = sim.Ms, X=sim.X, Mn = sim.Mn, N = sim.GRF, scaled.S = scale(sim.S),scaled.X = scale(sim.X), whitened.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
  }
}
#--------------------------------------------

##############
###############
# BRISK: generateData_v2 alters individual component in second dataset to be a little more sparse,
# which makes it more realistic. (The original scenario was a more pathological example where logis fail
# but JB succeeds...)
# IGAY: added centering to mj so already column-centered approximately

#' Simulate data from the SING model.
#'
#' \code{generateData_v3} Create two datasets. The first dataset correspond to 2D images, roughly similar to a brain activation map, the second dataset corresponds to symmetric matrices, roughly similar to brain networks.
#'
#' @param nsubject number of subjects in the simulated dataset.
#' @param snr vector of length two corresponding to the ratio of non-Gaussian (both joint and individual) to Gaussian components in each dataset.
#' @param vars the background variance of non-active voxels or edges within a non-Gaussian component. Setting this equal to zero results in true sparsity in a non-Gaussian components.
#' @return a list with matrices
#' ## dX: the first dataset, nsubject x nPixels, here nPixels=33*33=1089
#' ## dY: the second dataset, nsubject x nEdges, here there are 100 nodes, and the vectorized lower triangular of the 100x100 matrix is 100*99/2.
#' ## mj: the true subject scores, nsubject x rJ, here, rJ=2 (two true components
#' ## sjX: true non-Gaussian joint components (loadings) for the first dataset
#' ## sjY: true non-Gaussian joint components (loadings) for the second dataset
#' ## siX: true non-Gaussian individual components (loadings) for the first dataset
#' ## siY: true non-Gaussian individual components (loadings) for the second dataset
#' ## snr: snr specified in the input
#' ## R2x: proportion of joint signal variance/(total variance) in dataset X
#' ## R2y: proportion of joint signal variance/(total variance) in dataset Y
#' @export
generateData_v3 <- function(nsubject = 48, snr = c(0.2, 0.2), vars = c(0.01,0.01)){
  # Generate mixing matrices
  n1 = round(nsubject/2)
  mj1 = c(rep( 1, n1), rep(-1, nsubject - n1)) + rnorm(nsubject) #joint subject scores
  mj2 = c(rep(-1, n1), rep( 1, nsubject - n1)) + rnorm(nsubject)
  mj = cbind(mj1, mj2)
  # mj = mj - matrix(colMeans(mj), nsubject, 2, byrow = T)

  # Create X components:
  # grab the 1, 2, 3, 4 components snr doesn't matter here as just grab the components in S
  simData = newSimFMRI(var.inactive = vars[1])
  # joint and individual signal components:
  px = nrow(simData$S)
  simS = scale(simData$S)

  # Create joint structure for X
  sjX = t(simS[,1:2]) #joint structure of 1 and 2 in X
  djX = mj%*%sjX

  # Create individual structure for X

  siX = t(simS[,3:4])

  n4 = round(nsubject/4)
  miX1 = c(rep(-1,n4),rep(1,n4),rep(-1,n4),rep(1,nsubject-3*n4))+rnorm(nsubject) # X independent subject score
  miX2 = c(rep(1,n4),rep(-1,n4),rep(1,n4),rep(-1,nsubject-3*n4))+rnorm(nsubject)
  miX = cbind(miX1,miX2)
  # miX = miX - mean(miX)
  diX = miX%*%siX

  # Calculate Frobenius norm of the signal
  signalXF2 = sum((djX + diX)^2)

  # Generate noise
  nX = t(scale(matrix(rnorm((nsubject-4)*px),px)))
  mnX = matrix(rnorm((nsubject-4)*nsubject),nsubject) ## mnX=matrix(n*(n-rx-1))
  # mnX = mnX - matrix(colMeans(mnX), nsubject, nsubject - 3, byrow = T)
  dnX = mnX%*%nX

  # Adjust the noise with snr ratio
  # Wrt to Frobenius norm
  dnX = dnX * sqrt(signalXF2/(sum(dnX^2)*snr[1]))


  # Create data matrix X
  dX = djX + diX + dnX

  # Calculate R^2 values for X joint
  R2x = sum(djX^2)/sum(dX^2)


  # Create Y components:
  # components that will represent network communities:
  # use a block structure:
  temp1 = c(rep(1,10),numeric(90))
  temp1 = temp1%*%t(temp1)

  temp2 = c(numeric(10),rep(1,20),numeric(70))
  temp2 = temp2%*%t(temp2)

  # BRISK: edited this component to be more sparse
  temp3 = c(numeric(40),rep(1,30),numeric(30))
  temp3 = temp3%*%t(temp3)

  temp4 = c(numeric(80),rep(1,20))
  temp4 = temp4%*%t(temp4)

  # Add small noise within the block structure
  var.noise = vars[2]

  # Create joint structure for Y
  sjY = cbind(temp3[lower.tri(temp3)],temp4[lower.tri(temp4)])
  inactive = sum(sjY==0)
  sjY[sjY==0] = rnorm(inactive, mean = 0, sd = sqrt(var.noise))
  sjY = t(scale(sjY))
  scalemj = t(t(mj)*c(-5,2))
  djY = scalemj%*%sjY
  py = ncol(sjY)

  # Create individual structure for Y
  siY = cbind(temp1[lower.tri(temp1)],temp2[lower.tri(temp2)])
  inactive = sum(siY==0)
  siY[siY==0] = rnorm(inactive, mean = 0, sd = sqrt(var.noise))
  siY = t(scale(siY))
  n8 = round(nsubject/8)
  miY = cbind(c(rep(1,n8),rep(-1,n8),rep(1,n8),rep(-1,n8),rep(1,n8),rep(-1,n8),rep(1,n8),rep(-1,nsubject-7*n8))+rnorm(nsubject),c(rep(1,n1),rep(-1,nsubject-n1))+rnorm(nsubject))
  # miY = miY - matrix(colMeans(miY), nsubject, 2, byrow = T)
  diY = miY%*%siY

  # Calculate Frobenius norm of the signal
  signalYF2 = sum((djY + diY)^2)

  # Generate noise for Y
  nY = t(scale(matrix(rnorm((nsubject-4)*py),py)))
  mnY = matrix(rnorm((nsubject-4)*nsubject),nsubject)
  # mnY = mnY - matrix(colMeans(mnY), nsubject, nsubject - 4, byrow = T)
  dnY = mnY%*%nY

  # Adjust the noise with snr ratio
  # Wrt to Frobenius norm
  dnY = dnY * sqrt(signalYF2/(sum(dnY^2)*snr[2]))

  # Create data matrix Y
  dY = djY + diY + dnY

  # Calculate R^2 values for X joint
  R2y = sum(djY^2)/sum(dY^2)

  return(list(dX = dX, dY = dY, mj = mj, sjX = sjX, sjY = sjY, siX = siX, siY = siY, snr = snr, R2x = R2x, R2y = R2y))
}


################################################
#' Sim data generation for 123 graph
#'
#' @param snr the signal to noise ratio
#' @param noisyICA whether use noisyICA
#' @param nTR default value is 50
#' @param nImages the number of images, default value is 1
#' @param phi default value = 0.5
#' @param dim.data dimension of data matrix
#' @param var.inactive default value is 0.0001
#'
#' @return a list of simulation data
#' @export
#' @import neuRosim
#' @import steadyICA
SimFMRI123 = function(snr = 1, noisyICA=FALSE, nTR=50, nImages=1, phi=0.5, dim.data=c(33,33), var.inactive=0.0001) {
  ##ASSUME 1,000 samples

  m = nImages
  #Latent components are fixed for each simulation:
  x1 = rep(3,5)
  y1 = c(3:7)
  s1.coords = cbind(x1,y1)
  s1 = specifyregion(dim = dim.data, coord = s1.coords, form = "manual")
  s1[s1!=0] = seq(0.5,1,length=length(x1))
  x2 = c(8,8,8,9,10,9,10,10,10,9,8)
  y2 = c(15,14,13,13,13,15,15,16,17,17,17)

  s2.coords = cbind(c(x2,x2+7),c(y2,y2))
  s2 = specifyregion(dim=dim.data, coord = s2.coords, form = 'manual')
  s2[s2!=0] = seq(0.5,1,length=2*length(x2))

  x3 = c(13,14,15,15,15,14,13,15,15,14,13)
  y3 = c(19,19,19,20,21,21,21,22,23,23,23)

  s3.coords = cbind(c(x3,x3+7,x3+14),c(y3,y3,y3))
  s3 = specifyregion(dim=dim.data, coord = s3.coords, form = 'manual')
  s3[s3!=0] = seq(0.5,1,length=3*length(x3))

  sim.S = cbind(as.vector(s1),as.vector(s2),as.vector(s3))

  if(m>1) {
    t.sim.S = sim.S
    for(i in 1:(m-1)) t.sim.S = rbind(t.sim.S,sim.S)
    sim.S = t.sim.S
    rm(t.sim.S)
  }

  ## Add small amount of Gaussian noise to inactive voxels
  nInactive = sum(sim.S == 0)
  baseline = rnorm(nInactive,mean=0,sd=sqrt(var.inactive))
  sim.S[sim.S==0] = baseline

  ##For noise, simulate Gaussian random field. Unique for each simulation:
  if(noisyICA)  nscan = nTR else nscan = nTR-3
  sim.GRF = NULL
  for(k in 1:m) {
    t.sim.GRF <- spatialnoise(dim = dim.data, sigma=1, nscan = nscan, method = "gaussRF", FWHM = 6)
    dim(t.sim.GRF) <- c(prod(dim.data),nscan)
    sim.GRF = rbind(sim.GRF,t.sim.GRF)
  }

  ##Mixmat:
  #create timecourses for latent components:
  totaltime <- nTR
  nOnsets = 5+1
  onsets <- seq(from=1, to=totaltime, length=nOnsets)
  dur <- totaltime/10
  #s <- stimfunction(totaltime = totaltime, onsets = onsets, durations = dur, accuracy = 1)
  row1 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(1,3)]), durations = list(dur), effectsize = 1, TR = 1, conv = "gamma")
  row2 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(2,5)]), durations = list(dur), effectsize = 1, TR=1, conv='gamma')
  #NOTE: Time courses can not be identical.
  row3 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(2,4)]), durations=list(dur), effectsize=1, TR=1, conv='gamma')

  sim.Ms = matrix(c(row1,row2,row3),nrow=3,byrow=TRUE)
  sim.Xs = sim.S%*%sim.Ms

  if(noisyICA)  {
    sim.Mn = NULL
    sim.Xn = sim.GRF
    for(t in 2:nTR) sim.Xn[,t] = phi*sim.Xn[,t-1]+sim.Xn[,t]
  }  else {
    sim.Mn = matrix(rnorm(nscan*nTR,0,1),nrow=nscan,ncol=nTR)
    for(t in 2:nTR) sim.Mn[,t] = phi*sim.Mn[,t-1] + sim.Mn[,t]
    sim.Xn = sim.GRF%*%sim.Mn
  }
  #sim.Xs = sim.Xs/sqrt(mean(sim.Xs^2))
  #sim.Xn = sim.Xn/sqrt(mean(sim.Xn^2))
  sim.Xs = sim.Xs/sd(as.vector(sim.Xs)) #standardize so we can control SNR
  sim.Xn = sim.Xn/sd(as.vector(sim.Xn))
  sim.Xs = sqrt(snr)*sim.Xs
  sim.X = sim.Xs + sim.Xn
  sim.X.whitened = whitener(X=sim.X)

  if(noisyICA) {
    return(list(S = sim.S, Ms = sim.Ms, X=sim.X, Mn = sim.Mn, N = sim.Xn, scaled.S = scale(sim.S),scaled.X = scale(sim.X), whitened.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
  } else {
    return(list(S = sim.S, Ms = sim.Ms, X=sim.X, Mn = sim.Mn, N = sim.GRF, scaled.S = scale(sim.S),scaled.X = scale(sim.X), whitened.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
  }
}


###############
# BRISK: generateData_v2 alters individual component in second dataset to be a little more sparse,
# which makes it more realistic. (The original scenario was a more pathological example where logis fail
# but JB succeeds...)
# IGAY: added centering to mj so already column-centered approximately
generateData_v2 <- function(nsubject = 48, snr = c(0.2, 0.2), vars = c(0.01,0.01)){
  # Generate mixing matrices
  n1 = round(nsubject/2)
  mj1 = c(rep( 1, n1), rep(-1, nsubject - n1)) + rnorm(nsubject)
  mj2 = c(rep(-1, n1), rep( 1, nsubject - n1)) + rnorm(nsubject)
  mj = cbind(mj1, mj2)
  # mj = mj - matrix(colMeans(mj), nsubject, 2, byrow = T)

  # Create X components:
  # grab the 1, 2, 3 components used in the LCA paper; snr doesn't matter here as just grab the components in S
  simData = SimFMRI123(var.inactive = vars[1]) #simulates LCA model with 3 LCs
  # joint and individual signal components:
  px = nrow(simData$S)
  simS = scale(simData$S)

  # Create joint structure for X
  sjX = t(simS[,2:3])
  djX = mj%*%sjX

  # Create individual structure for X
  siX = t(simS[,1])
  n4 = round(nsubject/4)
  miX = c(rep(-1,n4),rep(1,n4),rep(-1,n4),rep(1,nsubject-3*n4))+rnorm(nsubject)
  # miX = miX - mean(miX)
  diX = miX%*%siX

  # Calculate Frobenius norm of the signal
  signalXF2 = sum((djX + diX)^2)

  # Generate noise
  nX = t(scale(matrix(rnorm((nsubject-3)*px),px)))
  mnX = matrix(rnorm((nsubject-3)*nsubject),nsubject)
  # mnX = mnX - matrix(colMeans(mnX), nsubject, nsubject - 3, byrow = T)
  dnX = mnX%*%nX

  # Adjust the noise with snr ratio
  # Wrt to Frobenius norm
  dnX = dnX * sqrt(signalXF2/(sum(dnX^2)*snr[1]))


  # Create data matrix X
  dX = djX + diX + dnX

  # Calculate R^2 values for X joint
  R2x = sum(djX^2)/sum(dX^2)


  # Create Y components:
  # components that will represent network communities:
  # use a block structure:
  temp1 = c(rep(1,10),numeric(90))
  temp1 = temp1%*%t(temp1)

  temp2 = c(numeric(10),rep(1,20),numeric(70))
  temp2 = temp2%*%t(temp2)

  # BRISK: edited this component to be more sparse
  temp3 = c(numeric(40),rep(1,30),numeric(30))
  temp3 = temp3%*%t(temp3)

  temp4 = c(numeric(80),rep(1,20))
  temp4 = temp4%*%t(temp4)

  # Add small noise within the block structure
  var.noise = vars[2]

  # Create joint structure for Y
  sjY = cbind(temp1[lower.tri(temp1)],temp2[lower.tri(temp2)])
  inactive = sum(sjY==0)
  sjY[sjY==0] = rnorm(inactive, mean = 0, sd = sqrt(var.noise))
  sjY = t(scale(sjY))
  scalemj = t(t(mj)*c(-5,2))
  djY = scalemj%*%sjY
  py = ncol(sjY)

  # Create individual structure for Y
  siY = cbind(temp3[lower.tri(temp3)],temp4[lower.tri(temp4)])
  inactive = sum(siY==0)
  siY[siY==0] = rnorm(inactive, mean = 0, sd = sqrt(var.noise))
  siY = t(scale(siY))
  n8 = round(nsubject/8)
  miY = cbind(c(rep(1,n8),rep(-1,n8),rep(1,n8),rep(-1,n8),rep(1,n8),rep(-1,n8),rep(1,n8),rep(-1,nsubject-7*n8))+rnorm(nsubject),c(rep(1,n1),rep(-1,nsubject-n1))+rnorm(nsubject))
  # miY = miY - matrix(colMeans(miY), nsubject, 2, byrow = T)
  diY = miY%*%siY

  # Calculate Frobenius norm of the signal
  signalYF2 = sum((djY + diY)^2)

  # Generate noise for Y
  nY = t(scale(matrix(rnorm((nsubject-4)*py),py)))
  mnY = matrix(rnorm((nsubject-4)*nsubject),nsubject)
  # mnY = mnY - matrix(colMeans(mnY), nsubject, nsubject - 4, byrow = T)
  dnY = mnY%*%nY

  # Adjust the noise with snr ratio
  # Wrt to Frobenius norm
  dnY = dnY * sqrt(signalYF2/(sum(dnY^2)*snr[2]))

  # Create data matrix Y
  dY = djY + diY + dnY

  # Calculate R^2 values for X joint
  R2y = sum(djY^2)/sum(dY^2)

  return(list(dX = dX, dY = dY, mj = mj, sjX = sjX, sjY = sjY, siX = siX, siY = siY, snr = snr, R2x = R2x, R2y = R2y))
}

