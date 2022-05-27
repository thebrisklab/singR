### Functions for graph.

# function to create network matrices from vectorized lower diagonals:

#' Create network matrices from vectorized lower diagonals
#' \code{vec2net} transfer the matrix vectorized lower diagonals into net to show the component image.
#' @param invector vectorized lower diagonals.
#' @param make.diag default value = 1.
#'
#' @return a net matrx
#' @export
#'
#' @examples net = vec2net(1:10)
vec2net = function(invector,make.diag=1) {
  #invector: choose(p,2) x 1, where p is the number of nodes
  nNode = (1 + sqrt(1+8*length(invector)))/2
  outNet = matrix(0,nNode,nNode)
  outNet[lower.tri(outNet)] = invector
  dim(outNet) = c(nNode,nNode)
  outNet = outNet + t(outNet)
  diag(outNet) = make.diag
  outNet
}


############
#‘ Function for plotting networks with ggplot
#' create graph dataset with netmat and mmp_order
#' a data.frame called with vectorization of reordered netmat by mmp_order.
#' @param gmatrix netmat
#' @param sort_indices mmp_order
#'
#' @return a data.frame with vectors:
#'  ## X1: vector of numerics.
#'  ## X2: vector of numerics.
#'  ## value: vectorization of reordered netmat by mmp_order.
#'
#' @export
#'
create.graph.long = function(gmatrix,sort_indices=NULL) {
  nnode = nrow(gmatrix)
  X1 = c(1:nnode)%x%rep(1,nnode) # row number comes from Kronecker product
  X2 =  rep(1,nnode)%x%c(1:nnode) # column number
  if (!is.null(sort_indices)) {
    gmatrix = gmatrix[sort_indices,sort_indices]
  }
  value = as.vector(as.matrix(gmatrix))
  data.frame(X1,X2,value)
}




#IGAY: adjust function so that the path is not hard-coded
#' PlotNetwork for components
#'
#' @param component component for the plot, px x rx
#' @param title title for plot
#' @param qmin default value = 0.005
#' @param qmax default value = 0.995
#' @param path default path
#' @param make.diag default value = NA
#'
#' @return a list as followed:
#' \describe{
#'       \item{\code{netmatfig}}{component loadings from rs correlation}
#'       \item{\code{loadingsfig}}{the sum of the absolute values of the rows of netmatfig}
#'       \item{\code{netmat}}{a matrix can be used with image() function}
#'       \item{\code{loadingsummary}}{loadings summary for each row of netmatfig}
#' }
#' @export
#' @import ggplot2
#' @import grid
#' @import scales
plotNetwork = function(component,title='',qmin=0.005, qmax=0.995, path = '~/Dropbox/JINGCA/Data/community_affiliation_mmpplus.csv',make.diag=NA) {
  # component:
  # vectorized network of length choose(n,2)


  # load communities for plotting:
  mmp_modules = read.csv(path)
  mmp_order = order(mmp_modules$Community_Vector)

  #check community labels:
  #table(mmp_modules$Community_Label)
  #table(mmp_modules$Community_Label,mmp_modules$Community_Vector)

  labels = c('VI','SM','DS','VS','DM','CE','SC')
  coords = c(0,70.5,124.5,148.5,197.5,293.5,360.5)


  zmin = quantile(component,qmin)
  zmax = quantile(component,qmax)

  netmat = vec2net(component,make.diag)

  meltsub = create.graph.long(netmat,mmp_order)
  #g2 = ggplot(meltsub, aes(X1,X2,fill=value))+ geom_tile()+ scale_fill_gradient2(low = "blue",  high = "red",limits=c(zmin,zmax),oob=squish)+labs(title = paste0("Component ",component), x = "Node 1", y = "Node 2")+coord_cartesian(clip='off',xlim=c(-0,390))

  g2 = ggplot(meltsub, aes(X1,X2,fill=value))+
    geom_tile()+
    scale_fill_gradient2(low = "blue",  high = "red",limits=c(zmin,zmax),oob=squish)+
    labs(title = title, x = "Node 1", y = "Node 2")+
    coord_cartesian(clip='off',xlim=c(-0,390))

  for (i in 1:7) {
    if (i!=3) {
      g2 = g2+geom_hline(yintercept = coords[i],linetype="dotted",size=0.5)+geom_vline(xintercept = coords[i],linetype="dotted",size=0.5)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),ymin = (coords[i]+10), ymax = (coords[i]+10), xmin = 385, xmax = 385)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),xmin = (coords[i]+10), xmax = (coords[i]+10), ymin = -7, ymax = -7)
    } else{
      g2 = g2+geom_hline(yintercept = coords[i],linetype="dotted",size=0.5)+geom_vline(xintercept = coords[i],linetype="dotted",size=0.5)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),ymin = (coords[i]+10), ymax = (coords[i]+10), xmin = 385, xmax = 385)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),xmin = (coords[i]+1), xmax = (coords[i]+1), ymin = -7, ymax = -7)
    }
  }
  # which nodes are prominent:
  loadingsummary = apply(abs(netmat),1,sum,na.rm=TRUE)
  loadingsum2 = loadingsummary[mmp_order]

  Community = factor(mmp_modules$Community_Label)[mmp_order]

  g3 = qplot(c(1:379),loadingsum2,col=Community,size=I(3))+xlab('MMP Index')+ylab('L1 Norm of the Rows')

  return(list(netmatfig = g2, loadingsfig = g3, netmat=netmat, loadingsummary = loadingsummary))
}
#########
#' Plot for the changed example of Y component.
#'
#' @param component component for the plot, px x rx
#' @param title title for plot
#' @param qmin default value = 0.005
#' @param qmax default value = 0.995
#' @param path default path
#' @param make.diag default value = NA
#'
#' @return a list as followed:
#' \describe{
#'       \item{\code{netmatfig}}{component loadings from rs correlation}
#'       \item{\code{loadingsfig}}{the sum of the absolute values of the rows of netmatfig}
#'       \item{\code{netmat}}{a matrix can be used with image() function}
#'       \item{\code{loadingsummary}}{loadings summary for each row of netmatfig}
#' }
#' @export
#' @import ggplot2
#' @import grid
#' @import scales
#'
plotNetwork_change = function(component,title='',qmin=0.005, qmax=0.995, path = '~/Dropbox/JINGCA/Data/community_affiliation_mmpplus.csv',make.diag=NA) {
  # component:
  # vectorized network of length choose(n,2)


  # load communities for plotting:
  mmp_modules = read.csv(path)
  mmp_order = order(mmp_modules$Community_Vector)

  #check community labels:
  #table(mmp_modules$Community_Label)
  #table(mmp_modules$Community_Label,mmp_modules$Community_Vector)

  #labels = c('VI','SM','DS','VS','DM','CE','SC')
  #coords = c(0,70.5,124.5,148.5,197.5,293.5,360.5)


  zmin = quantile(component,qmin)
  zmax = quantile(component,qmax)

  netmat = vec2net(component,make.diag)

  meltsub = create.graph.long(netmat,mmp_order)
  #g2 = ggplot(meltsub, aes(X1,X2,fill=value))+ geom_tile()+ scale_fill_gradient2(low = "blue",  high = "red",limits=c(zmin,zmax),oob=squish)+labs(title = paste0("Component ",component), x = "Node 1", y = "Node 2")+coord_cartesian(clip='off',xlim=c(-0,390))

  g2 = ggplot(meltsub, aes(X1,X2,fill=value))+
    geom_tile()+
    scale_fill_gradient2(low = "blue",  high = "red",limits=c(zmin,zmax),oob=squish)+
    labs(title = title, x = "Node 1", y = "Node 2")+
    coord_cartesian(clip='off',xlim=c(-0,100))

  #for (i in 1:7) {
  #  if (i!=3) {
  #    g2 = g2+geom_hline(yintercept = coords[i],linetype="dotted",size=0.5)+geom_vline(xintercept = coords[i],linetype="dotted",size=0.5)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),ymin = (coords[i]+10), ymax = (coords[i]+10), xmin = 385, xmax = 385)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),xmin = (coords[i]+10), xmax = (coords[i]+10), ymin = -7, ymax = -7)
  #  } else{
  #    g2 = g2+geom_hline(yintercept = coords[i],linetype="dotted",size=0.5)+geom_vline(xintercept = coords[i],linetype="dotted",size=0.5)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),ymin = (coords[i]+10), ymax = (coords[i]+10), xmin = 385, xmax = 385)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),xmin = (coords[i]+1), xmax = (coords[i]+1), ymin = -7, ymax = -7)
  #  }
  #}
  # which nodes are prominent:
  loadingsummary = apply(abs(netmat),1,sum,na.rm=TRUE)
  loadingsum2 = loadingsummary[mmp_order]

  Community = factor(mmp_modules$Community_Label)[mmp_order]

  g3 = qplot(c(1:100),loadingsum2,col=Community,size=I(3))+xlab('MMP Index')+ylab('L1 Norm of the Rows')

  return(list(netmatfig = g2, loadingsfig = g3, netmat=netmat, loadingsummary = loadingsummary))
}



