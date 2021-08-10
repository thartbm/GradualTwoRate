
# Convert Methods Figures -----

library('rsvg')
library(magick)

rePlotMethodsFigX <- function(target='inline', fig=1) {
  
  if (target=='svg') {
    # nothing to do here:
    return()
  }
  
  dpi <- c(300,1200)[fig]
  
  # read in the figure from the SVG source file:
  FigX <- magick::image_read_svg(sprintf('doc/Fig%d.svg',fig))
  
  # extract width and height properties:
  width <- magick::image_info(FigX)[['width']]
  height <- magick::image_info(FigX)[['height']]
  
  # read style settings for the project to get the standard width:
  styles <- getStyle()
  # get the output figure width
  fw <- styles$figwidth[1]
  # determine what the output figure height should be:
  fh <- fw * (height/width)
  
  
  # create the graphics device (if not inline):
  if (target == 'tiff') {
    tiff( filename=sprintf('doc/Fig%d.tiff',fig),
          width=fw*dpi, height=fh*dpi, units='px',
          type='cairo',compression='lzw',res=dpi)
  }
  if (target == 'pdf') {
    pdf( file = sprintf('doc/Fig%d.pdf',fig),
         width=fw*dpi, height=fh*dpi)
  }
  
  # output the figure:
  par(mai=c(0,0,0,0))
  plot(FigX)
  
  # close non-inline graphics device:
  if (target %in% c('tiff','pdf')) {
    dev.off()
  }
  
}

# performance plot: -----

plotExpBehavior <- function(target='inline', exp, version=1) {
  
  info <- getInfo()
  style <- getStyle()
  
  groups <- info$group[which(info$experiment == exp)]
  
  fw_i <- 7
  if (version == 1) {
    fh_i <- 2.5*length(groups)
  }
  if (version == 2) {
    fh_i <- 2.5*(length(groups)+1)
  }
  
  if (target=='pdf') {
    pdf(file=sprintf('doc/Fig%d.pdf',1+(exp*2)), width=fw_i, height=fh_i)
  }
  
  data <- getSelectedGroupsData(groups=groups)
  
  data <- parseData(data, FUN=baseline)
  data <- parseData(data, FUN=addtrial)
  blocked <- parseData(data, FUN=block)
  
  ngroups = length(data)
  
  if (version == 1) {
    layout( mat=matrix(data=c(1:(ngroups*2)), nrow = ngroups, ncol=2, byrow = TRUE), widths=c(1,0.5) )
  }
  if (version == 2) {
    m <- c(unlist(lapply( ((c(1:ngroups)-1) * 3), function(x){rep(c(1,2)+as.numeric(x),each=3)} )), rep(c(1:ngroups)*3, each=6/ngroups) )
    layout( mat=matrix(data=m, nrow=ngroups+1, byrow=TRUE ) )
  }
  
  par(mar=c(4,3,1,0.2))
  
  for (group in names(data)) {
    
    group_idx <- which(info$group == group)
    rotation <- info$rotation[group_idx]
    group_label <- info$label[group_idx]
    
    if (rotation == 45) {
      ylimits <- c(0,45)
    } else {
      ylimits <- c(-1,1)*rotation
    }
    
    ntrials <- max(data[[group]][['gradual']]$trial)
    
    for (condition in c('abrupt','gradual')) {
      
      if (version == 1 & condition=='gradual') {
        # skip making a new panel
      } else {
        plot(-1000,-1000,
             main=group_label,xlab='',ylab='',
             xlim=c(0,ntrials+1),
             ylim=ylimits+(c(-0.1,0.1)*rotation),
             bty='n',
             ax=F)
      }
      
      cond_idx <- which(style$condition == condition)
      
      df <- data[[group]][[condition]]
      
      cdf <- aggregate(reachdeviation_deg ~ trial, data=df, FUN=getConfidenceInterval, method='t')
      
      pX <- c(c(1:ntrials),rev(c(1:ntrials)))
      pY <- c(cdf$reachdeviation_deg[,1], rev(cdf$reachdeviation_deg[,2]))
      
      polygon(x=pX,y=pY,border=NA,col=as.character(style$color_t[cond_idx]))
      
      adf <- aggregate(reachdeviation_deg ~ trial, data=df, FUN=median)
      lines(adf$reachdeviation_deg, col=as.character(style$color_s[cond_idx]))
      
      if (version == 1 & condition=='gradual') {
        # skip making a new panel
      } else {
        axis(side=1, at=c(1,ntrials))
        axis(side=2, at=c(0,ylimits))
      }
      
    }
    
    blocks <- c('rotated','clamped')
    
    plot(-1000,-1000,
         main='',xlab='',ylab='',
         xlim=c(0,(3*length(blocks))+1),
         ylim=c(c(-0.1,1.1)*rotation),
         bty='n', ax=F)
    
    for (blockno in c(1:length(blocks))) {
      
      block <- blocks[blockno]
      
      for (conditionno in c(1,2)) {
        
        condition <- c('abrupt','gradual')[conditionno]
        
        cond_idx <- which(style$condition == condition)
        
        bdf <- blocked[[group]][[condition]]
        bd <- bdf[which(bdf$block == block),]
        mu <- mean(bd$reachdeviation_deg, na.rm=TRUE)
        CI <- getConfidenceInterval(data=bd$reachdeviation_deg, method='b')
        
        xo <- ((blockno-1)*3) + conditionno
        
        polygon(x=c(c(0,1,1,0)+xo),
                y=rep(CI,each=2),
                col=as.character(style$color_t[cond_idx]),
                border=NA)
        
        lines(x=c(0,1)+xo,
              y=rep(mu,2),
              col=as.character(style$color_s[cond_idx])
              )
        
      }
      
    }
    
    axis(side=1, at=(c(0,1)*3)+2, labels = blocks)
    axis(side=2, at=c(0,rotation))
    
  }
  
  if (target %in% c('pdf')) {
    dev.off()
  }
  
}