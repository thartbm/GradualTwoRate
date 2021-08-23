
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
  if (version == 3) {
    fw_i <- 4
    fh_i <- 7.5
  }
  
  if (target=='pdf') {
    pdf(file=sprintf('doc/Fig%d.pdf',0+(exp*3)), width=fw_i, height=fh_i)
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
  if (version == 3) {
    m <- c( rep( (c(1:ngroups)*2)-1, each=ngroups), c(1:ngroups)*2)
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
      
      if (version %in% c(1,3) & condition=='gradual') {
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
      
      adf <- aggregate(reachdeviation_deg ~ trial, data=df, FUN=mean, na.rm=TRUE)
      lines(adf$reachdeviation_deg, col=as.character(style$color_s[cond_idx]))
      
      if (version %in% c(1,3) & condition=='gradual') {
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


plotExpModelFits <- function(target='inline', exp, version=1) {
  
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
  if (version == 3) {
    fw_i <- 4
    fh_i <- 7.5
  }
  
  if (target=='pdf') {
    pdf(file=sprintf('doc/Fig%d.pdf',1+(exp*3)), width=fw_i, height=fh_i)
  }
  
  data <- getSelectedGroupsData(groups=groups)
  
  data <- parseData(data, FUN=baseline)
  data <- parseData(data, FUN=addtrial)

  ngroups = length(data)
  
  if (version == 1) {
    layout( mat=matrix(data=c(1:(ngroups*2)), nrow = ngroups, ncol=2, byrow = TRUE), widths=c(1,0.5) )
  }
  if (version == 2) {
    m <- c(unlist(lapply( ((c(1:ngroups)-1) * 3), function(x){rep(c(1,2)+as.numeric(x),each=3)} )), rep(c(1:ngroups)*3, each=6/ngroups) )
    layout( mat=matrix(data=m, nrow=ngroups+1, byrow=TRUE ) )
  }
  if (version == 3) {
    m <- c( rep( (c(1:ngroups)*2)-1, each=ngroups), c(1:ngroups)*2)
    layout( mat=matrix(data=m, nrow=ngroups+1, byrow=TRUE ) )
  }
  
  par(mar=c(4,3,1,0.2))
  
  for (group in names(data)) {
    
    group_idx <- which(info$group == group)
    rotation <- info$rotation[group_idx]
    group_label <- info$label[group_idx]
    
    if (rotation == 45) {
      ylimits <- c(-15,45)
    } else {
      ylimits <- c(-1,1)*rotation
    }
    
    ntrials <- max(data[[group]][['gradual']]$trial)
    
    for (condition in c('abrupt','gradual')) {
      
      if (version %in% c(1,3) & condition=='gradual') {
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
      
      adf <- aggregate(reachdeviation_deg ~ trial, data=df, FUN=mean, na.rm=TRUE)
      
      schedule <- df$rotation_deg[which(df$participant == df$participant[1])]
      
      fit <- twoRateFit(schedule=schedule, 
                        reaches=adf$reachdeviation_deg, 
                        gridpoints=14, 
                        gridfits=10, 
                        checkStability=TRUE)
      
      model <- twoRateModel(par=fit,
                            schedule=schedule)
      lines(model$total, col=as.character(style$color_s[cond_idx]), lty=1)
      lines(model$slow, col=as.character(style$color_s[cond_idx]), lty=2)
      lines(model$fast, col=as.character(style$color_s[cond_idx]), lty=3)
      
      if (version %in% c(1,3) & condition=='gradual') {
        # skip making a new panel
      } else {
        axis(side=1, at=c(1,ntrials))
        axis(side=2, at=c(0,ylimits))
      }
      
    }
    
    plot(-1000,-1000,
         main='',xlab='',ylab='',
         xlim=c(-5,35),
         ylim=c(0,0.4),
         bty='n', ax=F)
    
    MSEs <- read.csv(sprintf('data/%s/rampedMSEfromAbruptFit.csv', group), stringsAsFactors = F)
    
    X <- seq(-5,55,length.out=301)
    
    for (cond_idx in c(1,2)) {
      condition <- c('abrupt','gradual')[cond_idx]
      MSE <- as.numeric(MSEs[,sprintf('%s.MSE',condition)])
      Y <- density(x = MSE,
                   kernel='gaussian',
                   width=2,
                   n=301,
                   from=-5,
                   to=55)$y
      
      lines(X,Y,col=as.character(style$color_s[cond_idx]))
      
    }
    
    MSEd <- as.numeric(MSEs[,'abrupt.MSE']) - as.numeric(MSEs[,'gradual.MSE'])
    Y <- density(x = MSEd,
                 kernel='gaussian',
                 width=2,
                 n=301,
                 from=-5,
                 to=55)$y
    
    lines(X,Y,col='black')
    
    axis(side=1, at=c(0,10,20,30))
    #axis(side=2, at=c(0,rotation))
    
  }
  
  if (target %in% c('pdf')) {
    dev.off()
  }
  
}


plotExpModelParameters <- function(target='inline', exp, version=1) {
  
  info <- getInfo()
  style <- getStyle()
  
  groups <- info$group[which(info$experiment == exp)]
  
  if (version == 1) {
    fw_i <- 7
    fh_i <- (7/4)*length(groups)
  }

  if (target=='pdf') {
    pdf(file=sprintf('doc/Fig%d.pdf',2+(exp*3)), width=fw_i, height=fh_i)
  }
  
  ngroups = length(groups)
  
  if (version == 1) {
    layout( mat=matrix(data=c(1:(ngroups*4)), nrow = ngroups, ncol=4, byrow = TRUE) )
  }

  par(mar=c(4,3,1,1))
  
  for (group in groups) {
    
    group_idx <- which(info$group == group)
    group_label <- info$label[group_idx]
    
    fits <- read.csv(sprintf('data/%s/bootstrapped_TwoRateFits.csv',group),stringsAsFactors = F)
    
    for (parameter in c('Ls','Rs','Lf','Rf')) {
      
      ylims <- list('Ls'=c(0,40),
                    'Rs'=c(0,40),
                    'Lf'=c(0,15),
                    'Rf'=c(0,15))[[parameter]]
      
      main<-''
      if (parameter == 'Ls') {
        main<-group
      }
      xlab<-''
      if (group == groups[length(groups)]) {
        xlab<-parameter
      }
      print(xlab)
      plot(-1000,-1000,
           main=main,xlab=xlab,ylab='',
           xlim=c(0,1),
           ylim=ylims,
           bty='n',
           ax=F)
      
      for (condition in c('abrupt','gradual')) {
        
        cond_idx <- which(style$condition == condition)
        
        colname <- sprintf('%s.%s',condition,parameter)
        parvals <- as.numeric(unlist(fits[colname]))
        
        X <- seq(0,1,length.out = 101)
        Y <- density(parvals,
                     n=101,from=0,to=1,
                     na.rm=T,kernel='gaussian',width=0.05)$y
        
        lines(X,Y,col=as.character(style$color_s[cond_idx]))
        
      }
      
      axis(side=1,at=c(0,0.5,1))
      
    }

  }
  
  if (target %in% c('pdf')) {
    dev.off()
  }
  
}