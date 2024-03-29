library('svglite')

initialPerturbationSVGs <- function() {
  
  a <- seq(from=0, 
           to=2*pi, 
           length.out = 91)
  dot <- list()
  dot$x <- cos(a)
  dot$y <- sin(a)
  
  # 2 workspaces with:
  # targets + rotated reaches?
  
  # what size is each workspace?? 4:3? sure
  svglite(filename='doc/methods1_exp1_workspace.svg',
          width=4, height=3,
          bg="#999999",
          fix_text_size = FALSE)
  
  par(mai=c(0,0,0,0))
  
  home_target_distance = 10.5
  cursor_radius = 0.72/2 # scaled screen so it matches real tablet size
  
  # just make it all look the same:
  home_target_distance = 12
  cursor_radius = 1/2
  rotations <- (c(60,30,-30,-60)/180)*pi
  
  plot(-1000,-1000,
       main='',xlab='',ylab='',
       ylim=c(-5,16),xlim=c(-14,14),
       bty='n',ax=F,asp=1)
  
  target_degs <- c(40,50,130,140)
  target_rads <- (target_degs/180)*pi
  
  targetx <- c(0,cos(target_rads)*home_target_distance)
  targety <- c(0,sin(target_rads)*home_target_distance)
  
  for (t_idx in c(1:5)) {
    
    if (t_idx == 4) {
      polygon(x=(dot$x*cursor_radius)+targetx[t_idx],
              y=(dot$y*cursor_radius)+targety[t_idx],
              border=NA,
              col='#FFFFFF')
      
    }# else {
      lines(x=(dot$x*cursor_radius)+targetx[t_idx],
            y=(dot$y*cursor_radius)+targety[t_idx],
            #border=NA,
            col='#FFFFFF')
    #}
  }
  
  text(0,-2,'home')
  for (td in target_degs) {
    text(x=cos(pi*(td/180))*1.2*home_target_distance,
         y=sin(pi*(td/180))*1.2*home_target_distance,
         sprintf('%d°',td))
  }
  
  
  rot_target <- target_rads[3]
  lines(x = cos(rot_target)*(c(0.1,0.9)*home_target_distance),
        y = sin(rot_target)*(c(0.1,0.9)*home_target_distance),
        col='#FF0000')
  polygon(x = cos(rot_target+(c(0,0.035,0,-0.035,0)))*(c(0.92,0.85,0.87,0.85,0.92)*home_target_distance),
          y = sin(rot_target+(c(0,0.035,0,-0.035,0)))*(c(0.92,0.85,0.87,0.85,0.92)*home_target_distance),
          col='#FF0000', border=NA)
  
  for (rot in rotations) {
    lines(x = cos(rot_target+rot)*(c(0.1,0.9)*home_target_distance),
          y = sin(rot_target+rot)*(c(0.1,0.9)*home_target_distance),
          col='#FF0000',
          lty=2)
    polygon(x = cos(rot_target+(rot+c(0,0.035,0,-0.035,0)))*(c(0.92,0.85,0.87,0.85,0.92)*home_target_distance),
            y = sin(rot_target+(rot+c(0,0.035,0,-0.035,0)))*(c(0.92,0.85,0.87,0.85,0.92)*home_target_distance),
            col='#FF0000', border=NA)
    
  }
  
  
  dev.off()
  
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  svglite(filename='doc/methods1_exp2_workspace.svg',
          width=4, height=3,
          bg="#999999",
          fix_text_size = FALSE)
  
  par(mai=c(0,0,0,0))
  
  home_target_distance = 12
  cursor_radius = 1/2
  rotations <- (c(45)/180)*pi
  
  plot(-1000,-1000,
       main='',xlab='',ylab='',
       ylim=c(-5,16),xlim=c(-14,14),
       bty='n',ax=F,asp=1)
  
  target_degs <- c(50,60,120,130)
  target_rads <- (target_degs/180)*pi
  
  targetx <- c(0,cos(target_rads)*home_target_distance)
  targety <- c(0,sin(target_rads)*home_target_distance)
  
  for (t_idx in c(1:5)) {
    
    if (t_idx == 4) {
      polygon(x=(dot$x*cursor_radius)+targetx[t_idx],
              y=(dot$y*cursor_radius)+targety[t_idx],
              border=NA,
              col='#FFFFFF')
      
    }# else {
      lines(x=(dot$x*cursor_radius)+targetx[t_idx],
            y=(dot$y*cursor_radius)+targety[t_idx],
            #border=NA,
            col='#FFFFFF')
    #}    
  }
  text(0,-2,'home')
  for (td in target_degs) {
    text(x=cos(pi*(td/180))*1.2*home_target_distance,
         y=sin(pi*(td/180))*1.2*home_target_distance,
         sprintf('%d°',td))
  }
  
  rot_target <- target_rads[3]
  lines(x = cos(rot_target)*(c(0.1,0.9)*home_target_distance),
        y = sin(rot_target)*(c(0.1,0.9)*home_target_distance),
        col='#FF0000')
  polygon(x = cos(rot_target+(c(0,0.035,0,-0.035,0)))*(c(0.92,0.85,0.87,0.85,0.92)*home_target_distance),
          y = sin(rot_target+(c(0,0.035,0,-0.035,0)))*(c(0.92,0.85,0.87,0.85,0.92)*home_target_distance),
          col='#FF0000', border=NA)
  for (rot in rotations) {
    lines(x = cos(rot_target+rot)*(c(0.1,0.9)*home_target_distance),
          y = sin(rot_target+rot)*(c(0.1,0.9)*home_target_distance),
          col='#FF0000',
          lty=2)
    # x <- cos(rot_target+(rot+c(0,0.035,0,-0.035,0)))*(c(0.92,0.85,0.87,0.85,0.92)*home_target_distance)
    # x <- c(x, x[1])
    # y <- sin(rot_target+(rot+c(0,0.035,0,-0.035,0)))*(c(0.92,0.85,0.87,0.85,0.92)*home_target_distance)
    # y <- c(y, y[1])
    polygon(x = cos(rot_target+(rot+c(0,0.035,0,-0.035,0)))*(c(0.92,0.85,0.87,0.85,0.92)*home_target_distance),
            y = sin(rot_target+(rot+c(0,0.035,0,-0.035,0)))*(c(0.92,0.85,0.87,0.85,0.92)*home_target_distance),
            col='#FF0000', border=NA)
    
  }
  
  
  dev.off()
  
}

initialScheduleModelSVGs <- function() {
  
  # 2 schedules + example fit?
  style <- getStyle()
  
  svglite(filename='doc/methods2_exp1_schedules.svg',
          width=7, height=3,
          fix_text_size = FALSE)
  
  par(mar=c(4,3,0.1,0.1))
  
  sched <- getPlotSchedule(group='young30', condition='abrupt')
  plot(sched$x[1:7],
       sched$y[1:7],
       type='l',
       xlim=c(0,165),ylim=c(-1.2,1.2),
       main='',xlab='',ylab='',
       col='black',lty=1,
       bty='n',ax=F)
  
  lines(sched$x[7:8],
        sched$y[7:8],
         col='black',lty=3)
  
  for (idx in c(1:4)) {
    
    lab <- c('aligned','rotated','opposed','clamped')[idx]
    x <- mean(sched$x[(((idx-1)*2)+1):(idx*2)])
    y <- c(0.15, 1.15, -1.15, 0.15)[idx]
    text(x,y,lab)
    
  }
  
  tickx <- unique(sched$x)
  tickx[2:4] <- tickx[2:4]-1
  
  sched <- getPlotSchedule(group='young30', condition='gradual')
  
  lines(sched$x[2:3],
        sched$y[2:3],
        col='black',lty=2)
  
  
  axis(1, at=tickx, cex.axis=0.8)
  axis(2, at=c(-1,0,1), labels = c('± 30°/60°', '0°', '± 30°/60°'), cex.axis=0.8)
  
  dev.off()
  
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  svglite(filename='doc/methods2_exp2_schedules.svg',
          width=7, height=3,
          fix_text_size = FALSE)
  
  par(mar=c(4,3,0.1,0.1))
  
  sched <- getPlotSchedule(group='young45', condition='abrupt')
  plot(sched$x[1:7],
       sched$y[1:7],
       type='l',
       xlim=c(0,221),ylim=c(-.2,1.2),
       main='',xlab='',ylab='',
       col='black',lty=1,
       bty='n',ax=F)
  
  lines(sched$x[7:8],
        sched$y[7:8],
        col='black',lty=3)
  
  for (idx in c(1:4)) {
    
    lab <- c('aligned','rotated','washout','clamped')[idx]
    x <- mean(sched$x[(((idx-1)*2)+1):(idx*2)])
    y <- c(0.15, 1.15, -0.15, 0.15)[idx]
    text(x,y,lab)
    
  }
  
  tickx <- unique(sched$x)
  tickx[2:4] <- tickx[2:4]-1
  
  sched <- getPlotSchedule(group='young45', condition='gradual')
  lines(sched$x[2:3],
        sched$y[2:3],
        col='black',lty=2)
  
  axis(1, at=tickx, cex.axis=0.8)
  axis(2, at=c(0,1), labels = c('0°', '± 45°'), cex.axis=0.8)
  
  
  dev.off()
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  # get schedule for exp1:
  y30dat <- getSelectedGroupsData(groups=c('young30'))
  y30dat <- y30dat[['young30']][['abrupt']]
  participants <- unique(y30dat$participant)
  schedule <- y30dat$rotation_deg[which(y30dat$participant == participants[1])]
  
  # SMCL example data fit:
  pars <- c('Ls'=0.07305,
            'Lf'=0.42132,
            'Rs'=0.99859,
            'Rf'=0.68611)
  
  fit <- twoRateModel(par=pars, schedule=schedule)
  
  ai <- which(style$condition == 'abrupt')
  
  
  svglite(filename='doc/methods2_model_fit.svg',
          width=7, height=3,
          fix_text_size = FALSE)
  
  par(mar=c(4,3,0.1,0.1))
  
  sched <- getPlotSchedule(group='young30', condition='abrupt')
  
  plot(sched$x,
       sched$y*30,
       type='l',
       xlim=c(0,165),ylim=30*c(-1.1,1.1),
       main='',xlab='',ylab='',
       col='black',lty=1,
       bty='n',ax=F)
  
  lines(fit$total, lty=1, col=style$color_s[ai])
  lines(fit$slow, lty=2, col=style$color_s[ai])
  lines(fit$fast, lty=3, col=style$color_s[ai])
  
  legend(1,-5,legend=c('perturbation','total output','slow process','fast process'),col=c('#999999',rep(style$color_s[ai],3)),lty=c(1,1,2,3),bty='n')
  
  dev.off()
  
}


# Convert Methods Figures -----

library('rsvg')
library('magick')

rePlotMethodsFigX <- function(target='inline', fig=1) {
  
  if (target=='svg') {
    # nothing to do here:
    return()
  }
  
  dpi <- c(300,1200)[fig] #?????????????????????????????????????????????????????????????????
  
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

plotExpBehavior <- function(target='inline', exp, version=3) {
  
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
    fw_i <- 6
    fh_i <- 7.5
  }
  
  if (target=='pdf') {
    pdf(file=sprintf('doc/Fig%d.pdf',(exp*3)-1), width=fw_i, height=fh_i)
  }
  if (target=='png') {
    png(filename=sprintf('doc/Fig%d.png',(exp*3)-1), 
        width=fw_i, height=fh_i, units="in",
        res=200,
        type="cairo")
  }
  if (target == 'svg') {
    svglite::svglite( sprintf('doc/Fig%d.svg', (exp*3)-1),
                      width=fw_i, height=fh_i,
                      fix_text_size = FALSE)
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
  
  par(mar=c(4.5,4,1,0.2))
  
  for (group_no in c(1:length(names(data)))) {
    
    group <- names(data)[group_no]
    
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
        
        # add gradual schedule to figure?
        plotschedule <- getPlotSchedule(group=group, condition=condition)
        lines(x=plotschedule$x[2:3], 
              y=plotschedule$y[2:3]*rotation,
              lty=1,
              col='#999999')
      } else {
        plot(-1000,-1000,
             main=group_label,xlab='',ylab='',
             xlim=c(0,ntrials+1),
             ylim=ylimits+(c(-0.1,0.1)*rotation),
             bty='n',
             ax=F)
        title(xlab='trial', line=2.25)
        title(ylab='reach deviation [°]', line=2.5)
        # add abrupt schedule to figure?
        plotschedule <- getPlotSchedule(group=group, condition=condition)
        lines(x=plotschedule$x[1:7], 
              y=plotschedule$y[1:7]*rotation,
              lty=1,
              col='#999999')
        lines(x=plotschedule$x[7:8], 
              y=plotschedule$y[7:8]*rotation,
              lty=3,
              col='#999999')
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
        if (group_no == 1) {
          legend( x = c(  0, 180)[exp],
                  y = c( 30,  45)[exp],
                  legend = c('abrupt','gradual'),
                  bty='n',
                  lty=1,
                  col=style$color_s)
        }
      } else {
        axis(side=1, at=c(1,ntrials))
        axis(side=2, at=c(0,ylimits))
      }
      
    }
    
    blocks <- c('rotated','clamped')
    
    plot(-1000,-1000,
         main='',xlab='',ylab='',
         xlim=c(0.5,(3*length(blocks))+0.5),
         ylim=c(c(-0.1,1.1)*rotation),
         bty='n', ax=F)
    title(xlab='block', line=2.25)
    title(ylab='reach deviation [°]', line=2.5)
    title(main=group_label, font.main=1.8, cex.main=1.2)
    
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
  
  if (target %in% c('pdf','png','svg')) {
    dev.off()
  }
  
}

# model plots -----

plotExpModelFits <- function(target='inline', exp, version=4) {
  
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
  if (version == 4) {
    fh_i <- 3 * (length(groups))
  }
  
  if (target=='pdf') {
    pdf(file=sprintf('doc/Fig%d.pdf',0+(exp*3)), width=fw_i, height=fh_i)
  }
  if (target=='png') {
    png(filename=sprintf('doc/Fig%d.png',0+(exp*3)), 
        width=fw_i, height=fh_i, units="in",
        pointsize=12,
        res=250,
        type="cairo")
  }
  if (target=='svg') {
    svglite::svglite( filename = sprintf('doc/Fig%d.svg',0+(exp*3)),
                      width = fw_i, height = fh_i,
                      fix_text_size = FALSE)
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
  if (version == 4) {
    m <- c(1:ngroups)
    layout( mat = matrix(data=m,
                         nrow=ngroups,
                         byrow=TRUE) )
  }
  
  par(mar=c(4.5,4,1,0.2))
  
  for (group_no in c(1:length(names(data)))) {
    
    group <- names(data)[group_no]
    
    fits <- read.csv(sprintf('data/%s/twoRateFits.csv',group), stringsAsFactors = FALSE)
    
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
      
      fit <- c('Ls'=0,'Lf'=0,'Rs'=0,'Rf'=0)
      for (parameter in names(fit)) {
        column <- sprintf('%s.%s',condition,parameter)
        fit[parameter] <- fits[group_idx <- which(fits$participant == -1),column]
      }
      
      if (version %in% c(1,3,4) & condition=='gradual') {
        # skip making a new panel
        
        # add gradual schedule to figure?
        # add gradual schedule to figure?
        plotschedule <- getPlotSchedule(group=group, condition=condition)
        lines(x=plotschedule$x[2:3], 
              y=plotschedule$y[2:3]*rotation,
              lty=1,
              col='#999999')
        title(xlab='trial', line=2.25)
        title(ylab='reach deviation [°]', line=2.5)
      } else {
        plot(-1000,-1000,
             main=group_label,xlab='',ylab='',
             xlim=c(0,ntrials+1),
             ylim=ylimits+(c(-0.1,0.1)*rotation),
             bty='n',
             ax=F)
        # add abrupt schedule to figure?
        plotschedule <- getPlotSchedule(group=group, condition=condition)
        lines(x=plotschedule$x[1:7], 
              y=plotschedule$y[1:7]*rotation,
              lty=1,
              col='#999999')
        lines(x=plotschedule$x[7:8], 
              y=plotschedule$y[7:8]*rotation,
              lty=3,
              col='#999999')
      }
      
      cond_idx <- which(style$condition == condition)
      
      df <- data[[group]][[condition]]
      
      cdf <- aggregate(reachdeviation_deg ~ trial, data=df, FUN=getConfidenceInterval, method='t')
      
      pX <- c(c(1:ntrials),rev(c(1:ntrials)))
      pY <- c(cdf$reachdeviation_deg[,1], rev(cdf$reachdeviation_deg[,2]))
      
      polygon(x=pX,y=pY,border=NA,col=as.character(style$color_t[cond_idx]))
      
      adf <- aggregate(reachdeviation_deg ~ trial, data=df, FUN=mean, na.rm=TRUE)
      
      schedule <- df$rotation_deg[which(df$participant == df$participant[1])]
      
      # fit <- twoRateFit(schedule=schedule, 
      #                   reaches=adf$reachdeviation_deg, 
      #                   gridpoints=14, 
      #                   gridfits=10, 
      #                   checkStability=TRUE)
      
      model <- twoRateModel(par=fit,
                            schedule=schedule)
      lines(model$total, col=as.character(style$color_s[cond_idx]), lty=1)
      lines(model$slow, col=as.character(style$color_s[cond_idx]), lty=2)
      lines(model$fast, col=as.character(style$color_s[cond_idx]), lty=3)
      
      if (version %in% c(1,3,4) & condition=='gradual') {
        # skip making a new panel
        if (group_no == 1) {
          legend( x = c(  0, 180)[exp],
                  y = c( 30,  45)[exp],
                  legend = c('abrupt','gradual'),
                  bty='n',
                  lty=1,
                  col=style$color_s,
                  cex = c(0.8, 1.2)[exp] )
          legend( x = c(  0,   0)[exp],
                  y = c(-10,  45)[exp],
                  legend = c('slow','fast','total'),
                  bty='n',
                  lty=c(2,3,1),
                  col='#999999',
                  cex = c(0.8, 1.2)[exp] )
          
        }
      } else {
        axis(side=1, at=c(1,ntrials))
        axis(side=2, at=c(0,ylimits))
      }
      
    }
    
    if (version != 4) {
      plot(-1000,-1000,
           main='',xlab='',ylab='',
           xlim=c(-5,35),
           ylim=c(0,0.45),
           bty='n', ax=F)
      
      MSEs <- read.csv(sprintf('data/%s/modelMSEs.csv', group), stringsAsFactors = F)
      
      X <- seq(-5,35,length.out=301)
      
      for (cond_idx in c(1,2)) {
        condition <- c('abrupt','gradual')[cond_idx]
        MSE <- as.numeric(MSEs[,sprintf('%sfits2%sdata.MSE',substr(condition,1,2),substr(condition,1,2))])
        Y <- density(x = MSE,
                     kernel='gaussian',
                     #width=2,
                     n=301,
                     from=-5,
                     to=35)$y
        
        
        #lines(X,Y,col=as.character(style$color_s[cond_idx]))
        polX <- c(X,35,-5)
        polY <- c(Y,0,0)
        polygon(polX,polY,border=NA,col=as.character(style$color_t[cond_idx]))
        
      }
      
      for (cond_idx in c(1,2)) {
        condition <- c('abrupt','gradual')[cond_idx]
        MSE <- as.numeric(MSEs[,sprintf('abfits2%sdata.MSE',substr(condition,1,2))])
        Y <- density(x = MSE,
                     kernel='gaussian',
                     #width=2,
                     n=301,
                     from=-5,
                     to=35)$y
        
        
        lines(X,Y,col=as.character(style$color_s[cond_idx]))
        
      }
      
      MSEd <- as.numeric(MSEs[,'abfits2grdata.MSE']) - as.numeric(MSEs[,'grfits2grdata.MSE'])
      Y <- density(x = MSEd,
                   kernel='gaussian',
                   #width=2,
                   n=301,
                   from=-5,
                   to=55)$y
      
      #print(quantile(MSEd,probs=c(0.025,0.975)))
      
      lines(X,Y,col='black')
      
      axis(side=1, at=c(0,10,20,30))
      #axis(side=2, at=c(0,rotation))
      
    }
    
  }
  
  if (target %in% c('pdf','png','svg')) {
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
    pdf(file=sprintf('doc/Fig%d.pdf',1+(exp*3)), width=fw_i, height=fh_i)
  }
  if (target=='png') {
    png(filename=sprintf('doc/Fig%d.png',1+(exp*3)), 
        width=fw_i, height=fh_i, units="in",
        res=200,
        type="cairo")
  }
  if (target=='svg') {
    svglite::svglite( filename = sprintf('doc/Fig%d.svg',1+(exp*3)),
                      width = fw_i, height = fh_i,
                      fix_text_size = FALSE)
  }
  
  
  
  ngroups = length(groups)
  
  if (version == 1) {
    layout( mat=matrix(data=c(1:(ngroups*4)), nrow = ngroups, ncol=4, byrow = TRUE) )
  }

  par(mar=c(3.5,2,1,0.25))
  
  inset.figs <- list()
  inset.fig.idx <- 1
  
  for (group_no in c(1:length(groups))) {
    
    group <- groups[group_no]
    
    group_idx <- which(info$group == group)
    group_label <- info$label[group_idx]
    
    fits <- read.csv(sprintf('data/%s/bootstrapped_TwoRateFits.csv',group),stringsAsFactors = F)
    
    for (parameter in c('Ls','Rs','Lf','Rf')) {
      
      ylims <- list(  'Ls'=c(-8,40),
                      'Rs'=c(-8,40),
                      'Lf'=c(-3,20),
                      'Rf'=c(-3,20)  )[[parameter]]
      
      # position of inset plots:
      p <- list( 'Ls'=c(.50, .90, .70, .95),
                 'Rs'=c(.10, .50, .70, .95),
                 'Lf'=c(.50, .90, .70, .95),
                 'Rf'=c(.10, .50, .70, .95)  )[[parameter]]
      
      main <- ''
      ylab <- ''
      xlab <- ''
      if (group == groups[1]) {
        main <- parameter
      }

      plot(-1000,-1000,
           main=main,xlab=xlab,ylab=ylab,
           xlim=c(0,1),
           ylim=ylims,
           bty='n',
           ax=F)
      if (parameter == 'Ls') {
        title( ylab = info$label[which(info$group == group)], line=1 )
      }
      
      groupfits <- read.csv(sprintf('data/%s/twoRateFits.csv',group), stringsAsFactors=TRUE)
      
      for (condition in c('abrupt','gradual')) {
        
        column  <- sprintf('%s.%s',condition,parameter)
        parvals <- groupfits[which(groupfits$participant > -1),column]
        
        cond_idx <- which(style$condition == condition)
        
        points(x=parvals,
               y=rep(cond_idx/3*ylims[1],length(parvals)),
               col=as.character(style$color_t[cond_idx]),
               pch=16, cex=1.0)
        
        colname <- sprintf('%s.%s',condition,parameter)
        parvals <- as.numeric(unlist(fits[colname]))
        
        X <- seq(0,1,length.out = 201)
        Y <- density(parvals,
                     n=201,from=0,to=1,
                     na.rm=T,kernel='gaussian',width=0.05)$y
        
        lines(X,Y,col=as.character(style$color_s[cond_idx]))
        
      }
      
      axis(side=1,at=c(0,1))
      
      if (group_no == 1 & parameter == 'Ls') {
        legend( x = 0.55,
                y = 18,
                legend = c('abrupt','gradual'),
                lty = 1, bty='n', cex=0.8,
                col=style$color_s)
      }
      
      if (group_no == length(groups)) {
        title(xlab = 'parameter value', line=2)
      }
      
      
      inset.figs[[inset.fig.idx]] <- c(grconvertX(p[1:2], from="npc", to="ndc"),
                                       grconvertY(p[3:4], from="npc", to="ndc"))
      inset.fig.idx <- inset.fig.idx + 1
      
    }

  }
  
  inset.fig.idx <- 1
  
  
  recovered_fits <- list()
  recovered_fits[['abrupt']] <- read.csv('data/young45/recovered_fits_zero_abrupt.csv', stringsAsFactors=FALSE)
  recovered_fits[['ramped']] <- read.csv('data/young45/recovered_fits_zero_ramped.csv', stringsAsFactors=FALSE)
  
  recovered_CIs <- list()
  for (parameter in c('Ls','Rs','Lf','Rf')) {
    pardiffs <- as.numeric(unlist(recovered_fits[['ramped']][[parameter]] - recovered_fits[['abrupt']][[parameter]]))
    recovered_CIs[[parameter]] <- quantile(pardiffs, 
                                           probs=c(0.025,0.975),
                                           names=FALSE)
  }
  

  for (group in groups) {
    
    group_idx <- which(info$group == group)
    group_label <- info$label[group_idx]
    
    fits <- read.csv(sprintf('data/%s/bootstrapped_TwoRateFits.csv',group),stringsAsFactors = F)
    
    for (parameter in c('Ls','Rs','Lf','Rf')) {
      
      # create inset figure
      op <- par( fig=inset.figs[[inset.fig.idx]], 
                 new=TRUE, 
                 mar=rep(0, 4), 
                 cex=0.5)
      
      xlim <- c(-.1,.1)
      if (parameter %in% c('Lf','Rf')) {
        xlim <- c(-.3,.3)
      }
      if (parameter == 'Rs') {
        xlim <- c(-.02,.02)
      }
      
      plot(-1000,-1000,
           xlim=xlim,
           ylim=c(0,3),
           bty='n', ax=F)
      
      #lines(c(0,0),c(0.25,2.75),col='black')

      CI <- recovered_CIs[[parameter]]
      polX <- c(CI,rev(CI))
      polY <- c(0.25,0.25,2.75,2.75)
      polygon(polX,polY,col='#999999',border=NA)
      
      col1 <- sprintf('abrupt.%s',parameter)
      col2 <- sprintf('gradual.%s',parameter)
      pardiffs <- as.numeric(unlist(fits[col2])) - as.numeric(unlist(fits[col1]))
      CI <- quantile(pardiffs, probs=c(0.025,0.975), names = FALSE)
      
      polX <- c(CI,rev(CI))
      polY <- c(0.6,0.6,1,1)
      polygon(polX,polY,col=as.character(style$color_t[2]),border=NA)
      
      dX <- seq(xlim[1],xlim[2],length.out = 100)
      dY <- density(pardiffs,from = xlim[1], to=xlim[2], n=100)$y
      lines(dX,(dY/max(dY))+1.2,col=as.character(style$color_s[2]))
      

      axis(side=1,at=xlim,padj=-1.5)
      axis(side=1,at=0,labels=NA)
      
      # finalize inset figure
      inset.fig.idx <- inset.fig.idx + 1
      par(op)
    }
  }
  
  if (target %in% c('pdf','svg','png')) {
    dev.off()
  }
  
}

# plot parameter recovery simulation -----

plotParameterRecovery <- function(target='inline') {
  
  info <- getInfo()
  style <- getStyle()
  
  fw_i <- 7
  fh_i <- 2
  
  if (target=='pdf') {
    pdf(file=sprintf('doc/FigA1.pdf'), width=fw_i, height=fh_i)
  }
  
  
  layout( mat = matrix(data = c(1:8), 
                       nrow = 2, 
                       ncol = 4, 
                       byrow = TRUE) )
  
  par(mar=c(3,4,1,0.2))
  
  # first we retrieve the true parameters, just like in the simulation:
  # we want ground truth parameters in this vector:
  ground_truth_parameters <- c('Ls'=NA,
                               'Lf'=NA,
                               'Rs'=NA,
                               'Rf'=NA)
  
  # we'll use the young45 group's abrupt fit for the group:
  y45fits <- read.csv('data/young45/twoRateFits.csv', 
                      stringsAsFactors = FALSE)
  
  # and put those parameters in our vector:
  group_idx <- which(y45fits$participant == -1)
  for (par in names(ground_truth_parameters)) {
    column <- sprintf('abrupt.%s',par)
    ground_truth_parameters[par] <- as.numeric(unlist(y45fits[group_idx,column]))
  }
  
  # we also load all the recovered fits:
  recovered_fits <- list('zero' <- list('abrupt'=NA,
                                        'ramped'=NA),
                         'cntr' <- list('abrupt'=NA,
                                        'ramped'=NA))
  for (counter in c('zero','cntr')) {
    for(condition in c('abrupt','ramped')) {
      df <- read.csv(file=sprintf('data/young45/recovered_fits_%s_%s.csv',counter,condition), 
                     stringsAsFactors = FALSE)
      recovered_fits[[counter]][[condition]] <- df
    }
  }
  
  # now start going through data to actually plot it:
  for (counter_no in c(1:2)) {
    
    counter <- c('zero','cntr')[counter_no]
    counter_name <- c('washout','counter')[counter_no]
    
    for (parameter in c('Ls','Rs','Lf','Rf')) {
      
      if (counter_no == 1) {
        main <- parameter
      } else {
        main <- ''
      }
      if (parameter == 'Ls') {
        ylab <- counter_name
      } else {
        ylab <- ''
      }
      
      plot(-1000,-1000,
           main=main, ylab=ylab, xlab='',
           xlim=c(0,1),
           ylim=c(0,3),
           bty='n',
           ax=F)
      
      lines(x = rep(ground_truth_parameters[parameter],2),
            y = c(0.5,2.5),
            col='black')
      
      for (condition in c('abrupt','gradual')) {
        
        condition_name <- c('abrupt'='abrupt',
                            'gradual'='ramped')[condition]
        
        style_idx <- which(style$condition == condition)
        
        df <- recovered_fits[[counter]][[condition_name]]
        parvals <- df[[parameter]]
        CI <- quantile(parvals, probs=c(0.025, 0.975), names=FALSE)
        
        for (FUN in c(min,max)) {
          extreme <- FUN(parvals,na.rm=TRUE)
          lines(x = rep(extreme,2),
                y = c(-0.4,0.4)+style_idx,
                col=as.character(style$color_t[style_idx]))
        }
        
        polX <- c(CI,rev(CI))
        polY <- rep( c(-0.4,0.4)+style_idx, each=2)
        polygon(polX,polY,border = NA, col=as.character(style$color_t[style_idx]))
        
      }
      
      axis(side=1, at=c(0,1))
      
    }
    
  }
  
  if (target %in% c('pdf')) {
    dev.off()
  }
  
}

# support functions -----

getPlotSchedule <- function(group, condition) {
  
  # 30/60 tablet groups:
  #   
  #  32 aligned
  # 100 rotated
  #  12 opposed
  #  20 clamped
  # 
  # [40 trial ramp, still 0 on first rotated trial, max on 41st trial]
  # 
  # 45 munich groups:
  #   
  #  40 aligned
  # 120 rotated
  #  20 washout
  #  40 clamped
  # 
  # [60 trial ramp, still 0 on first rotated trial, max on 61st trial]
  
  if (group %in% c('young30', 'young60')) {
    if (condition == 'abrupt') {
      plotschedule <- list( 'x' = c(1, 33, 33, 133, 133, 145, 145, 164),
                            'y' = c(0,  0,  1,   1,  -1,  -1,   0,   0)   )
    }
    if (condition == 'gradual') {
      plotschedule <- list( 'x' = c(1, 33, 73, 133, 133, 145, 145, 164),
                            'y' = c(0,  0,  1,   1,  -1,  -1,   0,   0)   )
    }
  }
  
  if (group %in% c('young45', 'older45', 'mcebat45')) {
    if (condition == 'abrupt') {
      plotschedule <- list( 'x' = c(1, 41, 41, 161, 161, 181, 181, 220),
                            'y' = c(0,  0,  1,   1,   0,   0,   0,   0)   )
    }
    if (condition == 'gradual') {
      plotschedule <- list( 'x' = c(1, 41, 101, 161, 161, 181, 181, 220),
                            'y' = c(0,  0,   1,   1,   0,   0,   0,   0)   )
    }
  }
  
  return(plotschedule)
  
}