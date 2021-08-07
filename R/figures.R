
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

plotExpPerformance <- function(target='inline', exp) {
  
  info <- getInfo()
  
  groups <- info$group[which(info$experiment == exp)]
  
  
  
  
  
}