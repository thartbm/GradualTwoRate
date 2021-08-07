

# there will be 2 ways to select:

# 1: remove out of bounds trials
# 2: remove non-learning participants

getSelectedGroupsData <- function(groups=c()) {
  
  info <- getInfo()
  
  selectedData <- list()
  
  for (group in groups) {
    
    rotation <- info$rotation[which(info$group == group)]
    window <- info$window[which(info$group == group)]
    selectedData[[group]] <- list()
    
    # participants have to be selected based on performance in both conditions:
    for (condition in c('abrupt', 'gradual')) {
      
      df <- read.csv(sprintf('data/%s/%s_%s.csv', group, group, condition), stringsAsFactors = FALSE)
      
      df <- selectReaches(df=df, window=window)
      
      df <- selectParticipants(df=df, rotation=rotation)
      
      selectedData[[group]][[condition]] <- df
      
    }
    
    abrupt_participants  <- unique(selectedData[[group]][['abrupt' ]]$participant)
    gradual_participants <- unique(selectedData[[group]][['gradual']]$participant)
    
    pp_select <- intersect( gradual_participants, abrupt_participants)
    
    selectedData[[group]][['abrupt' ]] <- selectedData[[group]][['abrupt' ]][which(selectedData[[group]][['abrupt' ]]$participant %in% pp_select),]
    selectedData[[group]][['gradual']] <- selectedData[[group]][['gradual']][which(selectedData[[group]][['gradual']]$participant %in% pp_select),]
  }
  
  return(selectedData)
  
}


selectReaches <- function(df, window) {
  
  # reaches have to fall within the window:
  minimum <- rep(-1*window, dim(df)[1])
  maximum <- rep( 1*window, dim(df)[1])
  
  # but we want to adjust that for the rotation:
  rot <- df$rotation_deg * -1 # counter!
  rot[which(is.na(rot))] <- 0
  
  # adjust window down for positive rotations:
  min_idx <- which(rot < 0)
  minimum[min_idx] <- minimum[min_idx] + rot[min_idx]
  
  # adjust window up for negative rotations:
  max_idx <- which(rot > 0)
  maximum[max_idx] <- maximum[max_idx] + rot[max_idx]
  
  # points outside window:
  lo_idx <- which(df$reachdeviation_deg < minimum)
  hi_idx <- which(df$reachdeviation_deg > maximum)
  
  # set those points to NA
  df$reachdeviation_deg[c(lo_idx,hi_idx)] <- NA
  
  # return the result:
  return(df)
  
}

selectParticipants <- function(df, rotation) {
  
  # we want to look at individual participants, so let's just loop through them
  # first we need to know who they are:
  participants <- unique(df$participant)
  # and we will want to return a data frame with the accepted participants only:
  selected <- NA
  
  # here we loop through participants:
  for (participant in participants) {
    
    # we get the data for only that participant:
    pp_idx <- which(df$participant == participant)
    pp_df <- df[pp_idx,]
    
    # and index the trials for the first rotation:
    rot_idx <- sort(which(pp_df$rotation_deg == (-1 * rotation)), decreasing = TRUE)[1:40]
    
    # this is the average adaptation:
    adaptation <- mean(pp_df$reachdeviation_deg[rot_idx], na.rm=TRUE)
    
    if (adaptation > (rotation*(2/3))) {
      
      if (is.data.frame(selected)) {
        selected <- rbind(selected, pp_df)
      } else {
        selected <- pp_df
      }
      
    }
    
  }
  
  #print( setdiff( unique(df$participant), unique(selected$participant) )  )
  
  return(selected)
  
}