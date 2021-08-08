
# this returns parsed data
# using the parse function specified
# for every group and participant present

parseData <- function(data, FUN) {
  
  
  # data has to either be a list of data frames
  # or a data frame
  if (is.data.frame(data)) {
    # if it is a data frame, we process it
    return(FUN(data))
  }
  
  if (is.list(data)) {
    # if it is a list, we loop through it's entries:
    
    for (itemno in c(1:length(data))) {
      
      data[[itemno]] <- parseData(data=data[[itemno]], FUN=FUN)
      
    }
    
    return(data)
  }
  
  # this should never happen:
  print('can not parse this data, returning unchanged')
  return(data)
  
}

block <- function(df) {
  
  # we want to loop through participants:
  participants <- unique(df$participant)
  
  # we collect blocked data here:
  blocked <- NA
  
  # what are the blocks?
  ntrials <- length(which(df$participant == participants[1]))
  
  if (ntrials == 164) {
    rotated_idx <- c(123:132) # last ten?
    counter_idx <- c(141:144) # last four?
    clamped_idx <- c(155:164) # last ten?
  }
  if (ntrials == 220) {
    rotated_idx <- c(151:160) # last ten?
    counter_idx <- c(177:180) # last four?
    clamped_idx <- c(211:220) # last ten?
  }
  
  # loop through participants
  for (participant in participants) {
    
    # get the participants data
    pp_idx <- which(df$participant == participant)
    pp_df <- df[pp_idx,]
    
    # get their blocked data...
    rotated <- mean(pp_df$reachdeviation_deg[rotated_idx], na.rm=TRUE)
    counter <- mean(pp_df$reachdeviation_deg[counter_idx], na.rm=TRUE)
    clamped <- mean(pp_df$reachdeviation_deg[clamped_idx], na.rm=TRUE)
    
    pp_block <- data.frame( block=c('rotated', 'counter', 'clamped'),
                            reachdeviation_deg=c(rotated, counter, clamped))
    pp_block$participant <- participant
    
    if (is.data.frame(blocked)) {
      blocked <- rbind(blocked, pp_block)
    } else {
      blocked <- pp_block
    }
    
  }
  
  return(blocked)
  
}


baseline <- function(df) {
  
  # we want to loop through participants and targets:
  participants <- unique(df$participant)
  
  # loop through participants:
  for (participant in participants) {
    
    # get the participants' data
    pp_idx <- which(df$participant == participant) 
    pp_df <- df[pp_idx,]
    
    # targets for the participant:
    targets <- unique(pp_df$targetangle_deg)
    
    if (dim(pp_df)[1] == 164) {
      baseline_idx <- c(17:32)
    }
    if (dim(pp_df)[1] == 220) {
      baseline_idx <- c(21:40)
    }
    
    for (target in targets) {
      
      # which trials correspond to this target
      target_idx <- which(pp_df$targetangle_deg == target)
      
      # this is the average reach deviation for the target during aligned phase:
      baseline <- mean( pp_df$reachdeviation_deg[ intersect(target_idx, baseline_idx) ], na.rm=TRUE )
      
      # subtract from all reaches to the target:
      pp_df$reachdeviation_deg[target_idx] <- pp_df$reachdeviation_deg[target_idx] - baseline
      
    }
    
    # put back in main data frame
    df$reachdeviation_deg[pp_idx] <- pp_df$reachdeviation_deg
    
  }
  
  # return baselined data frame:
  return(df)
  
}

