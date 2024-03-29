
library('osfr')


# meta info -----

# this provides scripts with info

getInfo <- function() {
  
  group      <- c(   'young30',    'young60',               'mcebat45',  'older45',    'young45')
  
  label      <- c('younger 30°', 'younger 60°', 'mild cerebellar ataxia', 'older 45°', 'younger 45°')
  
  rotation   <- c(          30,           60,                       45,         45,           45)
  
  experiment <- c(           1,            1,                        2,          2,            2)
  
  window     <- c(          45,           45,                       45,         45,           45)
  
  return(data.frame(group, 
                    label, 
                    rotation, 
                    experiment, 
                    window))
  
}

getStyle <- function() {
  
  condition <- c('abrupt', 'gradual')
  
  label     <- c('abrupt', 'ramped')
  
  color_s   <- c('#ff8200ff', '#b400e4ff')
  
  color_t   <- c('#ff82002f', '#b400e42f')
  
  return(data.frame(condition,
                    label,
                    color_s,
                    color_t))
  
}

getOrder <- function(exp=NULL) {
  
  if (exp == 1) {
    
    # A: gradual, left, clockwise; abrupt, right, counterclockwise
    # B: gradual, left, counterclockwise; abrupt, right, clockwise
    # C: gradual, right, clockwise; abrupt, left, counterclockwise 
    # D: gradual, right, counterclockwise; abrupt, left, clockwise 
    # E: abrupt, left, clockwise; gradual, right, counterclockwise
    # F: abrupt, left, counterclockwise; gradual, right, clockwise
    # G: abrupt, right, clockwise; gradual, left, counterclockwise
    # H: abrupt, right, counterclockwise; gradual, left, clockwise
    
    version <- c('A','B','C','D','E','F','G','H')
    first_condition <- c('gradual','gradual','gradual','gradual','abrupt','abrupt','abrupt','abrupt')
    first_quadrant <- c(2,2,1,1,2,2,1,1)
    first_rotation <- c(-1,1,-1,1,-1,1,-1,1)
    
    return(data.frame(version,first_condition,first_quadrant,first_rotation))
    
  }
  
  if (exp==2) {
    
    version <- c(2,3,4,5)
    
    first_condition <- c('gradual','gradual','abrupt','abrupt')
    first_quadrant <- c(1,2,1,2)
    first_rotation <- c(1,-1,1,-1)
    
    return(data.frame(version,first_condition,first_quadrant,first_rotation))
    
  }
  
  warning('`exp` variable is NULL: set to 1 or 2')
  return(exp)
  
}

getVersion <- function(exp=NULL, group=NULL) {
  
  # since versions are coded as different types in the 2 experiments
  # (letter versus number: character / integer)
  # you can only get 1 group, or 1 experiment from this function
  # otherwise the integers get converted to characters
  
  info <- getInfo()
  
  if (is.null(group)) {
    groups <- info$group[which(info$experiment == exp)]
  } else {
    groups <- c(group)
  }
  
  
  df <- NA
  
  for (group in groups) {
    
    experiment <- info$experiment[which(info$group == group)]
    exp_order <- getOrder(exp=experiment)
    
    gdf <- read.csv(sprintf('data/%s_demographics.csv',group), stringsAsFactors = FALSE)
    gdf <- gdf[c('participant','version')]
    
    gdf <- merge(x=gdf, y=exp_order, by.x='version', by.y='version', all=TRUE, no.dups=FALSE)
    
    if (is.data.frame(df)) {
      df <- rbind(df,gdf)
    } else {
      df <- gdf
    }
    
  }
  
  return(df)
  
}

# download OSF data -----

# by default, this function:
# - downloads data of all 5 groups from OSF
# - but does not overwrite it if there already is data, and
# - does not remove any downloaded zip files
# for now, the OSF repo is private, so you need a personal access token set in the OSF_PAT environment variable

downloadOSFdata <- function(groups=c(getInfo()$group,'processed'), overwrite=TRUE, removezip=FALSE) {
  
  if (overwrite) {
    conflicts = 'overwrite'
  } else {
    conflicts = 'skip'
  }
  
  #osfr::osf_auth(Sys.getenv("OSF_PAT"))
  
  OSFnode <- osfr::osf_retrieve_node("c5ezv")
  
  files <- osfr::osf_ls_files(OSFnode)
  print(files)
  
  for (group in groups) {
    
    # get list of files:
    filenames <- c(sprintf('%s.zip',group))
    
    for (filename in filenames) {
      
      print(filename)
      
      # find which line corresponds to the file:
      idx <- which(files$name == filename)
      
      # check that the file exists on OSF, and is unique:
      # if not: abort
      if (length(idx) > 1) {
        return(NULL)
      }
      # no unique file found: aborting this one
      if (length(idx) == 0) {
        return(FALSE)
      }
      
      # download the file:
      if (!file.exists(file.path('data',files$name[idx])) | overwrite) {
        osfr::osf_download(x = files[idx,], 
                           path = 'data/', 
                           conflicts=conflicts)
      }
      
      # check if it is a zip file:
      if (grepl('\\.zip$', filename)) {
        
        # then unzip it there:
        unzip(file.path('data',files$name[idx]), exdir='data')
        
        # and remove the zip file, if that is wanted:
        if (removezip) {
          file.remove(file.path('data',files$name[idx]))
        }
        
      }
      
    }
    
  }
  
}

# data information -----


getAllGroupParticipants <- function(group) {
  
  df <- read.csv(sprintf('data/%s_demographics.csv',group), stringsAsFactors = FALSE)
  
  if ('participant' %in% names(df)) {
    return(df$participant)
  }
  
  if ('ID' %in% names(df)) {
    print('ID')
    return(df$ID)
  }
  
  # young30:  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 28 29 30 31
  # young60: 101 102 104 105 106 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 130 131 132
  
  # older45: 301 302 303 304 305 306 307 308 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216
  # mcebat45: 401 402 403 404 405 406 407 408 409 410 411 412 413 414 415 416 417 418 419 420 421 422
  # younger45: 500 501 502 503 504 505 506 507 508 509 510 511 512 513 514 515 516 517 518 519 520 521 522 523 524 525 526 527 528 529
  
  
}


checkGroupData <- function(group, report='list') {
  
  participants <- getAllGroupParticipants(group)
  notOK <- c()
  
  for (participant in participants) {
    
    for (condition in c('abrupt','gradual')) {
      
      filename <- sprintf('data/%s/%s_p%03d_%s.csv', group, group, participant, condition)
      
      if (!file.exists(filename)) {
        if (report == 'print') {
          print(filename)
        }
        if (report == 'list') {
          notOK <- c(notOK, participant)
        }
        
      }
      
    }
    
  }
  
  if (report == 'list') {
    if (length(notOK) > 0) {
      return(unique(notOK))
    } else {
      return(notOK)
    }
  }
  
}

getGroupParticipants <- function(group) {
  
  participants <- getAllGroupParticipants(group)
  notOK <- checkGroupData(group)
  
  participants <- setdiff(participants, notOK)
  
  return(participants)
  
}



# data selecting -----


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

getGoodGroupParticipants <- function(group) {
  
  info <- getInfo()
  rotation <- info$rotation[which(info$group == group)]
  
  gradual_pp <- unique(selectParticipants(df=read.csv(sprintf('data/%s/%s_gradual.csv',group,group), stringsAsFactors = FALSE), rotation)$participant)
  abrupt_pp <- unique(selectParticipants(df=read.csv(sprintf('data/%s/%s_abrupt.csv',group,group), stringsAsFactors = FALSE), rotation)$participant)
  
  goodParticipants <- intersect(gradual_pp, abrupt_pp)
  
  return(goodParticipants)
}

showGoodParticipants <- function() {
  
  info <- getInfo()
  
  group <- c()
  all_pp <- c()
  good_pp <- c()
  
  for (group_name in info$group) {
    
    group <- c(group, group_name)
    all_pp <- c(all_pp, length(getGroupParticipants(group=group_name)))
    good_pp <- c(good_pp, length(getGoodGroupParticipants(group=group_name)))
    
  }
  
  return(data.frame(group,all_pp,good_pp))
  
}


# data parsing -----

# this returns parsed data
# using the parse function specified
# for every group and participant present
# 
# data can be a (nested) list with data frames
# or a data frame (parseData is recursive)
#
# functions defined here:
# - block (get average reach deviations in a few specified blocks)
# - baseline (subtract average reach deviations from the second half of the aligned phase, per target)
# - addtrial (adds a column in the data frame with trial numbers)
# - addphase (adds a column in the data frame with phase numbers)



parseData <- function(data, FUN, ...) {
  
  
  # data has to either be a list of data frames
  # or a data frame
  if (is.data.frame(data)) {
    # if it is a data frame, we process it
    return(FUN(data, ...))
  }
  
  if (is.list(data)) {
    # if it is a list, we loop through it's entries:
    
    for (itemno in c(1:length(data))) {
      
      #print(names(data)[itemno])
      data[[itemno]] <- parseData(data=data[[itemno]], FUN=FUN, ...)
      
    }
    
    return(data)
  }
  
  # this should never happen:
  print('can not parse this data, returning unchanged')
  return(data)
  
}

block <- function(df, blocks=c('rotated', 'clamped')) {
  
  # by default, the average reach deviation in two blocks of trials are returned:
  # 'rotated': the last 10 trials of the initial rotation
  # 'clamped': the last 10 trials of the error-clamped phase
  # these are the ones we use in our analyses
  
  # two more are available:
  # 'early': the first 4 trials of the initial rotation (makes no sense for gradual)
  # 'counter': the last 4 trials of the counter rotation (pretty flaky, estimates)
  
  # more could be added, or the number of trials could be changed... 
  
  # we want to loop through participants:
  participants <- unique(df$participant)
  
  # we collect blocked data here:
  blocked <- NA
  
  # what are the blocks?
  ntrials <- length(which(df$participant == participants[1]))
  
  if (ntrials == 164) {
    early_idx   <- c(33:36) # first four?
    rotated_idx <- c(123:132) # last ten?
    counter_idx <- c(141:144) # last four?
    clamped_idx <- c(155:164) # last ten?
  }
  if (ntrials == 220) {
    early_idx   <- c(41:44) # first four?
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
    early   <- mean(pp_df$reachdeviation_deg[early_idx], na.rm=TRUE)
    rotated <- mean(pp_df$reachdeviation_deg[rotated_idx], na.rm=TRUE)
    counter <- mean(pp_df$reachdeviation_deg[counter_idx], na.rm=TRUE)
    clamped <- mean(pp_df$reachdeviation_deg[clamped_idx], na.rm=TRUE)
    
    reachdevs <- c('early'=early, 'rotated'=rotated, 'counter'=counter, 'clamped'=clamped)
    
    # pp_block <- data.frame( block=c('rotated', 'counter', 'clamped'),
    #                         reachdeviation_deg=c(rotated, counter, clamped))
    pp_block <- data.frame( block=blocks,
                            reachdeviation_deg=reachdevs[blocks])
    pp_block$participant <- participant
    
    if (is.data.frame(blocked)) {
      blocked <- rbind(blocked, pp_block)
    } else {
      blocked <- pp_block
    }
    
  }
  
  rownames(blocked) <- c()
  
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

addtrial <- function(df) {
  
  # we want to loop through participants and targets:
  participants <- unique(df$participant)
  
  output <- NA
  
  # loop through participants:
  for (participant in participants) {
    
    # get the participants' data
    pp_idx <- which(df$participant == participant) 
    pp_df <- df[pp_idx,]
    
    pp_df$trial <- c(1:dim(pp_df)[1])
    
    if (is.data.frame(output)) {
      output <- rbind(output, pp_df)
    } else {
      output <- pp_df
    }
    
  }
  
  # return data frame with trial numbers:
  return(output)
  
}


addphase <- function(df) {
  
  # we want to loop through participants and targets:
  participants <- unique(df$participant)
  
  output <- NA
  
  # loop through participants:
  for (participant in participants) {
    
    # get the participants' data
    pp_idx <- which(df$participant == participant) 
    pp_df <- df[pp_idx,]
    
    if (dim(pp_df)[1] == 164) {
      pp_df$phase <- c(rep(1,32), rep(2,100), rep(3,12), rep(4,20))
    }
    if (dim(pp_df)[1] == 220) {
      pp_df$phase <- c(rep(1,40), rep(2,120), rep(3,20), rep(4,40))
    }
    
    if (is.data.frame(output)) {
      output <- rbind(output, pp_df)
    } else {
      output <- pp_df
    }
    
  }
  
  # return data frame with phase added:
  return(output)
  
}

normalize <- function(df) {
  
  rotation <- max(abs(df$rotation_deg), na.rm = TRUE)
  
  df$reachdeviation_deg <- df$reachdeviation_deg / rotation
  
  return(df)
  
}

# statistical functions -----

getConfidenceInterval <- function(data, variance = var(data, na.rm=TRUE), conf.level = 0.95, method='t-distr', resamples=1000, FUN=mean) {
  
  data <- data[which(!is.na(data))]
  
  if (method %in% c('t-distr','t')) {
    
    z = qt((1 - conf.level)/2, df = length(data) - 1, lower.tail = FALSE)
    
    xbar = mean(data)
    sdx = sqrt(variance/length(data))
    
    return(c(xbar - z * sdx, xbar + z * sdx))
    
  }
  
  # add sample z-distribution?
  
  if (method %in% c('bootstrap','b')) {
    
    data <- data[which(is.finite(data))] #need is.finite due to NA values
    
    samplematrix <- matrix(sample(data, size = resamples*length(data), replace = TRUE), nrow = resamples)
    BS <- apply(samplematrix, c(1), FUN=FUN) 
    
    lo <- (1-conf.level)/2.
    hi <- 1 - lo
    
    return(quantile(BS, probs = c(lo,hi)))
    
  }
  
}

