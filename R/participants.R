

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


