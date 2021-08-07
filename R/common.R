
library('osfr')


# this provides scripts with info

getInfo <- function() {
  
  group <- c('young30', 'young60', 'mcebat45', 'older45', 'young45')
  
  label <- c('younger 30', 'younger 60', 'mild cerebellar ataxia', 'older 45', 'younger 45')
  
  rotation <- c(30, 60, 45, 45, 45)
  
  experiment <- c(1, 1, 2, 2, 2)
  
  
  return(data.frame(group, label, rotation, experiment))
  
  
}

# by default, this function:
# - downloads data of all 5 groups from OSF
# - but does not overwrite it if there already is data, and
# - does not remove any downloaded zip files
# for now, the OSF repo is private, so you need a personal access token set in the OSF_PAT environment variable

downloadOSFdata <- function(groups=getInfo()$group, overwrite=FALSE, removezip=FALSE) {
  
  if (overwrite) {
    conflicts = 'overwrite'
  } else {
    conflicts = 'skip'
  }
  
  osfr::osf_auth(Sys.getenv("OSF_PAT"))
  
  OSFnode <- osfr::osf_retrieve_node("c5ezv")
  
  # get a list of files for the year and semester that is requested:
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
        # no unique file found: aborting this one
        if (length(idx) == 0) {
          return(FALSE)
        }
        return(NULL)
      }
      
      # download the file:
      if (!file.exists(sprintf('data/%s',files$name[idx])) | overwrite) {
        osfr::osf_download(x = files[idx,], 
                           path = 'data/', 
                           conflicts=conflicts)
      }
      
      # check if it is a zip file:
      if (grepl('\\.zip$', filename)) {
        
        # then unzip it there:
        unzip(sprintf('data/%s',files$name[idx]), exdir='data/')
        
        # and remove the zip file, if that is wanted:
        if (removezip) {
          file.remove(sprintf('data/%s',files$name[idx]))
        }
        
      }
      
    }
    
  }
  
}