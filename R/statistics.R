
require('afex')
require('BayesFactor')


# order effects -----

testOrderEffects <- function(exp=NULL, groups=NULL, dv='lambda') {

  if (is.null(groups)) {
    # we need to know the groups in the experiment:
    info <- getInfo()
    groups <- info$group[which(info$experiment == exp)]
  }

  # we collect all data here for an ANOVA on all the data:
  df <- NA

  group_anovas <- list()

  for (group in groups) {

    # we collect group data here, so we can do group-based ANOVAs:
    g_df <- read.csv(sprintf('data/%s/abrupt_learningrates.csv',group))

    # add group variable to the data frame:
    g_df$group <- group

    # make factors:
    for (column in c('participant','first_condition','group')) {
      g_df[,column] <- as.factor(g_df[,column])
    }

    # do group ANOVA on group data:

    g_aov <- afex::aov_ez(id = "participant",
                          dv = dv,
                          data = g_df,
                          between = c("first_condition"),
                          type = 3,
    )

    group_anovas[[group]] <- g_aov


    # add group data to all data:

    if (is.data.frame(df)) {
      df <- rbind(df,g_df)
    } else {
      df <- g_df
    }

  }

  # do ANOVA on all groups:
  if (!is.null(exp)) {
    e_aov <- afex::aov_ez(id = "participant",
                          dv = dv,
                          data = df,
                          between = c("group","first_condition"),
                          type = 3 )
    
    cat(sprintf('\n===== Experiment: %d =====\n',exp ))
    print(summary(e_aov))
  } else {
    # print the group-based ANOVAs:
    for (group in names(group_anovas)) {
      cat(sprintf('\n===== %s =====\n',toupper(group) ))
      print(summary(group_anovas[[group]]))
    }
  }
  

}


# adaptation and rebounds -----

doExpANOVAs <- function(exp=NULL, groups=NULL) {
  
  if (is.null(groups)) {
    # we need to know the groups in the experiment:
    info <- getInfo()
    groups <- info$group[which(info$experiment == exp)]
  }
  
  # we get the data form these groups:
  data <- getSelectedGroupsData(groups=groups)
  
  # we pre-process it, so that it is ready for ANOVAs
  data <- parseData(data, FUN=baseline)
  data <- parseData(data, FUN=addtrial)
  data <- parseData(data, FUN=normalize)
  data <- parseData(data, FUN=block)
  
  # we collect all data here for an ANOVA on all the data:
  df <- NA
  
  group_anovas <- list()
  
  for (group in names(data)) {
    
    # we collect group data here, so we can do group-based ANOVAs:
    g_df <- NA
    
    for (condition in c('abrupt', 'gradual')) {
      
      # add group and condition variables to the data frame:
      gc_df <- data[[group]][[condition]]
      gc_df$group <- group
      gc_df$condition <- condition
      
      # collect data in group data frame:
      if (is.data.frame(g_df)) {
        g_df <- rbind(g_df, gc_df)
      } else {
        g_df <- gc_df
      }
      
    }
    
    # make factors:
    for (column in c('participant','condition','group')) {
      g_df[,column] <- as.factor(g_df[,column])
    }
    
    # do group ANOVA on group data:
    
    g_aov <- afex::aov_ez(id = "participant",
                          dv = "reachdeviation_deg",
                          data = g_df,
                          within = c("block","condition"),
                          #observed = 'reachdeviation_deg',
                          type = 3,
                          )
    
    group_anovas[[group]] <- g_aov
    
    
    # add group data to all data:
    
    if (is.data.frame(df)) {
      df <- rbind(df,g_df)
    } else {
      df <- g_df
    }
    
  }
  
  # # make factors:
  # for (column in c('participant','condition','group')) {
  #   df[,column] <- as.factor(df[,column])
  # }
  
  # do ANOVA on all groups:
  if (!is.null(exp)) {
    e_aov <- afex::aov_ez(id = "participant",
                          dv = "reachdeviation_deg",
                          data = df,
                          within = c("block","condition"),
                          between = "group",
                          type = 3 )
    
    cat(sprintf('\n===== Experiment: %d =====\n',exp ))
    print(summary(e_aov))
  } else {
    # print the group-based ANOVAs:
    for (group in names(group_anovas)) {
      cat(sprintf('\n===== %s =====\n',toupper(group) ))
      print(summary(group_anovas[[group]]))
    }
  }
  
}

conditionBayesFactors <- function(exp=NULL, groups=NULL) {
  
  if (is.null(groups)) {
    # we need to know the groups in the experiment:
    info <- getInfo()
    groups <- info$group[which(info$experiment == exp)]
  }
  
  # we get the data form these groups:
  data <- getSelectedGroupsData(groups=groups)  
  
  # we pre-process it, so that it is ready for BF stats
  data <- parseData(data, FUN=baseline)
  data <- parseData(data, FUN=addtrial)
  data <- parseData(data, FUN=normalize)
  data <- parseData(data, FUN=block)
  
  # we collect all data here for an ANOVA on all the data:
  df <- NA
  
  # ttestBF {BayesFactor}
  
  group_BFs <- list()
  
  for (group in names(data)) {
    
    # we collect group data here, so we can do group-based ANOVAs:
    g_df <- NA
    
    group_BFs[[group]] <- list()
    
    for (condition in c('abrupt', 'gradual')) {
      
      # add group and condition variables to the data frame:
      gc_df <- data[[group]][[condition]]
      gc_df$group <- group
      gc_df$condition <- condition
      
      # collect data in group data frame:
      if (is.data.frame(g_df)) {
        g_df <- rbind(g_df, gc_df)
      } else {
        g_df <- gc_df
      }
      
    }
    
    # do Bayes Factors for all group
    
    x = g_df$reachdeviation_deg[which(g_df$condition == 'abrupt' & g_df$block == 'rotated')]
    y = g_df$reachdeviation_deg[which(g_df$condition == 'gradual' & g_df$block == 'rotated')]
    group_BFs[[group]][['rotated']] <- BayesFactor::ttestBF(x=x, y=y, paired=TRUE)
    x = g_df$reachdeviation_deg[which(g_df$condition == 'abrupt' & g_df$block == 'clamped')]
    y = g_df$reachdeviation_deg[which(g_df$condition == 'gradual' & g_df$block == 'clamped')]
    group_BFs[[group]][['clamped']] <- BayesFactor::ttestBF(x=x, y=y, paired=TRUE)
    
    # add group data to all data:
    
    if (is.data.frame(df)) {
      df <- rbind(df,g_df)
    } else {
      df <- g_df
    }
    
  }
  
  if (!is.null(exp)) {
    cat(sprintf('\n===== Experiment: %d =====\n',exp ))
    x = df$reachdeviation_deg[which(df$condition == 'abrupt' & df$block == 'rotated')]
    y = df$reachdeviation_deg[which(df$condition == 'gradual' & df$block == 'rotated')]
    print(summary(BayesFactor::ttestBF(x=x, y=y, paired=TRUE)))
    x = df$reachdeviation_deg[which(df$condition == 'abrupt' & df$block == 'clamped')]
    y = df$reachdeviation_deg[which(df$condition == 'gradual' & df$block == 'clamped')]
    print(summary(BayesFactor::ttestBF(x=x, y=y, paired=TRUE)))
  
  } else {
    # print the group-based ANOVAs:
    for (group in names(group_BFs)) {
      cat(sprintf('\n===== %s =====\n',toupper(group) ))
      print(summary(group_BFs[[group]][['rotated']]))
      print(summary(group_BFs[[group]][['clamped']]))
    }
  }
  
  
}


# model quality and parameter comparisons -----

crossConditionMSEs <- function(exp=NULL, groups=NULL) {
  
  if (is.null(groups)) {
    # we need to know the groups in the experiment:
    info <- getInfo()
    groups <- info$group[which(info$experiment == exp)]
  }
  
  for (group in groups) {
    
    MSEs <- read.csv(sprintf('data/%s/modelMSEs.csv', group))
    
    cat(sprintf('\n--- %s: ---\n',toupper(group)))
    print( quantile(as.numeric(MSEs[,'abfits2grdata.MSE']) - as.numeric(MSEs[,'grfits2grdata.MSE']), probs=c(0.025, 0.975)) )
    
  }
  
}

conditionParameterDiffs <- function(exp=NULL, groups=NULL) {
  
  if (is.null(groups)) {
    # we need to know the groups in the experiment:
    info <- getInfo()
    groups <- info$group[which(info$experiment == exp)]
  }
  
  for (group in groups) {
    
    fits <- read.csv(sprintf('data/%s/bootstrapped_TwoRateFits.csv', group))

    cat(sprintf('\n--- %s ---\n',toupper(group)))
    
    for (parameter in c('Ls','Rs','Lf','Rf')) {
      column1 <- sprintf('abrupt.%s',parameter)
      column2 <- sprintf('gradual.%s',parameter)
      cat(sprintf('    %s:\n',parameter))
      print( quantile( as.numeric(fits[,column1]) - as.numeric(fits[,column2] ), 
                       probs=c(0.025, 0.975)) )
    }
    
    
    
  }
  
}