# Balance Plot Generation
# make a funciton that creates balance plots (not in strata)
# arguments include: 
# data = data set
# weights = column of weights 
# treatment = treatment variable
# list of covs

bal.plots <- function(df, weights, treatment, covs){
 
    # Variance 
    data.var <- df %>% group_by(.data[[treatment]]) %>% 
      summarize(across(covs, ~var(.x))) %>% as.data.frame() # if trouble with across, restart R session
    
    # Calculate pooled var in the overall unweighted data
    c.var <- as.numeric(data.var[1,])
    t.var <- as.numeric(data.var[2,])
    c.var <- c.var[-1]
    t.var <- t.var[-1]
    pooled.var <- sqrt((t.var + c.var)/2)
    
    # Calculate the mean in the unweighted data
    um.wt <- df %>% group_by(.data[[treatment]]) %>% 
      summarize(across(covs, ~mean(.x))) %>% as.data.frame()
    
    # Calculate the mean in the weighted data
    bal.st <- df %>% group_by(.data[[treatment]]) %>% 
      summarize(across(covs, ~ weighted.mean(.x, .data[[weights]]))) %>% as.data.frame()
    
    # Make table of unweighted means in treated, untreated, and SMD
    um.wt.tab <- matrix(NA, length(covs), 3)
    um.wt.tab[,1] <- unlist(um.wt[1,-1]) 
    um.wt.tab[,2] <- unlist(um.wt[2,-1])                        
    um.wt.tab[,3] <- (unlist(um.wt[2,-1]) - unlist(um.wt[1,-1]))/pooled.var
    
    # Make table of weighted means in treated, untreated, and SMD
    bal.wt.tab <- matrix(NA, length(covs), 3)
    bal.wt.tab[,1] <- unlist(bal.st[1,-1]) 
    bal.wt.tab[,2] <- unlist(bal.st[2,-1])                        
    bal.wt.tab[,3] <- (unlist(bal.st[2,-1]) - unlist(bal.st[1,-1]))/pooled.var 
    
    ## Rename columns using covs
    rownames(um.wt.tab) <- covs # <--- here give better names? 
    rownames(bal.wt.tab) <- covs

    n.covs <- length(covs)
    
    ## Total Imbalance Reduction
    um.wt.bias <- um.wt.tab[,3]
    bal.bias <- bal.wt.tab[,3] 
    pbr <- (1 - (mean(abs(bal.bias))/mean(abs(um.wt.bias))))*100
    message("Total Imbalance Reduction is: " , pbr)
    
    # Plots              			  
    data.plot <- c( bal.wt.tab[,3], um.wt.tab[,3])
    data.plot <- as.data.frame(data.plot)
    names(data.plot) <- "std.dif"
    data.plot$contrast <- c(rep(1, n.covs), rep(2, n.covs))
    data.plot$contrast <- factor(data.plot$contrast, levels = c(1,2), 
                                 labels = c("Weighted", "Unweighted"))
    data.plot$covariate <- as.factor(rownames( bal.wt.tab))
    
    
    plot <- ggplot(data = data.plot, aes(x = std.dif, y = covariate, 
                                         shape = factor(contrast), color = factor(contrast))) +
      geom_point(size = 2, shape = 16) +
      scale_shape_manual(name = "Contrast", values = c(1, 15)) +
      scale_color_manual(name = "Contrast", values = c("darkred", "darkblue")) +
      xlab("Standardized Difference") + ylab("Covariates") +
      ggtitle("Balance Plot") +
      scale_y_discrete(limits = rev(levels(data.plot$covariate))) +
      geom_vline(xintercept = 0) +
      geom_vline(xintercept = 0.1) +
      geom_vline(xintercept = -0.1) +
      theme_bw() +
      guides(shape = FALSE)  # Hide the shape legend
    
    print(plot)
    print(1)#


  return(list(plot = plot, data = data.plot))
}



