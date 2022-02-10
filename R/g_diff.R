
g.limma = function(g, deps=T){
  if(deps) g=g.run.deps(g, c("g.make.contrasts", "g.combine.met"))
  ##########################################################
  # Limma for DE
  # Use limma for differential analysis for each omics data set in omicsList.
  for(i in 1:length(g$omicsList)){ try({
    eset <- g$omicsList[[i]][["eSet"]]
    IsThereTime <- FALSE; 
    # Create a model matrix and lmFit based on group - For more complex models, edit the following code in this chunk.
    
    # This will try to fit a model for the TimeSeries case, where there are groups changing over time.
    if(g$time_index>0 ){ try({
      eset<- eset[,which(pData(eset)$Group %in% g$contrastgroups)]
      f <- factor(pData(eset)$Group, levels=unique(g$contrastgroups)); 
      Time <- splines::ns(as.numeric(pData(eset)$TimeSeries), df=(g$time_index))
      design <- stats::model.matrix(~f*Time);
      fit1 <- limma::lmFit( eset, design);
      fit1 <- limma::eBayes(fit1);
      IsThereTime <- TRUE
    }) }
      
    # If there is no TimeSeries (either not present or the model fails), we will next chech to see if we are doing a "Paired" analysis.
    # This is added with a "Pairs" entry on the second page of the annotation sheet to specify which samples are matched.
    if(IsThereTime == FALSE){
      if("Pairs" %in% colnames(pData(eset)) ){
        f <- factor(pData(eset)$Group)
        p <- factor(pData(eset)$Pairs)
        design <- stats::model.matrix(~p+f);
        fit1 <- limma::lmFit( eset, design);
        fit1 <- limma::eBayes(fit1);
        
        # We adjust these variables for the case where there are paired samples.
        g$loop_list <- (length(levels(p))+1): ncol(design)
        g$contrast_strings <- colnames(design)[g$loop_list]
        g$contrast_strings_file <- gsub("-","_",g$contrast_strings)
      } else {
        # This last case is for the model without TimeSeries or Pairs
        # This is the default case, for comparisons/contrasts between two or more groups.
        f <- factor(pData(eset)$Group)#, levels=unique(pData(eset)$Group)); 
        design <- stats::model.matrix(~ 0 + f);
        colnames(design) <- make.names(levels(f));
        fit1 <- limma::lmFit( eset, design);
        contrast_matrix <- limma::makeContrasts(contrasts=g$contrast_strings, levels=design);
      
        fit1 <- limma::contrasts.fit(fit1, contrast_matrix);
        fit1 <- limma::eBayes(fit1);
      }
    }
    # Summary of the fit
    #end_col <- ncol(topTable(fit1, adjust="BH"))
    #start_col <- end_col - (3+length(contrast_list))
    #print(topTable(fit1, adjust="BH")[,seq(start_col, end_col)])
    
    # Save the fit in the omicsList object for the corresponding data set
    g$omicsList[[i]][[7]] <- fit1
    names(g$omicsList[[i]])[7] <- "fit";
  }) } 

  g$calls = c(g$calls, "g.limma")
  g
}

g.fc = function(g, deps=T){
  if(deps) g=g.run.deps(g, c("g.make.contrasts", "g.combine.met"))
  ##########################################################
  # Calculate log fold change in case limma fails
  # For the case where the differential analysis fails (e.g., due to too few samples per group), we calculate an average value for each group,
  # and then calculate a simple difference between all groups.
  logfc_index <- c();
  for( i in 1:length(g$omicsList) ){ try({
    # This first loop calculates a mean value per group.
    for( j in 1:length(g$contrastgroups) ){
      if( sum(pData(g$omicsList[[i]][["eSet"]])$Group==g$contrastgroups[j])>1 ){
        fData(g$omicsList[[i]][["eSet"]])[,paste("mean_",g$contrastgroups[j],sep="")] <-
          rowMeans(exprs(g$omicsList[[i]][["eSet"]])[,pData(g$omicsList[[i]][["eSet"]])$Group==g$contrastgroups[j] ] )
      } else {
        fData(g$omicsList[[i]][["eSet"]])[,paste("mean_",g$contrastgroups[j],sep="")] <-
          exprs(g$omicsList[[i]][["eSet"]])[,pData(g$omicsList[[i]][["eSet"]])$Group==g$contrastgroups[j] ] 
      }
    }
    # This loop creates a log fold change between each two groups. (Log fold change for log transformed values is approximately the difference.)
    for( j in 1:(length(g$contrastgroups)-1) ){
      for( k in 2:(length(g$contrastgroups)) ){
        col_name <- paste("logfc_",g$contrastgroups[k],"_",g$contrastgroups[j],sep="");
        fData(g$omicsList[[i]][["eSet"]])[,col_name] <- (fData(g$omicsList[[i]][["eSet"]])[,paste("mean_",g$contrastgroups[k],sep="")] -
                                                       fData(g$omicsList[[i]][["eSet"]])[,paste("mean_",g$contrastgroups[j],sep="")]    ) 
        logfc_index <- c(logfc_index, col_name)
      }
    }
    # We calculate a maximum fold change value across all groups as a way to see which features are changing most across all groups.
    # This is analogous to the F-statistic from the differential analysis.
    fData(g$omicsList[[i]][["eSet"]])[,"logfc_Overall"] <- apply(fData(g$omicsList[[i]][["eSet"]])[,grep("mean", colnames(fData(g$omicsList[[i]][["eSet"]])))], 1, function(x) max(x) - min(x))
  }) }
  g$calls = c(g$calls, "g.fc")
  g$logfc_index <- unique(logfc_index);
  g
}

g.volcano = function(g, deps=T){
  if(deps) g=g.run.deps(g, "g.limma")
  output_links = ""

  for(i in 1:length(g$omicsList)){
    if(class(g$omicsList[[i]][["fit"]])=='MArrayLM'){
      subset_rows=FALSE;
      if( class(g$subset_genes)!="logical" ){ try({
        subset_rows <- vector("list", length(g$subset_genes))
        names(subset_rows) <- names(g$subset_genes)
        for( j in 1:length(g$subset_genes) ){
          if( "Gene" %in% colnames(fData(g$omicsList[[i]][["eSet"]])) ){
            subset_rows[[j]] <- rownames(g$omicsList[[i]][["eSet"]])[which(fData(g$omicsList[[i]][["eSet"]])$Gene %in% g$subset_genes[[j]] )]
          }
          if( length(subset_rows[[j]])==0 ){
            subset_rows[[j]] <- rownames(g$omicsList[[i]][["eSet"]])[ which(rownames(g$omicsList[[i]][["eSet"]]) %in% g$subset_genes[[j]] ) ]
          }
        }
      }, silent=TRUE) }
      
      for(c in 1:length(g$loop_list)){
        # For the data set and coefficient, make the summary table and name
        top_sum <- limma::topTable(g$omicsList[[i]][["fit"]], adjust="BH", n=Inf, sort.by='p', coef=g$loop_list[c]);
        type_name <- paste(g$omicsList[[i]][["dataType"]],g$contrast_strings[c], sep="_");
      
        # run Volcan function
        drawVolcano(dat=top_sum, type=type_name, subset_rows=subset_rows, top_values=adjpcutoff,top_fc=g$fc_cutoff,
                    outputpath=g$output_contrast_path);
      
        # Make output hyperlinks for markdown
        add_link <- paste("[ ",type_name, " ](", g$output_contrast_subdir,"/",type_name, "_volcano",".pdf)", sep="")
        output_links <- paste(output_links, add_link, sep=" | " )   
      }
    }
    output_links <- paste(output_links, "  \n", sep="" );
  }

  g$calls = c(g$calls, "g.volcano")
  g$volcano_output_links = output_links
  g
}

g.md = function(g, deps=T){
  if(deps) g=g.run.deps(g, "g.fc")
  output_links = g$volcano_output_links

  for(i in 1:length(g$omicsList)){
      subset_rows=FALSE;
      if( class(g$subset_genes)!="logical" ){ try({
        subset_rows <- vector("list", length(g$subset_genes))
        names(subset_rows) <- names(g$subset_genes)
        for( j in 1:length(g$subset_genes) ){
          if( "Gene" %in% colnames(fData(g$omicsList[[i]][["eSet"]])) ){
            subset_rows[[j]] <- rownames(g$omicsList[[i]][["eSet"]])[which(fData(g$omicsList[[i]][["eSet"]])$Gene %in% g$subset_genes[[j]] )]
          }
          if( length(subset_rows[[j]])==0 ){
            subset_rows[[j]] <- rownames(g$omicsList[[i]][["eSet"]])[ grep( paste(g$subset_genes[[j]], collapse="|"), 
                                                                        rownames(exprs(g$omicsList[[i]][["eSet"]])) ) ]
          }
        }
    }, silent=TRUE) }
    
    for(c in 1:length(g$logfc_index) ){ try({
      # For the data set and coefficient, make the summary table and name
      top_sum <- data.frame(logFC= fData(g$omicsList[[i]][["eSet"]])[,g$logfc_index[c]],
                            Mean=rowMeans(exprs(g$omicsList[[i]][["eSet"]])), 
                            feature_identifier =fData(g$omicsList[[i]][["eSet"]])[,"feature_identifier"] );
      if( "Gene" %in% colnames(fData(g$omicsList[[i]][["eSet"]])) ){ top_sum$Gene <- fData(g$omicsList[[i]][["eSet"]])$Gene }
      else if( "mz" %in% colnames(fData(g$omicsList[[i]][["eSet"]])) ){ top_sum$mz <- fData(g$omicsList[[i]][["eSet"]])$mz }
      
      type_name <- paste(g$omicsList[[i]][["dataType"]],g$logfc_index[c], sep="_");
      
      # run Volcano function
      if(g$fc_cutoff==0){ tmp <- 1 } else { tmp <- g$fc_cutoff }
      drawMDPlot(dat=top_sum, type=type_name, subset_rows=subset_rows, cutoff=tmp,
                 outputpath=g$output_contrast_path);
      
      # Make output hyperlinks for markdown
      add_link <- paste("[ ",type_name, " ](", g$output_contrast_subdir,"/",type_name,"_MDplot",".pdf)", sep="")
      output_links <- paste(output_links, add_link, sep=" | " )   
    }) }
    output_links <- paste(output_links, "  \n", sep="" );
  }

  g$calls = c(g$calls, "g.md")
  g$volcano_output_links = output_links
  g
}

g.de.static.heatmap = function(g, deps=T){
  if(deps) g=g.run.deps(g, c("g.fc", "g.limma"))
  output_links = ""
  for(i in 1:length(g$omicsList)){
    sig_cutoff <- round(sig_percent*nrow(g$omicsList[[i]][["eSet"]]), digits=0)
    if(class(g$omicsList[[i]][["fit"]])=='MArrayLM'){
    for(c in 1:length(g$loop_list)){ try({
      # For the data set and coefficient, make the summary table and name
      top_sum <- rownames( limma::topTable(g$omicsList[[i]][["fit"]], adjust="BH", n=Inf, sort.by='p', coef=g$loop_list[c], p.value=adjpcutoff) );
      title_add <- paste("FDR<", adjpcutoff, sep="")
      if (length(top_sum)< sig_cutoff ){ 
        top_sum <- rownames( limma::topTable(g$omicsList[[i]][["fit"]], adjust="BH", n=Inf, sort.by='p', coef=g$loop_list[c]))[1:sig_cutoff];
        title_add <- paste("Top ", (sig_percent*100), "%", sep="")
      }
      type_name <- paste(g$omicsList[[i]][["dataType"]],g$contrast_strings[c], sep="_");
      
      # run heatmap function
      drawHeatmaps(eset=g$omicsList[[i]][["eSet"]], limmaSig=top_sum, type=type_name, title_add=title_add,
                    outputpath=g$output_plots_path, outputcontrastpath=g$output_contrast_path);
      
      if(g$num_contrasts>1) {
          drawHeatmaps(eset=g$omicsList[[i]][["eSet"]][,pData(g$omicsList[[i]][["eSet"]])$Group %in% unlist(strsplit(g$contrast_strings[c], "-")) ], outputpath=g$output_plots_path, outputcontrastpath=g$output_contrast_path,
                      limmaSig=top_sum,type=paste(type_name, "_subgroup", sep="") );
      }
      
      # Make output hyperlinks for markdown
      add_link <- paste("[ ",type_name, " ](", g$output_contrast_subdir,"/",type_name,"_heatmaps",".pdf)", sep="")
      if(g$num_contrasts>1){ 
        add_link <- paste(add_link, " | [ ",type_name, "_subgroup ](", g$output_contrast_subdir,"/",type_name,"_subgroup_heatmaps.pdf) \n", sep=""); }
      output_links <- paste(output_links, add_link, sep=" | " ); 
      
    }) } } else {
    for(c in 1:length(g$logfc_index) ){ try({  
      top_sum <- rownames( exprs(g$omicsList[[i]][["eSet"]])[order( abs(fData(g$omicsList[[i]][["eSet"]])[,g$logfc_index[c]]), decreasing=TRUE),] )[1:sig_cutoff]
      type_name <- paste(g$omicsList[[i]][["dataType"]],g$contrast_strings[c], "logFC", sep="_");
      title_add <- paste("Top ", (sig_percent*100), "%", sep="")
      
      # run heatmap function
      drawHeatmaps(eset=g$omicsList[[i]][["eSet"]], limmaSig=top_sum, type=type_name, title_add=title_add, outputpath=g$output_plots_path, outputcontrastpath=g$output_contrast_path);
      
      # Make output hyperlinks for markdown
      add_link <- paste("[ ",type_name, " ](", g$output_contrast_subdir,"/",type_name,"_heatmaps",".pdf)", sep="")
      output_links <- paste(output_links, add_link, sep=" | " ); 
    }) }
    }
    output_links <- paste(output_links, "  \n", sep="" );
  }

  # Heatmaps for logfc_Overall and F statistic
  if(g$statistic_index=="F"){
  for(i in 1:length(g$omicsList)){
    sig_cutoff <- round(sig_percent*nrow(g$omicsList[[i]][["eSet"]]), digits=0)
    if(class(g$omicsList[[i]][["fit"]])=='MArrayLM' ){ try({
      top_sum <- rownames( limma::topTableF(g$omicsList[[i]][["fit"]], adjust="BH", n=Inf, sort.by='F', p.value=adjpcutoff) );
      title_add <- paste("FDR<", adjpcutoff, sep="")
      if (length(top_sum)<sig_cutoff){ 
        top_sum <- rownames( limma::topTableF(g$omicsList[[i]][["fit"]], adjust="BH", n=Inf, sort.by='F' ) )[1:sig_cutoff];
        title_add <- paste("Top ", (sig_percent*100), "%", sep="")
      }
      type_name <- paste(g$omicsList[[i]][["dataType"]], "_","F-statistic", sep="");
      
      # run heatmap function
      drawHeatmaps(eset=g$omicsList[[i]][["eSet"]], limmaSig=top_sum, type=type_name, title_add=title_add, outputpath=g$output_plots_path, outputcontrastpath=g$output_contrast_path);
      
      # Make output hyperlinks for markdown
      add_link <- paste("[ ",type_name, " ](", g$output_contrast_subdir,"/",type_name,"_heatmaps",".pdf)", sep="")
      output_links <- paste(output_links, add_link, sep=" | " ); 
    }) }
    if(class(g$omicsList[[i]][["fit"]])=='MArrayLM' & g$time_index>0){ try({
      top_sum <- rownames( limma::topTable(g$omicsList[[i]][["fit"]], adjust="BH", n=Inf,coef=g$time_start:g$time_end, p.value=adjpcutoff) );
      title_add <- paste("FDR<", adjpcutoff, sep="")
      if (length(top_sum)<sig_cutoff){
        top_sum <- rownames( limma::topTable(g$omicsList[[i]][["fit"]], adjust="BH", n=Inf, coef=g$time_start:g$time_end ) )[1:sig_cutoff]; 
        title_add <- paste("Top ", sig_percent, "%", sep="")
      }
      type_name <- paste(g$omicsList[[i]][["dataType"]], "_", "TimeCourse_Overall", sep="");
      
      # run heatmap function
      drawHeatmaps(eset=g$omicsList[[i]][["eSet"]], limmaSig=top_sum, type=type_name, title_add=title_add, outputpath=g$output_plots_path, outputcontrastpath=g$output_contrast_path);
      
      # Make output hyperlinks for markdown
      add_link <- paste("[ ",type_name, " ](", g$output_contrast_subdir,"/",type_name,"_heatmaps",".pdf)", sep="")
      output_links <- paste(output_links, add_link, sep=" | " ); 
    }) }
    try({  
      top_sum <- rownames( exprs(g$omicsList[[i]][["eSet"]])[order( abs(fData(g$omicsList[[i]][["eSet"]])[,"logfc_Overall"]), decreasing=TRUE),] )[1:sig_cutoff]
      type_name <- paste(g$omicsList[[i]][["dataType"]],"logfc_Overall", sep="_");
      title_add <- paste("Top ", (sig_percent*100), "%", sep="")
      
      # run heatmap function
      drawHeatmaps(eset=g$omicsList[[i]][["eSet"]], limmaSig=top_sum, type=type_name, title_add=title_add, outputpath=g$output_plots_path, outputcontrastpath=g$output_contrast_path);
      
      # Make output hyperlinks for markdown
      add_link <- paste("[ ",type_name, " ](", g$output_contrast_subdir,"/",type_name,"_heatmaps",".pdf)", sep="")
      output_links <- paste(output_links, add_link, sep=" | " ); 
    }) 
    output_links <- paste(output_links, "  \n", sep="" );
  }
  }

  g$calls = c(g$calls, "g.de.static.heatmap")
  g$heatmap_output_links = output_links
  g
}

g.interactive.volcano = function(g, deps=T){
  if(deps) g=g.run.deps(g, c("g.fc", "g.limma"))
  output_links<-"";
  for(i in 1:length(g$omicsList)){
    if(class(g$omicsList[[i]][["fit"]])=='MArrayLM' & g$time_index==0){
    for(c in 1:length(g$loop_list) ){ try({
      top_sum <- limma::topTable(g$omicsList[[i]][["fit"]], adjust="BH", n=Inf, sort.by='p', coef=c);
      type_name <- paste(g$omicsList[[i]][["dataType"]],g$contrast_strings[c], sep="_");

      interactiveVolcano(eset=g$omicsList[[i]][["eSet"]], fit=g$omicsList[[i]][["fit"]], dt=decideTests(g$omicsList[[i]][["fit"]]),
                        limmaSig=top_sum, type=type_name, col=g$loop_list[c], outputcontrastpath=g$output_contrast_path);

      add_link <- paste("[ ",type_name, " ](", g$output_contrast_subdir,"/InteractivePlots/Volcano-Plot_",type_name,".html)", sep="");
      output_links <- paste(output_links, add_link, sep=" | " ); 
    }) }
    } else {
    for(c in 1:length(g$logfc_index) ){ try({
      type_name <- paste(g$omicsList[[i]][["dataType"]],g$logfc_index[c], sep="_");
      top_sum <- data.frame(logFC= fData(g$omicsList[[i]][["eSet"]])[,g$logfc_index[c]],
                            Mean=rowMeans(exprs(g$omicsList[[i]][["eSet"]])), 
                            feature_identifier =fData(g$omicsList[[i]][["eSet"]])[,"feature_identifier"] );
      if( "Gene" %in% colnames(fData(g$omicsList[[i]][["eSet"]])) ){ top_sum$Gene <- fData(g$omicsList[[i]][["eSet"]])$Gene }
      if( "mz" %in% colnames(fData(g$omicsList[[i]][["eSet"]])) ){ top_sum$mz <- fData(g$omicsList[[i]][["eSet"]])$mz }
      
      interactiveVolcano(eset=g$omicsList[[i]][["eSet"]],limmaSig=top_sum, type=type_name, col=c, outputcontrastpath=g$output_contrast_path);
      add_link <- paste("[ ",type_name, " ](", g$output_contrast_subdir,"/InteractivePlots/Volcano-Plot_",type_name,".html)", sep="");
      output_links <- paste(output_links, add_link, sep=" | " ); 
    }) }
    }
    output_links <- paste(output_links, "  \n", sep="" ); 
  }

  g$calls = c(g$calls, "g.interactive.volcano")
  g$intvolc_output_links = output_links
  g
}

ttable=function(ol,coef){
  limma::topTable(ol, adjust="BH", n=Inf, sort.by='p', coef=coef)
}

fttable=function(ol,coef){
  top_sum <- limma::topTableF(ol, adjust="BH", n=Inf, sort.by='F'); 
  top_sum[,"logFC"] <- top_sum[,"F"]
  top_sum
}

timetable=function(ol,coef){
  top_sum <- limma::topTable(ol, adjust="BH", n=Inf, sort.by='F', coef=coef); 
  top_sum[,"logFC"] <- top_sum[,"F"]
  top_sum
}

genfile=function(g,oListi,contrast_name,coef,tsum,fd){
  top_sum <- FALSE;

  try({ top_sum <- tsum(oListi[["fit"]], coef); })
  if( class(top_sum)=="logical" ){ try({
    top_sum <- data.frame(logFC=fd(fData(oListi[["eSet"]]),contrast_name),
                        Mean=rowMeans(exprs(oListi[["eSet"]])), 
                        feature_identifier=fData(oListi[["eSet"]])[,"feature_identifier"] );
    if( "Gene" %in% colnames(fData(oListi[["eSet"]])) ){ top_sum$Gene <- fData(oListi[["eSet"]])$Gene }
    if( "mz" %in% colnames(fData(oListi[["eSet"]])) ){ top_sum$mz <- fData(oListi[["eSet"]])$mz }
  })}

  saveFiles(data=oListi, limmaRes=top_sum, type=oListi[["dataType"]], contrast_name=contrast_name, outputpath=g$output_plots_path, outputcontrastpath=g$output_contrast_path_files)
  
  if( !( grepl("Human", species) ) & "Gene" %in% colnames(fData(oListi[["eSet"]]))   ){ 
    top_sum$Gene <- toupper(top_sum$Gene)
    saveFiles(oListi, limmaRes=top_sum, type=paste(oListi[["dataType"]], "_Uppercase", sep=""), contrast_name=contrast_name, saveRDS=F, outputpath=g$output_plots_path, outputcontrastpath=g$output_contrast_path_files)
  }
}

g.savedata = function(g, deps=T){
  if(deps) g=g.run.deps(g, c("g.fc", "g.limma"))

  iters=do.call("rbind",
    lapply(1:length(g$omicsList),function(oi){
      t(sapply(1:length(g$loop_list),function(li){
        c(oi,li)
      }))
    })
  )
  colnames(iters)=c("omics","loop")

  BiocParallel::bplapply(1:nrow(iters),FUN=function(ii){ try({
    i=iters[ii,"omics"]
    j=iters[ii,"loop"]
    genfile(g,g$omicsList[[i]],g$contrast_strings[j],g$loop_list[j],ttable,function(ol,cn){ol[,paste("logfc_",gsub("-","_",cn),sep="")]});
  })})

  if(g$statistic_index=='F' || g$time_index>0){
    BiocParallel::bplapply(1:length(g$omicsList),FUN=function(i){ try({
    if(g$statistic_index=='F'){
      genfile(g,g$omicsList[[i]],"Overall",NULL,fttable,function(ol,cn){ol[,"logfc_Overall"]});
    }
    if(g$time_index>0){
      genfile(g,g$omicsList[[i]],"Timecourse_Overall",g$time_start:g$time_end,timetable,NULL);
    }
  })})}

  g$calls = c(g$calls, "g.savedata")
  g
}
