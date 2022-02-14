g.norm = function(g, deps=T){
  if(deps) g=g.run.deps(g, "g.make.omicsList")

  for(i in 1:length(g$omicsList)){
    if(length(zero_percent)==1){use_zero<-zero_percent}else{use_zero<-zero_percent[i]}
    if(length(norm_method)==1){use_norm<-norm_method}else{use_norm<-norm_method[i]}

    #if(g$debug_opt){ index<-14 } else { index <- 5 }
    index=5
    g$omicsList[[i]][[5]] <-intensityNorm(eset=g$omicsList[[i]][[index]],
                                        type=g$omicsList[[i]][[1]], 
                                        data_format=g$omicsList[[i]][["dataFormat"]],
                                        norm=norm_method, zero_cutoff=zero_percent,
                                        min_feature=g$min_feature_per_sample,
                                        annotate=g$annot, outputpath=g$output_plots_path,
                                        norm_by_batch=g$norm_batches)

    # Batch correction
    if("Batch" %in% colnames(pData(g$omicsList[[i]][["eSet"]])) & !("Met" %in% g$omicsList[[i]][["dataFormat"]])) { try({
      suppressPackageStartupMessages({require(sva) });
      suppressMessages({ suppressWarnings({
        g$omicsList[[i]][[12]] <- g$omicsList[[i]][["eSet"]]
        names(g$omicsList[[i]])[12] <- "prebatch_eset"

        pheno = pData(g$omicsList[[i]][["eSet"]])
        edata = exprs(g$omicsList[[i]][["eSet"]])
        batch = pheno$Batch

        modcombat = model.matrix(~1, data=pheno)
        combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
        combat_eset = ExpressionSet(assayData=combat_edata)
        pData(combat_eset) = pData(g$omicsList[[i]][["eSet"]])
        fData(combat_eset) = fData(g$omicsList[[i]][["eSet"]])
      })})
      g$omicsList[[i]][["eSet"]] <- combat_eset
      tmp<-drawPCA(eset=g$omicsList[[i]][["prebatch_eset"]], type=paste(g$omicsList[[i]][["dataType"]],"_precorrection", sep=""),
                   outputpath=g$output_plots_path, outputfile=g$output_files_path);

      capture.output( intensityNorm(eset=g$omicsList[[i]][[index]],
                                    type=paste(g$omicsList[[i]][[1]],"_postBatchCor",sep=""), 
                                    data_format=g$omicsList[[i]][["dataFormat"]],
                                    norm=g$norm_post_batch, zero_cutoff=zero_percent,
                                    annotate=g$annot, outputpath=g$output_plots_path,
                                    min_feature=g$min_feature_per_sample, norm_by_batch=F) )

    }) }
  }

  g$calls = c(g$calls, "g.norm")
  g
}

g.norm.vis = function(g, deps=T){
  if(deps) g=g.run.deps(g, "g.norm")
  pca_output <- vector("list", length(g$omicsList) )
  output_links<-"";

  for(i in 1:length(g$omicsList)){

    # Make output hyperlinks for markdown
    add_link <- paste(g$omicsList[[i]][["dataType"]], ": [ Boxplot/Dist. ](", g$output_plots_subdir,
                      "/boxplot_histogram_",g$omicsList[[i]][["dataType"]],".pdf)", sep="")
    output_links <- paste(output_links, add_link, sep="" )


    if(g$remove_group!=""){
      g$omicsList[[i]][[5]] <- g$omicsList[[i]][[5]][, !grepl(paste(remove_group, collapse="|"),pData(g$omicsList[[i]][[5]])$Group) ]
    }

    try({
      # PCA plots
      pca_output[[i]] <- drawPCA(eset=g$omicsList[[i]][["eSet"]], type=g$omicsList[[i]][["dataType"]], outputpath=g$output_plots_path, outputfile=g$output_files_path)

      # Make output hyperlinks for markdown
      add_link <- paste( "[ PCA ](", g$output_plots_subdir,"/PCAplots_",g$omicsList[[i]][["dataType"]],".pdf)", sep="")
      output_links <- paste(output_links, add_link, sep=" | " )
    })
    try({ 
      drawUMAP(eset=g$omicsList[[i]][["eSet"]], type=g$omicsList[[i]][["dataType"]], outputpath=g$output_plots_path)
      add_link <- paste( "[ UMAP ](", g$output_plots_subdir,"/UMAPplot_",g$omicsList[[i]][["dataType"]],".pdf)", sep="")
      output_links <- paste(output_links, add_link, sep=" | " )
    })

    output_links <- paste(output_links, "  \n", sep="" )
  }

  g$calls = c(g$calls, "g.norm.vis")
  g$pca_output = pca_output
  g$output_links = output_links
  g
}

g.venn = function(g, deps=T){
  if(deps) g=g.run.deps(g, "g.make.omicsList")
  output_links <-"";

  draw_vdata=function(data_index,var,name){
    v_list <- vector("list", length(data_index))
    for(i in 1:length(data_index)){
      v_list[[i]] <- fData(g$omicsList[[ data_index[i] ]][["eSet"]])[,var]
      if(class(v_list[[i]])=="numeric"){
        v_list[[i]] <- stats::na.omit(round(v_list[[i]],digits=2))
      }
      names(v_list)[i] <- g$omicsList[[ data_index[i] ]][["dataType"]]
    }
    drawVenn(item_list=v_list, item_name=name, outputpath=g$output_plots_path)
    paste0("[ ",name," ](", g$output_plots_subdir,"/VennDiagram_",name,".pdf)")
  }
  output_links <- paste(c(
                          if( length(g$gene_data_index) > 1 ) draw_vdata(g$gene_data_index,"Gene","Genes"),
                          if( length(g$prot_data_index) > 1 ) draw_vdata(g$prot_data_index,"Protein","Proteins"),
                          if( length(g$mz_data_index) > 1 ) draw_vdata(g$mz_data_index,"mz","Metabolites_mz")
                          ), collapse=" | ")

  g$calls = c(g$calls, "g.venn")
  g$output_links = output_links
  g
}

drawcorplot <- function(g, item_name,data_index,dataname=NULL){ try({
  output_links <-"";
  if(is.null(dataname)) dataname <- rep("eSet",length(data_index))
  fprefix=paste0(dataname[length(data_index)],"_")
  outname <- paste0(sub("eSet_","",fprefix),g$contrastgroups)

  for (j in 1:length(g$contrastgroups) ){
    group_name <- paste0(item_name,"_", outname[j])
    gene_values <- vector("list", (length(data_index)) )
    for (i in 1:(length(data_index)) ){
      group_columns <- grepl(g$contrastgroups[j], pData(g$omicsList[[ data_index[i] ]][[ dataname[i] ]])$Group)
      gene_values[[i]] <- data.frame( exprs(g$omicsList[[ data_index[i] ]][[ dataname[i] ]])[,group_columns] ); # Make a data frame of values and genes
      gene_values[[i]][,item_name] <- fData(g$omicsList[[ data_index[i] ]][[ dataname[i] ]])[,item_name] ;                                                                         
      gene_values[[i]] <- tidyr::gather(gene_values[[i]], sample, value, -tidyselect::all_of(item_name)); # Gather values and get average by gene
      gene_values[[i]] <- dplyr::group_by_at(gene_values[[i]], tidyselect::all_of(item_name)) ;
      gene_values[[i]] <- dplyr::summarize(gene_values[[i]], mean=mean(value));
      gene_values[[i]] <- gene_values[[i]][gene_values[[i]][,item_name]!="",]
      names(gene_values)[i] <- paste0(g$omicsList[[ data_index[i] ]][["dataType"]],"_",outname[j])
    }
    drawXYCorr(item_list=gene_values, item_name=item_name, file_name = outname[j],
               subset_genes=g$subset_genes, outputpath=g$output_plots_path);
    add_link <- paste0("[ ",group_name," Correlation Plots ](",g$output_plots_subdir,"/Correlation_Plots_",group_name,".pdf)")
    output_links <- paste(output_links, add_link, sep=" | " )
  }
  output_links
}) }


g.corplot = function(g, deps=T){
  if(deps) g=g.run.deps(g, "g.norm")
  g$output_links <- paste(c(
                          if(length(g$gene_data_index)>1) drawcorplot(g,item_name="Gene",data_index=g$gene_data_index),
                          if(length(g$prot_data_index)>1) drawcorplot(g,item_name="Protein",data_index=g$prot_data_index)
                          ), collapse=" | ")
  g$calls = c(g$calls, "g.corplot")
  g
}

g.cor.across.groups = function(g, deps=T){
  if(deps) g=g.run.deps(g, "g.norm")

  output_links = BiocParallel::bplapply(1:length(g$omicsList),FUN=function(i){try({
    if(length(unique(pData(g$omicsList[[i]][["eSet"]])$Group))>1){
      gene_values <- data.frame(t(exprs(g$omicsList[[i]][["eSet"]])))
      gene_values$Group <- pData(g$omicsList[[i]][["eSet"]])$Group
      gene_values <- aggregate(. ~ Group, gene_values, mean)
      rownames(gene_values)<- gene_values$Group
      gene_values$Group <- NULL
      gene_values <- data.frame(t(gene_values))
      gene_values$Groups<- row.names(gene_values)

      item_list <- vector("list", ncol(gene_values)-1)
      for( j in 1:(ncol(gene_values)-1) ){
        item_list[[j]] <- tibble::as_tibble(cbind(gene_values[,"Groups"], gene_values[,j] ), .name_repair="minimal")
        names(item_list[[j]]) <-  c("Groups", colnames(gene_values)[j])
        item_list[[j]][,2] <- as.double(unlist(item_list[[j]][,2]))
      }
      names(item_list)<- colnames(gene_values)[1:(ncol(gene_values)-1)]
      file_name <- g$omicsList[[i]][["dataType"]]
      item_name <- "Groups"
      drawXYCorr(item_list=item_list, item_name=item_name, file_name=file_name, outputpath=g$output_plots_path);
      paste("[ ",file_name," ](",g$output_plots_subdir,"/Correlation_Plots_",item_name,"_",file_name,".pdf)", sep="");
    }
  }  ) })

  g$calls = c(g$calls, "g.cor.across.groups")
  g$output_links = paste0(output_links, collapse=" | " )
  g
}

g.norm.to.first = function(g, deps=T){
  if(deps) g=g.run.deps(g, "g.norm")
  output_links<-"";
  try({

    # Take first data set for normalization
    if(length(g$prot_data_index)>1){
      dataname=c("eSet",rep("siteNorm",length(g$prot_data_index)-1))
      item_name="Protein"
      data_index=g$prot_data_index
    } else {
      dataname=c("eSet",rep("siteNorm",length(g$gene_data_index)-1))
      item_name="Gene"
      data_index=g$gene_data_index
    }

    gene_values <- g$omicsList[[ data_index[1] ]][["eSet"]]; 
    for(i in 2:length(data_index)){
      if( identical(gene_values$SampleName, g$omicsList[[ data_index[i] ]][["eSet"]]$SampleName) ) {

        norm_values <- g$omicsList[[ data_index[i] ]][["eSet"]] ;
        norm_values <- norm_values[ fData(norm_values)[,item_name] %in% fData(gene_values)[,item_name] ,];

        nv = fData(norm_values)[,item_name]
        gv = fData(gene_values)[,item_name]
        gi = gv %in% nv
        gv = gv[gi]
        ag = aggregate(exprs(gene_values)[gi,],by=list(gv),FUN="mean")
        nm = merge(data.frame(Group.1=nv,ind=1:length(nv)),ag,by="Group.1")
        nm = nm[order(nm$ind),-2]
        rm = rowMeans(nm[,-1])
        temp_norm_values = as.matrix((exprs(norm_values) - nm[,-1]) + rm)
        rownames(temp_norm_values) = rownames(exprs(norm_values))
        exprs(norm_values) = temp_norm_values

        g$omicsList[[ data_index[i] ]][[10]] <- norm_values;
        names( g$omicsList[[ data_index[i] ]] )[10] <- "siteNorm";

        g$data_norm_index <- c(g$data_norm_index, data_index[i])
      }
    }

    output_links = drawcorplot(g,item_name=item_name,data_index=data_index,dataname=dataname)

    for(i in 1:length(g$data_norm_index)){
      type_name <- paste0(g$omicsList[[ g$data_norm_index[i] ]][["dataType"]], "_NormTo",g$omicsList[[ data_index[1] ]][["dataType"]])
      tmp<-drawPCA(eset=g$omicsList[[ g$data_norm_index[i] ]][["siteNorm"]], type=type_name, show_sample_names=TRUE, outputpath=g$output_plots_path, outputfile=g$output_files_path)

      add_link <- paste0("[ PCA:",type_name, " ](", g$output_plots_subdir,"/PCAplots_",type_name,".pdf)" )
      output_links <- paste(output_links, add_link, sep=" | " )
    }

  })

  g$calls = c(g$calls, "g.norm.to.first")
  g$output_links = output_links
  g
}

g.use.site.norm = function(g, deps=T){
  if(deps) g=g.run.deps(g, "g.make.omicsList")
  g$unnorm_gene_index <- g$gene_data_index
  g$gene_data_index <- g$data_index[1]

  for(i in 1:length(g$data_norm_index)){try({
    start_length <- length(g$omicsList) + 1
    g$omicsList[[start_length]] <- g$omicsList[[ g$data_norm_index[i] ]]
    g$omicsList[[start_length]][["eSet"]] <- g$omicsList[[start_length]][["siteNorm"]] 
    g$omicsList[[start_length]][["dataType"]] <- paste(g$omicsList[[start_length]][["dataType"]],
                                                     "_NormTo",g$omicsList[[ g$data_index[1] ]][["dataType"]],sep="")
    g$gene_data_index <- c(g$gene_data_index, start_length)
  }) }

  g$calls = c(g$calls, "g.use.site.norm")
  g
}

g.param.norm = function(g, deps=T){
  if(deps){
    g=g.run.deps(g, "g.norm.to.first")
    if(use_site_norm) g=g.run.deps(g, "g.use.site.norm")
  }
  g$calls = c(g$calls, "g.param.norm")
  g
}

g.varplot = function(g, deps=T){
  if(deps) g=g.run.deps(g, "g.param.norm")
  output_links<-"";

  # Pairs, Variation, and Corrplot
  for(i in 1:length(g$omicsList)){ try({
    g$omicsList[[i]][[6]] <- variationPlot(eset=g$omicsList[[i]][["eSet"]], type=g$omicsList[[i]][["dataType"]], outputpath=g$output_plots_path);
    names(g$omicsList[[i]])[6] <- "topVariable";

    # Make output hyperlinks for markdown
    type_name <- g$omicsList[[i]][["dataType"]]
    add_link <- paste("[ ",type_name, "-Variation ](", g$output_plots_subdir,"/variation_",type_name,".pdf) | ",
                      "[ ",type_name, "-Correlation ](", g$output_plots_subdir,"/corrplot_",type_name,".pdf) | ",
                      "[ ",type_name, "-Pairs ](",g$output_plots_subdir,"/pairs_",type_name,".png)",
                      sep="");

    output_links <- paste(g$output_links, add_link, sep="  \n" );
  }) }

  # Coefficient of Variation
  try({
    CVs <- data.frame( CV=apply( exprs(g$omicsList[[1]][["eSet"]]), 1, function(x) { (sd(x)/mean(x))*100 } ) )
    CVs$Dataset <- g$omicsList[[1]][["dataType"]]
    if( length(g$omicsList)>1 ){ for(i in 2:length(g$omicsList) ){
      tmp <- data.frame( CV=apply( exprs(g$omicsList[[i]][["eSet"]]), 1, function(x) { (sd(x)/mean(x))*100 } ) )
      tmp$Dataset <- g$omicsList[[i]][["dataType"]]
      CVs <- rbind(CVs, tmp)
    }}

    plot <- ggplot(CVs, aes(x=CV, fill=Dataset)) + geom_density(alpha=0.2)+ theme_bw() + 
      labs(title=paste("Coefficient of Variation Plot \n ", sep=''),
           x="Coefficient of Variation (%)", y="Frequency");

    output_filename<-file.path(g$output_plots_path, paste("CV_Plot",".pdf", sep=""))
    pdf(output_filename, width=3, height=3)
    print(plot+theme(legend.position="none"))
    suppressWarnings(print(plot+ scale_x_continuous(limits=c(0, 20))+theme(legend.position="none")) )
    suppressWarnings(print(plot+ scale_x_continuous(limits=c(0, 10))+theme(legend.position="none")) )
    gridExtra::grid.arrange(g_legend(plot))
    dev.off()

    avg <- aggregate( CVs[,"CV"], list(CVs[,"Dataset"]), mean) 
    colnames(avg) <- c("Dataset", "Average CV")
    print(avg)
    write.table(avg, file=file.path(g$output_files_path, "Average_CV.txt"), quote=F, sep="\t")

    add_link <- paste("[ CV Plot ](", g$output_plots_subdir,"/CV_Plot",".pdf) ", sep="");
    output_links <- paste(output_links, add_link, sep="  \n" );
  })

  g$calls = c(g$calls, "g.varplot")
  g$output_links = output_links
  g
}

g.combine.met = function(g, deps=T){
  if(deps) g=g.run.deps(g, "g.param.norm")
  if( length(g$mz_data_index) > 1 ){ # make combined eset for metabolomics
    intersect_rbind=function(xm,ym){
      inds=intersect(colnames(xm),colnames(ym))
      rbind(xm[,inds],ym[,inds])
    }
    merge_eset=function(xl,yl){
      x=xl$eSet; y=yl$eSet
      rn=c(paste(xl$dataType,rownames(exprs(x))),paste(yl$dataType,rownames(exprs(y))))
      rd=intersect_rbind(exprs(x),exprs(y))
      rownames(rd)=rn
      ret=ExpressionSet(assayData=rd)
      pData(ret)=pData(x)
      fr=intersect_rbind(fData(x),fData(y))
      rownames(fr)=rn
      fData(ret)=fr
      ret
    }
    newi=length(g$omicsList)+1
    g$mz_data_index=c(g$mz_data_index,newi)
    g$metab_data_index=c(g$metab_data_index,newi)
    g$omicsList[[newi]]=list(eSet=merge_eset(g$omicsList[[g$mz_data_index[1]]],g$omicsList[[g$mz_data_index[2]]]),dataType="met_combined",dataFormat="Metabolomics (combined)")
    g$omicsList[[newi]][["topVariable"]]=variationPlot(eset=g$omicsList[[newi]][["eSet"]], type=g$omicsList[[newi]][["dataType"]], outputpath=g$output_plots_path);
  }

  g$calls = c(g$calls, "g.combine.met")
  g
}


g.heatmap.static = function(g, deps=T){
  if(deps) g=g.run.deps(g, c("g.varplot","g.combine.met"))
  output_links<-"";
  # Heatmaps
  BiocParallel::bplapply(1:length(g$omicsList),FUN=function(i){ try({
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
    },silent=T) }
    # Draw the heatmaps and sace the Complex Heatmap object
    map_out <- drawHeatmaps(eset=g$omicsList[[i]][["eSet"]], emat_top=g$omicsList[[i]][["topVariable"]],
                            type=g$omicsList[[i]][["dataType"]], subset=subset_rows, k_clust=knn_heatmap,
                            outputpath=g$output_plots_path, outputcontrastpath=g$output_contrast_path);
    saveRDS(map_out, file=file.path(output_files_path, paste("Heatmap_",g$omicsList[[i]][["dataType"]],".RDS",sep="")) )

  }) })

  for(i in 1:length(g$omicsList)){
    # Make output hyperlinks for markdown
    add_link <- paste("[ ",g$omicsList[[i]][["dataType"]], " ](", g$output_plots_subdir,"/heatmaps_",g$omicsList[[i]][["dataType"]],".pdf)", sep="")
    output_links <- paste(output_links, add_link, sep=" | " )
  }

  g$calls = c(g$calls, "g.heatmap.static")
  g$output_links = output_links
  g
}


g.boxplots = function(g, deps=T){
  if(deps) g=g.run.deps(g, "g.param.norm")
  for(i in 1:length(g$omicsList)){ try({
    subset_rows=FALSE;
    if( class(g$subset_genes)!="logical" ){ try({
      g$subset_genes <- unique(unlist(g$subset_genes, use.names = F))
      subset_rows <- c()
      if( "Gene" %in% colnames(fData(g$omicsList[[i]][["eSet"]])) ){
        subset_rows <- rownames(g$omicsList[[i]][["eSet"]])[which(fData(g$omicsList[[i]][["eSet"]])$Gene %in% g$subset_genes )]
      }
      if( length(subset_rows)==0 ){
        subset_rows <- rownames(g$omicsList[[i]][["eSet"]])[ which(rownames(g$omicsList[[i]][["eSet"]]) %in% g$subset_genes ) ]
      }
    },silent=T) }
    # Draw the boxplots
    if(class(subset_rows)!="logical"){ 
      output_filename<-file.path(g$output_plots_path, paste(g$omicsList[[i]][["dataType"]],"_SelectFeatures",".pdf", sep=""))
      pdf(output_filename, width=3, height=3)
      for(j in 1:length(subset_rows)){ try({
        dat <- data.frame(Intensity=exprs(g$omicsList[[i]][["eSet"]])[subset_rows[j],],
                          Group = pData(g$omicsList[[i]][["eSet"]])$Group  )
        print(ggplot(data=dat, aes(color=Group, x=Group, y=Intensity)) + geom_boxplot() + geom_point() +
              theme_bw() + labs(title=paste( g$omicsList[[i]][["dataType"]],"\n", subset_rows[j], sep="")) +
              theme(legend.position = "none") + ylab("Log2 Intensity") +
              theme(axis.text.x=element_text(angle=45, hjust=1))
        )

      }) }
      dev.off()
    }
  }) }

  g$calls = c(g$calls, "g.boxplots")
  g
}

g.interactive.heatmap = function(g, deps=T){
  if(deps) g=g.run.deps(g, "g.param.norm")
  output_links<-"";

  for(i in 1:length(g$omicsList)){
    # Make the interactive heatmap
    interactiveHeatmap(eset=g$omicsList[[i]][["eSet"]], type=g$omicsList[[i]][["dataType"]], outputpath=g$output_plots_path);
    
    # Make output hyperlinks for markdown
    add_link <- paste("[ ",g$omicsList[[i]][["dataType"]], "-All ](", g$output_plots_subdir,"/heatmap_all_",
                      g$omicsList[[i]][["dataType"]],".html)", sep="");
    output_links <- paste(output_links, add_link, sep=" | " );
  }

  g$calls = c(g$calls, "g.interactive.heatmap")
  g$output_links = output_links
  g
}

g.save.data.qc = function(g, deps=T){
  if(deps) g=g.run.deps(g, "g.varplot")
  # save for omics integrator, expression matrix, expression set RDS, ranked list for GSEA
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
    saveFiles(data=g$omicsList[[i]], type=g$omicsList[[i]][["dataType"]], subset=subset_rows,
              outputpath=g$output_files_path, outputcontrastpath=g$output_contrast_path_files);
  }  

  # Save excel file summary for collaborators
  if(saveXlsx==TRUE){    
    wbOut <- openxlsx::createWorkbook()

    for(i in 1:length(g$omicsList)){
      if(g$omicsList[[i]][["dataType"]]!="met_combined"){
        eSet=g$omicsList[[i]][["eSet"]][,order(pData(g$omicsList[[i]][["eSet"]])$Group)]
        writeDataToSheets(wb=wbOut, eset=eSet, type=g$omicsList[[i]][["dataType"]], data_format=g$omicsList[[i]][["dataFormat"]]);
      }
    } 

    output_filename=file.path(g$output_path, paste(gsub("\\.","",make.names(project_name)), "_Summary", ".xlsx", sep=""))
    openxlsx::saveWorkbook(wbOut, file=output_filename, overwrite=TRUE)
  }

  g$calls = c(g$calls, "g.save.data.qc")
  g
}

