gsea_ranked_data=function(g, contrast_name, coef, rcols){
  do.call("rbind", lapply(1:length(g$gene_data_index),FUN=function(i){
    top_sum <- FALSE;
    try({ 
      top_sum <- limma::topTable(g$omicsList[[ g$gene_data_index[i] ]][["fit"]], adjust.method="BH", n=Inf, sort.by=if(is.null(coef)) 'B' else 'p', coef=coef);
      ranked <- cbind(top_sum[,"Gene"], rcols(top_sum)) 
    })
    if( class(top_sum)=="logical" ){ try({
      top_sum <- data.frame(logFC= fData(g$omicsList[[ g$gene_data_index[i] ]][["eSet"]])[,paste("logfc_",gsub("-","_",contrast_name),sep="")],
                        Gene=fData(g$omicsList[[ g$gene_data_index[i] ]][["eSet"]])$Gene, 
                        feature_identifier =fData(g$omicsList[[ g$gene_data_index[i] ]][["eSet"]])[,"feature_identifier"] );
      ranked <- top_sum[,c("Gene","logFC")]
    })}
    ranked
  }))
}

write_gsea_data=function(g, ranked, contrast_name, secondary_outfile){
  colnames(ranked)<-c("GeneName", "rank")
  ranked <- ranked[ranked[,"GeneName"]!="", ]
  ranked <- ranked[order(abs(as.numeric(ranked[,"rank"])), decreasing=TRUE),]
  ranked <- ranked[!duplicated(ranked[,"GeneName"]),]
  ranked <- ranked[order(as.numeric(ranked[,"rank"]), decreasing=TRUE),]
  output_filename <- file.path(g$output_contrast_path_files,paste("GSEA_",contrast_name,".rnk", sep=""));
  write.table(x=ranked, file=output_filename, sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE)
    
  if(is.function(secondary_outfile)){
    secondary_outfile(g, ranked,contrast_name);
  }
   
  if(!( grepl("Human", species)) ){ 
    output_filename <- file.path(g$output_contrast_path_files,paste0("GSEA_Uppercase_",contrast_name,".rnk"));
    ranked[,"GeneName"] <- toupper(ranked[,"GeneName"])
    write.table(x=ranked, file=output_filename, sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    if(is.function(secondary_outfile)){
      secondary_outfile(g, ranked, paste0("Uppercase_",contrast_name));
    }
  }
}

gen_gsea_file=function(g, contrast_name, coef, rcols, secondary_outfile){
  write_gsea_data(g, gsea_ranked_data(g, contrast_name, coef, rcols), paste0("combined_",contrast_name), secondary_outfile)
}
 
contrast_secondary=function(g, ranked, contrast_name){
    output_filename <- file.path(g$output_contrast_path_files,paste0("GSEA_combined_",contrast_name,"_AbsVal.rnk"));
    ranked_pval<- ranked
    ranked_pval[,"rank"] <- abs(as.numeric(ranked[,"rank"]))
    write.table(x=ranked_pval,file=output_filename, sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE);
}
 
#' GSEA Combined Contrast
#'
#' Save combined ranked lists for GSEA, for each contrast
#' 
#' @param g omics notebook state, created by g.notebook.setup()
#' @param deps automatically run dependencies?
#'
#' @return updated omics notebook state
#' 
#' @export
g.gsea.combined.contrast = function(g, deps=T){
  if(deps) g=g.run.deps(g, c("g.fc", "g.limma"))

  cnames=gsub("-","_",g$contrast_strings)
  for (j in 1:length(g$loop_list) ){ try({
    gen_gsea_file(g, cnames[j], g$loop_list[j], function(top_sum){sign(top_sum[,"logFC"]) * -log10(top_sum[,"P.Value"])},contrast_secondary)
  })}
  if (g$time_index>0) gen_gsea_file(g, "Timecourse_Overall", g$time_start:g$time_end, function(top_sum){top_sum[,"F"]},NULL)
  if (g$statistic_index=="F") gen_gsea_file(g, "Overall", NULL, function(top_sum){top_sum[,"F"]},NULL)

  g$calls = c(g$calls, "g.gsea.combined.contrast")
  g
}

#' Abundance GSEA
#'
#' Generate ranked lists for GSEA based on absolute abundance instead of differential analysis
#' 
#' @param g omics notebook state, created by g.notebook.setup()
#' @param deps automatically run dependencies?
#'
#' @return updated omics notebook state
#' 
#' @export
g.abundance.gsea = function(g, deps=T){
  if(deps) g=g.run.deps(g, c("g.fc", "g.limma"))

  g$loop_list <- c(1);
  g$output_contrast_subdir <- file.path(g$output_subdir,"Abundance")
  g$output_contast_subdir_files <- file.path(g$output_contrast_subdir, "files")
  g$output_contrast_path <- file.path(g$working_dir, g$output_contrast_subdir)
  g$output_contrast_path_files <- file.path(g$output_contrast_path, "files")
  if( dir.exists(g$output_contrast_path) == FALSE ) { dir.create(g$output_contrast_path) }
  if( dir.exists(g$output_contrast_path_files) == FALSE ) { dir.create(g$output_contrast_path_files) }
  #g$num_contrasts <- 1
  g$contrast_strings <- "Abundance"
  g$contrast_strings_file <- gsub("-","_",g$contrast_strings)

  for(i in 1:length(g$omicsList) ){
    if( "Gene" %in% colnames(fData(g$omicsList[[i]][["eSet"]])) | "mz" %in% colnames(fData(g$omicsList[[i]][["eSet"]])) ){
      if( "Gene" %in% colnames(fData(g$omicsList[[i]][["eSet"]])) ){
        top_sum <- data.frame(logFC= rowMeans(exprs(g$omicsList[[i]][["eSet"]])),
                              Gene=fData(g$omicsList[[ g$gene_data_index[i] ]][["eSet"]])$Gene, 
                              feature_identifier =fData(g$omicsList[[ g$gene_data_index[i] ]][["eSet"]])[,"feature_identifier"] );
        ranked <- top_sum[,c("Gene","logFC")]
      }
      if( "mz" %in% colnames(fData(g$omicsList[[i]][["eSet"]])) ){
        top_sum <- data.frame(logFC= rowMeans(exprs(g$omicsList[[i]][["eSet"]])),
                              Gene=fData(g$omicsList[[i]][["eSet"]])$mz, 
                              feature_identifier =fData(g$omicsList[[i]][["eSet"]])[,"feature_identifier"] );
        ranked <- top_sum[,c("Gene","logFC")]

      }

      species_backup=species
      ## XXX should this check be done when there are multiple contrasts too?
      if(!("Gene" %in% colnames(fData(g$omicsList[[i]][["eSet"]])))) species="xx";
      write_gsea_data(ranked, paste(g$omicsList[[i]][["dataType"]],g$contrast_strings,sep="_"), NULL)
      species=species_backup
    } 
  }

  g$calls = c(g$calls, "g.abundance.gsea")
  g
}

#' GSEA Data by Paramters
#'
#' Convenience function to generate the GSEA data appropriate for the current analysis parameters
#' 
#' @param g omics notebook state, created by g.notebook.setup()
#' @param deps automatically run dependencies?
#'
#' @return updated omics notebook state
#' 
#' @export
g.param.gsea.data = function(g, deps=T){
  if(deps){
    g=g.run.deps(g, "g.make.contrasts")
    if(length(g$contrastgroups)>1){
      dep = "g.savedata" #differential ranked files
      if(length(g$gene_data_index)>1){ # combined ranked files
        dep = c(dep, "g.gsea.combined.contrast")
      }
    } else {
      dep = "g.abundance.gsea" # abundance based ranked files
    }
    g=g.run.deps(g, dep)
  }
  g$calls = c(g$calls, "g.param.gsea.data")
  g
}

#' GSEA Prep
#'
#' Preprocessing the data in preparation for GSEA
#' 
#' @param g omics notebook state, created by g.notebook.setup()
#' @param deps automatically run dependencies?
#'
#' @return updated omics notebook state
#' 
#' @export
g.gsea.prep = function(g, deps=T){
  if(deps) g=g.run.deps(g, "g.param.gsea.data")

  case_modifier <- if(!grepl("Mouse|Human", species)) "_Uppercase_" else "_"

  analysis_names <- vector("list", (length(g$loop_list)*length(g$gene_data_index)) ) ;
  counter<-1;                       
  for (j in 1:length(g$gene_data_index)) {
    for (k in 1:length(g$loop_list)) {
      analysis_names[[counter]] <- paste0(g$omicsList[[ g$gene_data_index[j] ]][["dataType"]], "_", g$contrast_strings_file[k]) ;
      counter <- counter+1;
    }
  }
  rnk_files <- vector("list", (length(g$loop_list)*length(g$gene_data_index)) ) ;
  counter<-1;                       
  for (j in 1:length(g$gene_data_index)) {
    for (k in 1:length(g$loop_list)) {
      rnk_files[[counter]] <- file.path(paste0(g$output_contrast_path_files,"/GSEA_",
                                                g$omicsList[[ g$gene_data_index[j] ]][["dataType"]],case_modifier,
                                                g$contrast_strings_file[k],".rnk"))
      counter <- counter+1;
    }
  }

  if(length(g$gene_data_index)>1 & length(g$contrastgroups)>1){ 
    for (i in 1:length(g$loop_list)) {
      analysis_names[[length(analysis_names)+1]] <- paste0("Combined_",g$contrast_strings_file[i])
    }
    for (i in 1:length(g$loop_list)){
      rnk_files <- append(rnk_files, file.path(g$output_contrast_path_files,
                                                paste0("GSEA",case_modifier,"combined_",g$contrast_strings_file[i],".rnk") ))
    }
  }

  g$calls = c(g$calls, "g.gsea.prep")
  g$analysis_names = analysis_names
  g$rnk_files = rnk_files
  g
}


#' GSEA
#'
#' Run GSEA
#' 
#' @param g omics notebook state, created by g.notebook.setup()
#' @param gmt path to GMT files to 
#' @param gmt_name concise names for the GMT files, for labeling results
#' @param working_dir the name of the directory to write outputs to
#' @param deps automatically run dependencies?
#'
#' @return updated omics notebook state
#' 
#' @export
g.gsea = function(g, gmt=NULL, gmt_name=NULL, working_dir="3_GSEA", deps=T){
  if(deps) g=g.run.deps(g, "g.gsea.prep")

  try({
    gsea_working_path <- file.path(g$output_path, working_dir)
    if( dir.exists(gsea_working_path) == FALSE ) { dir.create(gsea_working_path) }

    if(!is.null(gmt)){
      dest_gmt_file <- gmt
    } else if( species!="Other" ) {
      g$species_gmt <- fetchGMT(g$gmt_path, species)
      dest_gmt_file <- g$species_gmt
    } else {
      return(g)
    }

    if(is.null(gmt_name)){ gmt_name = dest_gmt_file}

    for (i in 1:length(g$analysis_names)){
      analysis_name <- paste(g$analysis_names[[i]], basename(gmt_name), sep="_")
      rnk_file <- g$rnk_files[[i]]
      runGSEA(rnk=rnk_file, gmt=dest_gmt_file, analysisName=analysis_name, pathway_plots=T,out_path=gsea_working_path )
    }
  })  

  g$calls = c(g$calls, "g.gsea")
  g$gsea_working_dir = working_dir
  g
}

#' Custom GSEA
#'
#' Run GSEA with custom gene sets from the Gene_Sets directory
#' 
#' @param g omics notebook state, created by g.notebook.setup()
#' @param working_dir the name of the directory to write outputs to
#' @param run.GSEA actually run GSEA (T) or just adjust 'g' with the custom parameters (F)?
#' @param deps automatically run dependencies?
#'
#' @return updated omics notebook state
#' 
#' @export
g.gsea.custom = function(g, working_dir="3_GSEA_Custom", run.GSEA=F, deps=T){
  if(run.GSEA && deps) g=g.run.deps(g, "g.gsea.prep")

  gmt_dirs <- c("Gene_Sets",file.path("..","Gene_Sets"))
  g$gmt_files_cust <- list.files(gmt_dirs,".*.GMT",full.names=T,ignore.case=T)
  g$gmt_names_cust <- sub(".gmt$","",basename(g$gmt_files_cust),ignore.case=T)

  if(run.GSEA){
    for(i in 1:length(g$gmt_files_cust)){
      g = g.gsea(g, gmt=g$gmt_files_cust[[i]], gmt_name=g$gmt_names_cust[[i]], working_dir=working_dir)
    }
  }

  g$calls = c(g$calls, "g.gsea.custom")
  g
}

#' GSEA MOMENTA
#'
#' MOMENTA analysis
#' 
#' @param g omics notebook state, created by g.notebook.setup()
#' @param working_dir the name of the directory to write outputs to
#' @param deps automatically run dependencies?
#'
#' @return updated omics notebook state
#' 
#' @export
g.gsea.momenta = function(g, working_dir="3_GSEA_MOMENTA", deps=T){
  if(deps) g=g.run.deps(g, "g.gsea.prep")

  geneset_lookup <-  read.delim(get.data.fileconn(file.path("matched","MatchedGeneSets.txt")))
  
  g$gmt_files <- as.vector(geneset_lookup[which(geneset_lookup[,"Species"]==species),"GeneSet"])
  g$gmt_names <- as.vector(geneset_lookup[which(geneset_lookup[,"Species"]==species),"Name"])
  g$gmt_files = lapply(file.path("matched",g$gmt_files),FUN=get.data.fileconn)

  for(i in 1:length(g$gmt_files)){
    g = g.gsea(g, gmt=g$gmt_files[[i]], gmt_name=g$gmt_names[i], working_dir=working_dir)
  }

  g$calls = c(g$calls, "g.gsea.momenta")
  g
}

#' EnrichR
#'
#' Run EnrichR
#' 
#' @param g omics notebook state, created by g.notebook.setup()
#' @param working_dir the name of the directory to write outputs to
#' @param deps automatically run dependencies?
#'
#' @return updated omics notebook state
#' 
#' @export
g.enrichr = function(g, working_dir="3_EnrichR", deps=T){
  if(deps) g=g.run.deps(g, c("g.limma"))

  enrichR:::.onAttach() # load webservice
  run_enrichr_seperate <- TRUE;

  enrichr_working_path <- file.path(g$output_path, working_dir)
  if( dir.exists(enrichr_working_path) == FALSE ) { dir.create(enrichr_working_path) }

  # Run Enrichr for each data type
  for(i in 1:length(g$gene_data_index) ){
    sig_cutoff <- round(sig_percent*nrow(g$omicsList[[ g$gene_data_index[i] ]][["eSet"]]), digits=0)
    for(j in 1:length(g$loop_list) ){ try({
      top_sum<- FALSE;
      contrast_name <- g$contrast_strings[j]
      try({
        top_sum <- limma::topTable(g$omicsList[[ g$gene_data_index[i] ]][["fit"]], adjust="BH", n=Inf, sort.by='p', coef=g$loop_list[j], p.value=adjpcutoff);
        if (nrow(top_sum)<sig_cutoff){ top_sum <- limma::topTable(g$omicsList[[ g$gene_data_index[i] ]][["fit"]], adjust="BH", n=sig_cutoff, sort.by='p', coef=g$loop_list[j]); }
      })
      if( class(top_sum)=="logical" ){
        top_sum <- data.frame(logFC= fData(g$omicsList[[ g$gene_data_index[i] ]][["eSet"]])[,paste0("logfc_",g$contrast_strings_file[j])],
                              Gene=fData(g$omicsList[[ g$gene_data_index[i] ]][["eSet"]])$Gene, 
                              feature_identifier =fData(g$omicsList[[ g$gene_data_index[i] ]][["eSet"]])[,"feature_identifier"] );
        top_sum <- top_sum[top_sum[,"Gene"]!="", ]
        top_sum <- top_sum[order(abs(as.numeric(top_sum[,"logFC"])), decreasing=TRUE),]
        if( nrow(top_sum)>sig_cutoff ){  top_sum <- top_sum[1:sig_cutoff,] }
        contrast_name <- gsub("logfc_","",g$logfc_index[j])
      }
      type_name <- paste(g$omicsList[[ g$gene_data_index[i] ]][["dataType"]], contrast_name, sep="_");
      
      runEnrichR(genes=top_sum[, c("Gene", "logFC") ], type=type_name, run_seperate=run_enrichr_seperate,
                 search_dat=g$search_databases, outputpath=enrichr_working_path) 

    }) } 
    
    # Timecourse
    if(g$time_index>0){try({
      top_sum<- FALSE;
      contrast_name <- "TimeCourse_Overall"
      top_sum <- limma::topTable(g$omicsList[[ g$gene_data_index[i] ]][["fit"]], adjust="BH", n=Inf, coef=g$time_start:g$time_end, p.value=adjpcutoff);
      if (nrow(top_sum)<sig_cutoff){ top_sum <- limma::topTable(g$omicsList[[ g$gene_data_index[i] ]][["fit"]], adjust="BH", n=sig_cutoff, coef=g$time_start:g$time_end); }
      top_sum[,"logFC"] <- top_sum[,'F']
      type_name <- paste(g$omicsList[[ g$gene_data_index[i] ]][["dataType"]], "TimeCourse_Overall", sep="_");
      
      runEnrichR(genes=top_sum[, c("Gene", "logFC") ], type=type_name, run_seperate=run_enrichr_seperate,
                 search_dat=g$search_databases, outputpath=enrichr_working_path) 

    }) } 
  } 
    
  # Run enrichr for combined data
  try({
  if( length(g$gene_data_index)>1 ){
    for(j in 1:length(g$loop_list)){
      type_name <- paste("Combined",g$contrast_strings[j], sep="_");

      gene_table = do.call("rbind",lapply(1:length(g$gene_data_index),FUN=function(i){
        sig_cutoff <- round(sig_percent*nrow(g$omicsList[[ g$gene_data_index[1] ]][["eSet"]]), digits=0)
        gene_table <- limma::topTable(g$omicsList[[ g$gene_data_index[i] ]][["fit"]], adjust="BH", n=Inf, sort.by='p', coef=g$loop_list[j], p.value=adjpcutoff)
        if (nrow(gene_table)<sig_cutoff){ gene_table<-limma::topTable(g$omicsList[[ g$gene_data_index[i] ]][["fit"]], adjust="BH",
                                                              n=sig_cutoff, sort.by='p', coef=g$loop_list[j]) }
        gene_table[, c("Gene", "logFC") ];
      }))
          
      runEnrichR(gene_table[, c("Gene", "logFC") ], type=type_name, run_seperate=run_enrichr_seperate,
                 search_dat=g$search_databases, outputpath=enrichr_working_path )
    } 

    if(g$statistic_index=='F'){ 
      type_name <- paste("Combined","Overall", sep="_");

      gene_table = do.call("rbind",lapply(1:length(g$gene_data_index),FUN=function(i){
        sig_cutoff <- round(sig_percent*nrow(g$omicsList[[ g$gene_data_index[1] ]][["eSet"]]), digits=0)
        gene_table <- limma::topTable(g$omicsList[[ g$gene_data_index[1] ]][["fit"]], adjust="BH", n=Inf, sort.by='B', p.value=adjpcutoff)
        if (nrow(gene_table)<sig_cutoff){ gene_table<-limma::topTable(g$omicsList[[ g$gene_data_index[1] ]][["fit"]], adjust="BH", n=sig_cutoff, sort.by='B') }
        gene_table[, c("Gene", "F") ];
      }))
          
      gene_table[,"logFC" ] <- gene_table[,"F" ] 
      runEnrichR(gene_table[, c("Gene", "logFC") ], type=type_name, run_seperate=run_enrichr_seperate,
                 search_dat=g$search_databases, outputpath=enrichr_working_path )
    }
  }
  })

  # Add enrichr based on F stat
  for(i in 1:length(g$gene_data_index) ){ 
    if(g$statistic_index =='F'){ try({
      top_sum<- FALSE;
      contrast_name <- "F-statistic"
      sig_cutoff <- round(sig_percent*nrow(g$omicsList[[ g$gene_data_index[1] ]][["eSet"]]), digits=0)
      top_sum <- limma::topTable(g$omicsList[[ g$gene_data_index[i] ]][["fit"]], adjust="BH", n=Inf, sort.by='B', p.value=adjpcutoff);
      if (nrow(top_sum)<sig_cutoff){ top_sum <- limma::topTable(g$omicsList[[ g$gene_data_index[i] ]][["fit"]], adjust="BH", n=sig_cutoff, sort.by='B'); }
      top_sum[,"logFC"] <- top_sum[,'F']
      type_name <- paste(g$omicsList[[ g$gene_data_index[i] ]][["dataType"]], 'F-statistic', sep="_");
      
      runEnrichR(genes=top_sum[, c("Gene", "logFC") ], type=type_name, run_seperate=run_enrichr_seperate,
                 search_dat=g$search_databases, outputpath=enrichr_working_path) 
    }) }
  } 

  g$calls = c(g$calls, "g.enrichr")
  g$enrichr_working_dir = working_dir
  g
}

#' EnrichmentMap
#'
#' Run EnrichmentMap (currently unused)
#' 
#' @param g omics notebook state, created by g.notebook.setup()
#' @param deps automatically run dependencies?
#'
#' @return updated omics notebook state
#' 
#' @export
g.enrichmentmap = function(g, deps=T){
  #unused??
  g$calls = c(g$calls, "g.enrichmentmap")
  g
}

#' Kinase-Substrate Enrichment Analysis (KSEA)
#'
#' Run KSEA for phosphoproteomics
#' 
#' @param g omics notebook state, created by g.notebook.setup()
#' @param working_dir the name of the directory to write outputs to
#' @param deps automatically run dependencies?
#'
#' @return updated omics notebook state
#' 
#' @export
g.ksea = function(g, working_dir="3_KSEA", deps=T){
  if(deps) g=g.run.deps(g, c("g.limma"))
  output_links <- "";

  dir_path <- file.path(g$output_path, working_dir)
  if( dir.exists(dir_path) == FALSE ) {dir.create(dir_path)}

  try({
  for( i in 1:length(g$phos_data_index) ){
    
    for( j in 1:length(g$loop_list) ){
      PX2 <- FALSE
      try({
        top_sum <- limma::topTable(g$omicsList[[ g$phos_data_index[i] ]][["fit"]], adjust="BH", n=Inf, sort.by='p', coef=j);
        PX2 <- top_sum[,c("Protein","Gene","Sequence.window")]
        colnames(PX2) <- c("Protein", "Gene", "Peptide")
        PX2[,"Residue.Both"] <- paste0(top_sum[,"Amino.acid"],top_sum[,"Position"])
        PX2[,"p"] <- top_sum[,"P.Value"]
        PX2[,"FC"] <- 2^top_sum[,"logFC"]
      })
      if( class(PX2)=="logical" ){
        PX2 <- fData(g$omicsList[[ g$phos_data_index[i] ]][["eSet"]])[,c("Protein","Gene","Sequence.window")]
        colnames(PX2) <- c("Protein", "Gene", "Peptide")
        PX2[,"Residue.Both"] <- paste0(fData(g$omicsList[[ g$phos_data_index[i] ]][["eSet"]])[,"Amino.acid"],
                                      fData(g$omicsList[[ g$phos_data_index[i] ]][["eSet"]])[,"Position"])
        PX2[,"p"] <- "NULL";
        contrast_grab <- grep( paste0("logfc_", g$contrast_strings_file[j]), colnames(fData(g$omicsList[[ g$phos_data_index[i] ]][["eSet"]]) ) )
        PX2[,"FC"] <- 2^(fData(g$omicsList[[ g$phos_data_index[i] ]][["eSet"]])[,contrast_grab] )    
      }
      
      dir_name <- paste0("KSEA_", g$contrast_strings[j])
      dir_path2 <- file.path(dir_path, dir_name)
      if( dir.exists(dir_path2) == FALSE ) {dir.create(dir_path2)}
      
      try({
        KSEAapp::KSEA.Complete(KSEAapp::KSData, PX2, NetworKIN=TRUE, NetworKIN.cutoff=5, m.cutoff=3, p.cutoff=0.05, output_dir=dir_path2)
      
      output_links <- paste0(output_links,paste("[ ",dir_name," ](",g$output_subdir,"/3_KSEA/",dir_name,"/KSEA Bar Plot.png)  |  ",sep=""))
      })
    }
  }
  })

  g$calls = c(g$calls, "g.ksea")
  g$ksea_output_links = output_links
  g
}

#' Metabolomics Enrichment
#'
#' Run MetaboanalystR's Peak Set Enrichment Analysis (PSEA) with mummichog identification
#' 
#' @param g omics notebook state, created by g.notebook.setup()
#' @param working_dir the name of the directory to write outputs to
#' @param deps automatically run dependencies?
#'
#' @return updated omics notebook state
#' 
#' @export
g.metabo.enrich = function(g, working_dir="3_Metabo_Enrichment", deps=T){
  if(deps) g=g.run.deps(g, c("g.make.contrasts","g.limma"))

  output_mummichog_reldir <- file.path(g$output_subdir, working_dir)
  output_mummichog_absdir <- file.path(g$working_dir, output_mummichog_reldir)
  if( dir.exists(output_mummichog_absdir) == FALSE ) { dir.create(output_mummichog_absdir) }

  geneset_lookup <-  read.delim(get.data.fileconn("metabolite_model_table.txt"))
  mumm_libs <- as.vector(geneset_lookup[which(geneset_lookup[,"Species"]==species),"MetLib"])

  suppressWarnings({ suppressMessages({
  for(i in 1:length(g$metab_data_index)){ try({
    for(mumm_lib_i in 1:length(mumm_libs)) {
      for (j in 1:length(g$loop_list)){ try({
        if(g$omicsList[[ g$metab_data_index[i] ]][["dataType"]] == "met_combined") next;

        analysis_name <- paste0(mumm_libs[mumm_lib_i],"_",g$omicsList[[ g$metab_data_index[i] ]][["dataType"]],
                              "_",g$contrast_strings_file[j])
        
        if(g$contrast_strings[1]!="Abundance"){ 
          output_mummichog_reldir_tmp <- file.path(output_mummichog_reldir, analysis_name)
          output_mummichog_path_tmp <- file.path(g$working_dir, output_mummichog_reldir_tmp)
          if( dir.exists(output_mummichog_path_tmp) == FALSE ) { dir.create(output_mummichog_path_tmp) }
          setwd(output_mummichog_path_tmp)
          
          input_data <- na.omit(limma::topTable(g$omicsList[[ g$metab_data_index[i] ]][["fit"]],
                                        adjust="BH", n=5000, sort.by='p', coef=j)[,c("mz", "P.Value", "t")] );
          colnames(input_data) <- c("m.z", "p.value", "t.score")
          output_filename <- file.path(output_mummichog_absdir,paste0("mummichog_",g$omicsList[[ g$metab_data_index[i] ]][["dataType"]],
                                                                  "_",g$contrast_strings_file[j],".txt"));
          write.table(x=input_data, file=output_filename, sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE)
        
          if( grepl("Negative.Ion", g$omicsList[[ g$metab_data_index[i] ]][["dataFormat"]] ) ) { ion_mode<-"negative" } else { ion_mode<-"positive"}
        
          mSet <- MetaboAnalystR::InitDataObjects("mass_all", "mummichog", FALSE);
          #SetPeakFormat("mpt")
          mSet <- MetaboAnalystR::UpdateInstrumentParameters(mSet, 0.1, ion_mode)
          mSet <- MetaboAnalystR::Read.PeakListData(mSet, output_filename);
          mSet <- MetaboAnalystR::SanityCheckMummichogData(mSet);
          mSet <- MetaboAnalystR::SetPeakEnrichMethod(mSet, "integ");
          mSet <- MetaboAnalystR::SetMummichogPval(mSet, 0.05);
          
          mSet <- MetaboAnalystR::PerformPSEA(mSet, mumm_libs[mumm_lib_i], "current", permNum=1000);
        
          mSet <- MetaboAnalystR::PlotIntegPaths(mSet,imgName="mummichog_integ_paths");
          #PreparePDFReport(mSet, usrName=analysis_name);
        
          tmp_table <- read.delim("mummichog_fgsea_pathway_enrichment.csv", sep=",", stringsAsFactors=FALSE)
          file.remove("mummichog_fgsea_pathway_enrichment.csv")
          tmp_table[,1] <- toupper(tmp_table[,1])
          out_table <- tmp_table[,c(1,1, 4, 5)]
          names(out_table) <- c("Term", "Description", "p.Val", "FDR")
          out_table[,"Phenotype"] <- "+1"
          out_table[tmp_table[,"NES"]<0,"Phenotype"] <- "-1"
          out_table[,"Gene_Hits"] <- tmp_table[,"Hits"]
          out_table[,"Gene_Total"] <- tmp_table[,"Pathway_Total"]
          out_table[,"Gene_Ratio"] <- out_table[,"Gene_Hits"]/out_table[,"Gene_Total"]
          out_table[,"NES"] <- tmp_table[,"NES"]
          out_table[,"Ionization"] <- ion_mode
          write.table(out_table,file="emap_fgsea_pathway_enrichment.txt", sep="\t", quote=FALSE, row.names=FALSE)
        
          if(nrow(out_table)>15){ out_table <- out_table[1:15,] }
          drawEnrichment(out_table,type=analysis_name,label="fgsea", outputpath=output_mummichog_path_tmp)
        
          tmp_table <- read.delim("mummichog_pathway_enrichment.csv", sep=",", stringsAsFactors=FALSE)
          file.remove("mummichog_pathway_enrichment.csv")
          tmp_table[,1] <- toupper(tmp_table[,1])
          out_table <- tmp_table[,c(1,1, 7)]
          names(out_table) <- c("Term", "Description", "p.Val")
          #out_table[,"Phenotype"] <- "+1"
          out_table[,"Gene_Hits"] <- tmp_table[,"Hits.total"]
          out_table[,"Gene_Total"] <- tmp_table[,"Pathway.total"]
          out_table[,"Gene_Ratio"] <- out_table[,"Gene_Hits"]/out_table[,"Gene_Total"]
          out_table[,"Ionization"] <- ion_mode
          write.table(out_table,file="emap_mummi_pathway_enrichment.txt", sep="\t", quote=FALSE, row.names=FALSE)
        
          #if(nrow(out_table)>15){ out_table <- out_table[1:15,] }
          #drawEnrichment(out_table,type=analysis_name,label="mummichog", outputpath=output_mummichog_path_tmp)
          #mset<-PlotIntegPaths(mSet, "png", 400, width=10)
        
        } else { # Run based on enrichment of most intense/least intense
          input_data <- data.frame(m.z= fData(g$omicsList[[ g$metab_data_index[i] ]][["eSet"]])$mz,
                                  rank= rowMeans(exprs(g$omicsList[[ g$metab_data_index[i] ]][["eSet"]])) )
          input_data <- input_data[order(input_data[,"rank"], decreasing=T), ]
          output_filename <- file.path(output_mummichog_absdir,paste0("mummichog_",g$omicsList[[ g$metab_data_index[i] ]][["dataType"]],
                                                                  "_",g$contrast_strings_file[j],".txt"));
          write.table(x=input_data[,"m.z"], file=output_filename, sep='\t',row.names=FALSE, col.names=FALSE, quote=FALSE)
          
        }
      }) }
    }
  }) }
  }) })
      
  setwd(g$working_dir)

  g$calls = c(g$calls, "g.metabo.enrich")
  g$mumm_libs = mumm_libs
  g$mumm_working_dir = output_mummichog_reldir
  g
}
