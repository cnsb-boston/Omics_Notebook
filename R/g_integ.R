#' Active Pathways
#'
#' Run Active Pathways multiomics pathway enrichment analysis
#' 
#' @param g omics notebook state, created by g.notebook.setup()
#' @param working_dir the name of the directory to write outputs to
#' @param deps automatically run dependencies?
#'
#' @return updated omics notebook state
#' 
#' @export
g.activepath = function(g, working_dir="4_ActivePathways", deps=T){
  if(deps) g=g.run.deps(g, c("g.gsea.custom","g.limma","g.gsea.momenta"))

  gsea_working_path <- file.path(g$output_path, working_dir)
  if( dir.exists(gsea_working_path) == FALSE ) { dir.create(gsea_working_path) }

  if(is.null(g$species_gmt) && species!="Other"){ g$species_gmt <- fetchGMT(g$gmt_path, species) }
  a_pathways_gmt <- g$species_gmt # large generic gene set
  if( length(g$gmt_files_cust) > 0 ){ a_pathways_gmt <- c(a_pathways_gmt, gmt_files_cust) } # Add Custom gene sets if there
  a_paths_gmt_names <- gsub("(.*/\\s*(.*$))","\\2", gsub(".gmt|.GMT","",a_pathways_gmt) )
  if( length(g$gmt_files) > 0 ){ # Add MOMENTA gene sets
    a_pathways_gmt <- c(a_pathways_gmt, g$gmt_files)
    a_paths_gmt_names  = c(a_paths_gmt_names, g$gmt_names)
  }

  write_ap_file=function(ttable,outname){
    merge_scores <- function(x,y) merge(x, y, by="Gene", all=T)

    scores <- Reduce(merge_scores, lapply(1:length(g$gene_data_index),function(i){
      dat <- ttable(g$omicsList[[ g$gene_data_index[i] ]][["fit"]])
      colnames(dat) <- c("Gene", g$omicsList[[ g$gene_data_index[i] ]][["dataType"]] )
      dat <- dat[order(dat[,2], decreasing = F),]
      dat[!duplicated(dat[,1]),]
    }))

    scores[is.na(scores)] <- 1
    row.names(scores)<- scores[,"Gene"]
    scores <- scores[,-1]
    scores <- as.matrix(scores)

    # run active pathways for each gmt_file
    for( gmt_i in 1:length(a_pathways_gmt) ){
      out_table <- ActivePathways::ActivePathways(scores, a_pathways_gmt[[gmt_i]],geneset.filter = c(10, 500))
      out_table <- out_table[ order(out_table$adjusted.p.val, decreasing = F),]
      out_table$overlap <- apply(out_table, 1, function(x) paste(as.character(unlist(x["overlap"])), collapse=","))
      out_table$evidence <- apply(out_table, 1, function(x) paste(as.character(unlist(x["evidence"])), collapse=","))
      gene_columns <- colnames(out_table)[grepl("Genes_", colnames(out_table))]
      for(g_column in gene_columns){
        out_table[,g_column ] <- apply(out_table, 1, function(x) paste(as.character(unlist(x[ g_column ])), collapse=","))
      } 
      out_filename <- file.path(gsea_working_path,paste0(outname,"_",substr(a_paths_gmt_names[gmt_i],1,20),gmt_i,".txt") )
      write.table(out_table, file=out_filename, row.names = F, col.names = T, quote = F, sep = "\t")
    }
  }

  #perform for each contrast
  for(i in 1:length(g$loop_list)) { try({
    contrast_ttable=function(dat) limma::topTable(dat, adjust="BH", n=Inf, sort.by='p', coef=g$loop_list[i])[,c("Gene", "P.Value")]
    write_ap_file(contrast_ttable,g$contrast_strings_file[i])
  }) }

  #repeat for F statistic
  if(g$statistic_index=='F'){ try({
    fstat_ttable=function(dat) limma::topTable(dat, adjust="BH", n=Inf,sort.by='F')[,c("Gene", "P.Value")]
    write_ap_file(fstat_ttable,"Fstatistic")
  }) }

  g$calls = c(g$calls, "g.activepath")
  g$gsea_working_dir = working_dir
  g
}

#' Mummichog to omicsList
#'
#' Add mummichog ID columns to the state variable 'g', and write them to _ID.txt files
#' 
#' @param g omics notebook state, created by g.notebook.setup()
#' @param deps automatically run dependencies?
#'
#' @return updated omics notebook state
#' 
#' @export
g.mumm2omicsList = function(g, deps=T){
  if(deps) g=g.run.deps(g, c("g.metabo.enrich"))

  iters=do.call("rbind",
    lapply(1:length(g$mumm_libs),function(mum){
      t(sapply(1:length(g$metab_data_index),function(mi){
        c(mum,mi)
      }))
    })
  )
  colnames(iters)=c("mumm_lib_index","met_ind")

  BiocParallel::bplapply(1:nrow(iters),FUN=function(ii){
    mumm_lib_i=iters[ii,"mumm_lib_index"]
    i=iters[ii,"met_ind"]
    try({
        analysis_name <- paste(g$mumm_libs[mumm_lib_i],"_",g$omicsList[[ g$metab_data_index[i] ]][["dataType"]], sep="" )
        
        output_mummichog_reldir_tmp <- file.path(g$mumm_working_dir, paste0(analysis_name, "_",g$contrast_strings_file[1]))
        output_mummichog_absdir_tmp <- file.path(working_dir, output_mummichog_reldir_tmp)
        met_ids <- read.delim(file.path(output_mummichog_absdir_tmp, "mummichog_matched_compound_all.csv"),
                              header=TRUE, stringsAsFactors = F, sep=",")
        
        id_column_name <- paste("mummichogID_", mumm_libs[mumm_lib_i], sep="")
        fData(g$omicsList[[ g$metab_data_index[i] ]][["eSet"]])[,id_column_name] <- ""
        
        for(n in 1:nrow(fData(g$omicsList[[ g$metab_data_index[i] ]][["eSet"]])) ){
          mz_value <- fData(g$omicsList[[ g$metab_data_index[i] ]][["eSet"]])[n,"mz"]
          met_id_list <- paste(met_ids[which(met_ids[,"Query.Mass"] == mz_value), "Matched.Compound"], collapse=",")
          fData(g$omicsList[[ g$metab_data_index[i] ]][["eSet"]])[n,id_column_name] <- met_id_list
        }
        
        output_filename <- file.path(output_files_path, paste(analysis_name,"_IDs.txt",sep="") )
        write.table(fData(g$omicsList[[ g$metab_data_index[i] ]][["eSet"]])[,c("feature_identifier", "mz", id_column_name)],
                    file=output_filename, quote=F, sep="\t", row.names=F)
      
    })
  })

  g$calls = c(g$calls, "g.mumm2omicsList")
  g
}

#' Combine Rank
#'
#' Combine the rank files for gene and metabolite based GSEA analyses, for later integrative analysis.
#' 
#' @param g omics notebook state, created by g.notebook.setup()
#' @param deps automatically run dependencies?
#'
#' @return updated omics notebook state
#' 
#' @export
g.combine.rnk = function(g, deps=T){
  if(deps) g=g.run.deps(g, c("g.metabo.enrich", "g.param.gsea.data"))

  rnk_files <-vector("list", g$num_contrasts)

  for (ci in 1:g$num_contrasts){
    for(mumm_lib_i in 1:length(g$mumm_libs)) { try({
      met_ids_out <- do.call("rbind",lapply(1:length(g$metab_data_index),function(i){
        if(g$omicsList[[ g$metab_data_index[i] ]][["dataType"]] == "met_combined") return(NULL);

        analysis_name <- paste0(g$mumm_libs[mumm_lib_i],"_",g$omicsList[[ g$metab_data_index[i] ]][["dataType"]],
                              "_", g$contrast_strings_file[ci])
        
        output_mummichog_reldir_tmp <- file.path(g$mumm_working_dir, analysis_name)
        output_mummichog_absdir_tmp <- file.path(working_dir, output_mummichog_reldir_tmp)
        met_ids <- read.delim(file.path(output_mummichog_absdir_tmp, "mummichog_matched_compound_all.csv"),
                              header=TRUE, stringsAsFactors = F, sep=",")
    
        met_ranked <- stats::na.omit(limma::topTable(g$omicsList[[ g$metab_data_index[i] ]][["fit"]], adjust="BH", n=Inf,
                                      sort.by='p', coef=g$loop_list[ci])[,c("mz", "P.Value", "logFC")] );
        met_ranked[,"rank"] <- ( -log(met_ranked[,"P.Value"]) * sign(met_ranked[,"logFC"]) )
      
        for(n in 1:nrow(met_ids)){
          met_ids[n,"rank"] <- met_ranked[which(met_ids[n,"Query.Mass"] == met_ranked[,"mz"])[1], "rank"]
        }

        met_ids <- met_ids[,c("Matched.Compound", "rank")]
        colnames(met_ids) <- c("GeneName", "rank")
        met_ids
      }))

      output_filename <- paste0("Met_Ranked_", g$mumm_libs[mumm_lib_i],"_", g$contrast_strings_file[ci], "_Combined.rnk")
 
      write.table(met_ids_out, file=file.path(g$output_mummichog_absdir, output_filename), sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE)
      rnk_files[[ci]] <- c(rnk_files[[ci]], file.path(output_mummichog_path, output_filename) )
      names(rnk_files[[ci]])[length(rnk_files[[ci]])] <- g$mumm_libs[mumm_lib_i]

      if(length(g$gene_data_index)==0){
        next;
      }
      
      if(length(g$gene_data_index)>1){
        gene_file <- "GSEA_combined_"
      } else {
        gene_file <- paste0("GSEA_",g$omicsList[[ g$gene_data_index[1] ]][["dataType"]],"_")
      }
      if(!grepl("Human", species)){ 
        gene_file <- paste0(gene_file,"Uppercase_")
      }

      fname_suffix <- paste0(g$contrast_strings_file[ci], ".rnk")
      gene_file <- file.path(g$output_contrast_path_files,paste0(gene_file,fname_suffix))
        
      gene_rnk <- read.delim(gene_file, stringsAsFactors = F, header=T)
      
      output_filename <- paste0("Integrated_Ranked_", g$mumm_libs[mumm_lib_i], "_", fname_suffix)
      out_path=file.path(g$mumm_working_dir, output_filename)
      write.table(rbind(gene_rnk,met_ids_out), file=out_path,
                  sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE)
      
      rnk_files[[ci]] <- c(rnk_files[[ci]], out_path)
      names(rnk_files[[c]])[length(rnk_files[[ci]])] <- g$mumm_libs[mumm_lib_i]
      
      # # Randomize for Control
      # met_ids_out["GeneName"] <- met_ids_out[sample(nrow(met_ids_out)),"GeneName"]
      # output_filename <- paste("MetRanked_Control_", mumm_libs[mumm_lib_i],"_",
      #                          gsub("-","_",g$contrast_strings[ci]), "_Combined.rnk", sep="" )
      # 
      # write.table(met_ids_out, file=file.path(output_mummichog_path, output_filename), sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE)
      # 
      # output_filename <- paste("Integrated_Ranked_Control_", mumm_libs[mumm_lib_i],"_",
      #                          gsub("-","_",g$contrast_strings[ci]), ".rnk", sep="" )
      # 
      # write.table(rbind(gene_rnk,met_ids_out), file=file.path(output_mummichog_path, output_filename),
      #             sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE)
      # 
      # rnk_files[[ci]] <- c(rnk_files[[ci]], file.path(output_mummichog_path, output_filename) )
      # names(rnk_files[[ci]])[length(rnk_files[[ci]])] <- mumm_libs[mumm_lib_i]
    }) }
  }

  g$combine_rnk_files = rnk_files
  g$calls = c(g$calls, "g.combine.rnk")
  g
}

#' MOMENTA
#'
#' Run MOMENTA multi-omic metabolic enrichment and network analysis
#' 
#' @param g omics notebook state, created by g.notebook.setup()
#' @param deps automatically run dependencies?
#'
#' @return updated omics notebook state
#' 
#' @export
g.momenta = function(g, working_dir="4_MOMENTA_Integrated", deps=T){
  if(deps) g=g.run.deps(g, c("g.combine.rnk"))
  gsea_working_path <- file.path(g$output_path, working_dir)
  if( dir.exists(gsea_working_path) == FALSE ) { dir.create(gsea_working_path) }

  geneset_lookup <-  read.delim(get.data.fileconn(file.path("feature","MatchedFeatureSets.txt")))

  gmt_files <- as.vector(geneset_lookup[which(geneset_lookup[,"Species"]==species),"GeneSet"])
  gmt_names <- as.vector(geneset_lookup[which(geneset_lookup[,"Species"]==species),"Name"])

  for(gmt_i in 1:length(gmt_files)){

    dest_gmt_file2 <- get.data.fileconn(file.path("feature",gmt_files[gmt_i]))

    analysis_names <- vector("list", length(g$loop_list))
    for (i in 1:length(g$loop_list)) {
      analysis_names[[i]] <- paste0("Combined_",gmt_names[gmt_i],"_",g$contrast_strings_file[i])
    }

    rnk_files_select <-vector("list", length(g$loop_list))
    for (i in 1:length(g$loop_list)){ 
      rnk_files_select[[i]] <- g$combine_rnk_files[[i]][stringr::str_detect(gmt_names[gmt_i], names(g$combine_rnk_files[[i]]))]
    }

    for (i in 1:length(analysis_names)){
      for(j in 1:length(rnk_files_select[[i]])){
        analysis_name <- paste(gsub(".*/","", gsub(".rnk", "", rnk_files_select[[i]][[j]])), "_", gmt_names[gmt_i], sep="" ) #analysis_names[[i]]
        rnk_file <- rnk_files_select[[i]][[j]]

        runGSEA(rnk=rnk_file, gmt=dest_gmt_file2, analysisName=analysis_name, out_path=gsea_working_path )
      }
    }
    close(dest_gmt_file2)
  }

  g$calls = c(g$calls, "g.momenta")
  g
}
