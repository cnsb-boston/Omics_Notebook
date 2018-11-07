
# Import annotation file
annot <- data.frame(openxlsx::read.xlsx(annotation_filename, 1, colNames=FALSE)); # read annotation
contrastgroups <- gsub("\\.","", make.names(na.omit(t(annot[1,-1:-3]))));
annot<-annot[c(-1,-4),];

# Make omicsList object
omicsList <- vector("list", (nrow(annot)-2) )
for (i in 3:nrow(annot)){
  omicsList[[i-2]] <- vector("list", 14);
  omicsList[[i-2]][1] <- make.names(annot[i,1]);
  names(omicsList[[i-2]])[1] <- gsub("\\.","", make.names("dataType"));
  omicsList[[i-2]][2] <- make.names(annot[i,2]);
  names(omicsList[[i-2]])[2] <- "dataFormat";
  omicsList[[i-2]][[3]] <- annot[i,3];
  names(omicsList[[i-2]])[3] <- "filename";
}

# Format annotation table
annot[3:nrow(annot),3] <- annot[3:nrow(annot),1];
annot <- t(annot[,-1:-2]);
rownames(annot) <- NULL;
colnames(annot)<-make.unique(annot[1,]);
annot<-data.frame((annot[-1,]));
annot[,"Group"]<-gsub("\\.","",make.names(annot[,"Group"]))
annot[,"SampleName"]<-gsub("\\.","",make.names(make.unique(as.character(annot[,"SampleName"]))))

# Get batch information if there
try({
  annot2 <- openxlsx::read.xlsx(annotation_filename, 2, colNames=FALSE, rowNames=TRUE)
  annot2 <- data.frame(t(annot2))
  rownames(annot2) <- NULL;
  annot <- data.frame(cbind(annot, annot2))
}, silent=TRUE)

# Changes for generating txt report
if(txtFolder==TRUE) {
  for (i in 1:length(omicsList)){
    omicsList[[i]][[3]] <- file.path(working_dir, "txt", omicsList[[i]][[3]])
  }
} else {
  for (i in 1:length(omicsList)){
    omicsList[[i]][[3]] <- file.path(working_dir, omicsList[[i]][[3]])
  }
}

# Import data and make eset objects
if(!newcontrastonly){
  for (i in 1:length(omicsList)){
    omicsList[[i]][[4]] <- data.frame(read.delim(omicsList[[i]][[3]], header=TRUE, stringsAsFactors=FALSE, check.names=FALSE));
    names(omicsList[[i]])[4] <- "RawData";
    
    omicsList[[i]][[5]] <- makeEset(data=omicsList[[i]][[4]], annotate=annot, type=omicsList[[i]][[1]], 
                                    data_format=omicsList[[i]][[2]], cutoff=0, uniprot_annotation=FALSE)
    names(omicsList[[i]])[5] <- "eSet";
    
    print(paste(omicsList[[i]][[1]], ": ", nrow(omicsList[[i]][[4]]), " features detected in raw data.", sep=""))
    print(paste(omicsList[[i]][[1]], ": ", nrow(omicsList[[i]][[5]]), " features detected in filtered data.", sep=""))
  }
}



# Make list to capture normalization paramters
normParam <- vector("list", length(omicsList) )
for (i in 1:length(normParam)) {
  normParam[[i]] <- vector("list", 2)
  names(normParam[[i]]) <- c("norm_method", "zero_percent")
}




library(shiny)

ui <- fluidPage(
  # App Title
  titlePanel("Normalization method"), 
  sidebarLayout(
    sidebarPanel(
      sliderInput(inputId = "Zero Filter", min=0,max=1, value=0.7)))
  
  
  
  
  )