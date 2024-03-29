#-------------------------------------------------
#' Draw Venn Diagram
#'
#' This function generates venn diagrams/similar plots
#' 
#' @param item_list a list of item lists
#' @param item_name a name for the items
#' @param outputpath output file path for plots
#' 
#' @examples
#' 
#' @export
drawVenn <- function(item_list, item_name, outputpath=output_plots_path){
  if (length(item_list) > 5){
    item_list <- item_list[1:5]
  }
  fill_col <- grDevices::rainbow(length(item_list))
  futile.logger::flog.threshold(futile.logger::ERROR, name="VennDiagramLogger");
  venn <- VennDiagram::venn.diagram( x=item_list, filename=NULL,lty="blank",# height=2000, width=2000, 
                        cat.default.pos="outer", fill=fill_col, main=paste("Overlap: ", item_name, sep=""), 
                        fontfamily="sans", ext.text=FALSE, disable.logging=T);
  
  output_filename <- file.path(outputpath, paste("VennDiagram_", item_name,".pdf", sep="") );
  pdf(output_filename, width=4, height=3);
  print(UpSetR::upset(UpSetR::fromList(item_list), order.by = "freq", text.scale=1.2) );
  grid::grid.newpage();
  grid::pushViewport(grid::viewport(width=unit(0.8, "npc"), height=unit(0.8, "npc")));
  grid::grid.draw(venn);
  tmp<-dev.off();  
}

