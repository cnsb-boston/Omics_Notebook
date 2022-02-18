args <- commandArgs(trailingOnly=TRUE)

# https://github.com/molgenis/molgenis-pipelines/wiki/How-to-source-another_file.R-from-within-your-R-script
LocationOfThisScript = function() {
  this.file = NULL
  # This file may be 'sourced'
  for (i in -(1:sys.nframe())) {
    if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
  }

  if (!is.null(this.file)) return(dirname(this.file))

  # But it may also be called from the command line
  cmd.args = commandArgs(trailingOnly = FALSE)
  cmd.args.trailing = commandArgs(trailingOnly = TRUE)
  cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
  res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)

  # If multiple --file arguments are given, R uses the last one
  res = tail(res[res != ""], 1)
  if (0 < length(res)) return(dirname(res))

  # Both are not the case. Maybe we are in an R GUI?
  return(NULL)
}
current.dir = LocationOfThisScript()

notebook_path=if(is.null(current.dir)) "/project/cnsbomic/Omics_Notebook" else current.dir
source(paste0(notebook_path,"Config.R"))
source(paste0(notebook_path,"NotebookGUI.R"))

native.cmd = function(notebook="/home", analysis="/data", parameters="Parameters.R"){
  c("Rscript", paste0(notebook,"/src/Pipeline.R"), notebook, analysis, parameters)
}

run.native = function(args){
  cmdargs=list(notebook=args[1],analysis=args[2])
  if(length(args)>2) cmdargs$parameters=args[3]
  command = do.call("native.cmd",cmdargs)
  system2(command[1], command[-1], env=character(R_LIBS="/usr/local/R/local-library"))
}

run.generic.container = function(cmd, image, bindopt, extras=list(pre="",post=""), parameters){
  command = c(cmd, extras$pre,
              bindopt, paste0(notebook_path, ":/home:rw"),
              bindopt, paste0(analysis_dir, ":/data:rw"),
              bindopt, paste0(libdir, ":/usr/local/lib/R/local-library"), extras$post,
              image, native.cmd(parameters=parameters))
  command=command[command!=""]
  system2(command[1], command[-1])
}

run.docker = function(parameters,image){
  run.generic.container(cmd=c("docker","run","-it","--rm","-u","docker"), image="cnsbboston/omicsnotebook", bindopt="-v", parameters=parameters)
}

run.singularity = function(parameters,image){
  run.generic.container(cmd=c("singularity","run"), image=image, bindopt="--bind", parameters=parameters)
}

getopt = function(args, choices){
  ret=as.list(choices[,3])
  names(ret)=choices[,1]
 
  i=1
  while(i<=length(args)){
    if(substr(args[i],1,1)=="-"){
      opt=sub("^-*","",args[i])
      r_i=choices[,1]==opt | choices[,2]==opt
      if(any(r_i)){
        r_i=which(r_i)[1]
        if(i==length(args) || substr(args[i+1],1,1)=="-"){
          val=TRUE
        } else {
          val=args[i+1]
          i=i+1
        }
        ret[[r_i]]=val
      }
    }
    i=i+1
  }

  for(i in 1:nrow(choices)){
    class(ret[[i]])=choices[i,4]
  }
  ret
}

choices=matrix(c(
  'param', 'p', 'Parameters.R', 'character',
  'container', 'c', singularity_img, 'character',
  'nogui', 'g', F, 'logical',
  'help', 'h', F, 'logical'
  ), byrow=T, ncol=4)
opts=getopt(args,choices)

if(opts$help){
  for(i in 1:nrow(choices)){
    print(paste0("--",choices[i,1]," -",choices[i,2]," (default: ",choices[i,3],")"))
  }
  quit('no')
}

if(opts$nogui){
  analysis_dir = getwd()
} else {
  analysis_dir = make.gui(startdir=startdir)
}

parameters = if(length(args)>1) args[length(args)] else "Parameters.R"
switch(args[1],
       "Docker"=run.docker(opts$param,opts$container),
       "Singularity"=run.singularity(opts$param,opts$container),
       "GUI"=write(analysis_dir,stdout()),
       run.native(args))
