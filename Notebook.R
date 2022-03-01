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
script.dir = LocationOfThisScript()

notebook_path=if(is.null(script.dir)) "/project/cnsbomic/Omics_Notebook" else script.dir
source(paste0(notebook_path,"/Config.R"))

native.cmd = function(notebook="/home", analysis="/data", parameters="Parameters.R"){
  c("Rscript", paste0(notebook,"/src/Pipeline.R"), notebook, analysis, parameters)
}

run.native = function(args){
  cmdargs=list(notebook=args[1],analysis=args[2])
  if(length(args)>2) cmdargs$parameters=args[3]
  command = do.call("native.cmd",cmdargs)
  #write(paste0(command,collapse=" "),stderr())
  system2(command[1], command[-1], env=paste0(env_vars,collapse=";"))
}

run.generic.container = function(cmd, image, bindopt, extras=list(pre="",post=""), parameters, libdir, outputlog=F){
  native=native.cmd(parameters=parameters)
  if(outputlog){
    logfile=paste0("/data/",format(Sys.time(), "%Y-%m-%d.%H:%M:%S"),".log.txt")
    native= paste0(c(native, "2>&1", "|", "tee", logfile), collapse=" ")
    native=c("sh","-c",paste0("'",native,"'"))
  }
  command = c(cmd, extras$pre, paste0("--env=",env_vars),
              bindopt, paste0("'", notebook_path, ":/home:rw'"),
              bindopt, paste0("'", analysis_dir, ":/data:rw'"),
              bindopt, paste0("'", libdir, ":/usr/local/lib/R/local-library'"),
              extras$post,
              image, native)
  command=command[command!=""]
  #write(paste0(command,collapse=" "),stderr())
  system2(command[1], command[-1])
}

run.docker = function(param,container,libdir,outputlog,...){
  if(length(container)!=1 || nchar(container)<1) container=docker_img
  run.generic.container(cmd=c("docker","run","-it","--rm","-u","docker"), image=container, bindopt="-v",
                        extras=list(pre="--workdir=/data",post=""),parameters=param, libdir=libdir, outputlog=outputlog)
}

run.singularity = function(param,container,libdir,outputlog,...){
  if(length(container)!=1 || nchar(container)<1) container=singularity_img
  run.generic.container(cmd=c("singularity","run"), image=container, bindopt="--bind",
                        extras=list(pre="--pwd=/data",post=""), parameters=param, libdir=libdir, outputlog=outputlog)
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
  'libdir', 'l', libdir, 'character',
  'container', 'c', "", 'character',
  'nogui', 'g', F, 'logical',
  'outputlog', 'o', F, 'logical',
  'help', 'h', F, 'logical'
  ), byrow=T, ncol=4)
opts=getopt(args,choices)

if(opts$help){
  print("[GUI | [[R] [Docker | Singularity]]] {options ...}")
  print("options:")
  for(i in 1:nrow(choices)){
    print(paste0("--",choices[i,1]," -",choices[i,2]," (default: ",choices[i,3],")"))
  }
  quit('no')
}

if(opts$nogui){
  analysis_dir = getwd()
} else {
  source(paste0(notebook_path,"/NotebookGUI.R"))
  analysis_dir = make.gui(startdir=startdir)
}

#short circuits
switch(args[1],
       "GUI"={write(analysis_dir,stdout()); quit("no");},
       "R"={args=args[-1]; native.cmd=function(...)c("R","--quiet");})

#normal analysis
switch(args[1],
       "Docker"=do.call(run.docker, opts),
       "Singularity"=do.call(run.singularity, opts),
       run.native(args))
