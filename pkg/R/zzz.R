".onLoad" <-function(lib, pkg)
{
  library.dynam("mixvarselect", package = pkg, lib.loc = lib)
  return(invisible(0)) 
}
