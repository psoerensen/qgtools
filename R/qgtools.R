# .onLoad <- function(libname, pkgname) {
#   lib <- switch(Sys.info()[["sysname"]],
#                 "Linux"   = "linux/libqgtools.so",
#                 "Darwin"  = "mac/libqgtools.dylib",
#                 "Windows" = "windows/libqgtools.dll"
#   )
#   dyn.load(system.file("libs", lib, package = pkgname))
# }
