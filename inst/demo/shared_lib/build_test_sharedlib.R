demo_dir <- system.file("demo/shared_lib", package = "qgtools")
test_dir <- "C:\Users\au223366\Documents\GitHub\qgtools\inst\demo/shared_lib"
setwd(demo_dir)
list.files(demo_dir)

system("make")
system("Rscript test.R")
libpath <- system.file("demo/shared_lib/libgrsbed_cpp.so", package = "qgtools")
dyn.load(libpath)

# Load shared libraries
dyn.load("libgrsbed_cpp.so")
dyn.load("libgrsbed_fortran.so")

is.loaded("mtgrsbed_core_c")
is.loaded("mtgrsbed_core_f")

n  <- 1000
m  <- 5000
nt <- 3

b_snp <- rnorm(m * nt)

grs_cpp <- double(n * nt)
grs_f   <- double(n * nt)

dummy_file <- ""   # not used in demo
cls <- as.integer(rep(1, m))
af  <- rep(0.2, m)

cat("Running C++...\n")
system.time({
  .C("mtgrsbed_core_c",
     as.character(dummy_file),
     as.integer(n),
     as.integer(cls),
     as.double(af),
     as.integer(m),
     as.integer(nt),
     as.integer(1),
     as.double(b_snp),
     grs_cpp,
     as.integer(4),
     as.integer(64),
     as.integer(1024),
     as.integer(32)
  )
})

cat("Running Fortran...\n")
system.time({
  .C("mtgrsbed_core_f",
     as.character(dummy_file),
     as.integer(n),
     as.integer(cls),
     as.double(af),
     as.integer(m),
     as.integer(nt),
     as.logical(TRUE),
     as.double(b_snp),
     grs_f,
     as.integer(4),
     as.integer(64),
     as.integer(1024),
     as.integer(32)
  )
})

cat("Equal:", all.equal(grs_cpp, grs_f), "\n")
