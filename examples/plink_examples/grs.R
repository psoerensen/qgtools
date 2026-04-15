Glist <- readRDS(file = "C:/Users/au223366/Dropbox/Mohammed/Simulations/Glist_EUR_10k.RDS")


# Example inputs
file <- Glist$bedfiles        # bedfile
n <- Glist$n                  # number of individuals
m <- Glist$mchr               # number of SNPs
nt <- 3                       # number of traits

cls <- sample(1:m, m)         # SNP indices (1-based)
af  <- runif(m, 0.05, 0.5)    # allele frequencies

# SNP x trait matrix
S <- matrix(rnorm(m * nt), nrow = m, ncol = nt)

# Call your function
grs <- qgtools:::mtgrsbed_matrix(
  file = file,
  n = n,
  cls = cls,
  af = af,
  scale = TRUE,
  S = S,
  nthreads = 1
)
