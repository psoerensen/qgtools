set.seed(123)

# ------------------------------------------------------------
# Dimensions
# ------------------------------------------------------------
n_animals <- 300
n_dams    <- 120
years     <- 2015:2020

# ------------------------------------------------------------
# Pedigree-like structure
# ------------------------------------------------------------
Dam <- sample(1:n_dams, n_animals, replace = TRUE)

# Litter number per dam (1–3)
Litter <- ave(Dam, Dam, FUN = function(x) sample(1:3, length(x), replace = TRUE))

# Litter ID within dam (P in DMU)
P <- paste(Dam, Litter, sep = "_")

# ------------------------------------------------------------
# Fixed-effect covariates (INTEGER*4 style)
# ------------------------------------------------------------
Month  <- sample(1:12, n_animals, replace = TRUE)
Damage <- sample(2:7, n_animals, replace = TRUE)
Sex    <- sample(c(1, 2), n_animals, replace = TRUE)     # 1=male, 2=female
HY     <- sample(years, n_animals, replace = TRUE)
A      <- seq_len(n_animals)                             # Animal ID

# ------------------------------------------------------------
# Random effects (latent)
# ------------------------------------------------------------
# Direct genetic effect
A_ge <- rnorm(n_animals, sd = 2.0)

# Maternal genetic effect (shared within dam)
Dam_ge <- rnorm(n_dams, sd = 1.5)[Dam]

# Environmental litter-within-dam
L_Dam <- rnorm(length(unique(P)), sd = 1.0)
names(L_Dam) <- unique(P)
L_Dam_eff <- L_Dam[P]

# ------------------------------------------------------------
# Generate weights (REAL*4 style)
# ------------------------------------------------------------
# Birth weight
V0 <- 3.5 +
  0.1 * Month +
  0.2 * Damage +
  0.3 * (Sex == 1) +
  A_ge + Dam_ge + L_Dam_eff +
  rnorm(n_animals, sd = 1.0)

# Weight at 2 months
V1 <- V0 + rnorm(n_animals, mean = 8, sd = 1.5)

# Weight at 4 months
V2 <- V1 + rnorm(n_animals, mean = 10, sd = 2.0)

# ------------------------------------------------------------
# Gains (derived, as in DMU examples)
# ------------------------------------------------------------
G1 <- V1 - V0          # Gain 0–2 months
G2 <- V2 - V0          # Gain 0–4 months
G3 <- V2 - V1          # Gain 2–4 months

# ------------------------------------------------------------
# Final dataset (DMU-like)
# ------------------------------------------------------------
data <- data.frame(
  # INTEGER variables
  Month  = as.integer(Month),
  Damage = as.integer(Damage),
  Litter = as.integer(Litter),
  Sex    = as.integer(Sex),
  HY     = as.integer(HY),
  A      = as.integer(A),
  Dam    = as.integer(Dam),
  P      = P,

  # REAL variables
  V0 = V0,
  V1 = V1,
  V2 = V2,
  G1 = G1,
  G2 = G2,
  G3 = G3
)

# Inspect
str(data)
head(data)


formulas <- list(
  V0 = V0 ~ Month + Damage + Litter + Sex + HY +
    (1 | P) + (1 | Dam) + (1 | A),
  V1 = V1 ~ Month + Damage + Litter + Sex + HY +
    (1 | P) + (1 | Dam) + (1 | A)
)

qg_parse_formulas(formulas)
