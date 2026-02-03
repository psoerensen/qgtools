# ------------------------------------------------------------------------------
# Variable names used in formulas refer to columns in the data frame / file
# ------------------------------------------------------------------------------

# Define formulas for each trait
formulas <- list(
  V0 = V0 ~ Month + DamAge + Litter + Sex + HY +
    (1 | L_Dam) +
    (1 | Dam) +
    (1 | Animal),

  V1 = V1 ~ Month + DamAge + Litter + Sex + HY +
    (1 | L_Dam) +
    (1 | Dam) +
    (1 | Animal)
)




# ------------------------------------------------------------------------------
# Define kernels (data sources)
# ------------------------------------------------------------------------------

pedigree <- makePEDlist(fnPED = "pedigree.txt")
grm      <- makeGRMlist(fnGRM = "grm.txt")


# ------------------------------------------------------------------------------
# Define variance components
# ------------------------------------------------------------------------------

vcs <- varcomp(

  litter_env = vc(
    index  = "L_Dam",
    traits = c("V0", "V1")
  ),

  maternal_genetic = vc(
    index  = "Dam",
    traits = c("V0", "V1"),
    kernel = ped_kernel(pedigree),
    prior  = prior_iw(df = 4, S = diag(2))
  ),

  direct_genetic = vc(
    index  = "Animal",
    traits = c("V0", "V1"),
    kernel = ped_kernel(pedigree),
    prior  = prior_iw(
      df = 4,
      S  = matrix(c(4, 0.8,
                    0.8, 3), 2, 2, byrow = TRUE)
    )
  ),
  residual = vc(
    index  = "Residual",
    traits = c("V0", "V1"),
    prior  = prior_diag(c(5, 4))
  )
)

qg_parse_formulas(formulas)
qg_parse_kernels(formulas)
qg_parse_variance_components(vcs)



obj_written <- write_qg_yaml(
  formulas = formulas,
  vcov     = vcov,
  Klist    = Klist,
  file     = tmp
)

obj_read <- read_qg_yaml(tmp)
