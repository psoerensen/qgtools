# Mouse data files
mouse <- readRDS(url("https://github.com/psoerensen/bgcourse/raw/main/data/mouseqtl.rds"))
mouse$id <- rownames(mouse)

pedigree <- readRDS(url("https://github.com/psoerensen/bgcourse/raw/main/data/pedigree.rds"))
pedigree$date_of_birth <- 1:nrow(pedigree)

genotypes <- readRDS(url("https://github.com/psoerensen/bgcourse/raw/main/data/genotypes.rds"))
colnames(genotypes) <- paste0("M",1:ncol(genotypes))

genotypes <- cbind(id=rownames(genotypes),genotypes)

write.csv2(mouse[,c("id","sire", "dam", "sex", "reps", "Gl", "BW", "M227", "M1139")],
           file = "inst/examples/mouse/mouse.csv",
           row.names = FALSE)

write.csv2(pedigree[,c("id","sire","dam","date_of_birth")],
           file = "inst/examples/mouse/pedigree.csv",
           row.names = FALSE)

write.csv2(genotypes,
           file = "inst/examples/mouse/genotypes.csv",
           row.names = FALSE)


# Human plink data files
url <- "https://github.com/psoerensen/qgdata/raw/main/simulated_human_data/human.bed"
download.file(url = url, mode = "wb", destfile = "inst/examples/human/human.bed")
url <- "https://github.com/psoerensen/qgdata/raw/main/simulated_human_data/human.bim"
download.file(url = url, destfile = "inst/examples/human/human.bim")
url <- "https://github.com/psoerensen/qgdata/raw/main/simulated_human_data/human.fam"
download.file(url = url, destfile = "inst/examples/human/human.fam")

url <- "https://github.com/psoerensen/qgdata/raw/main/simulated_human_data/human.pheno"
download.file(url = url, destfile = "inst/examples/human/human.pheno")
url <- "https://github.com/psoerensen/qgdata/raw/main/simulated_human_data/human.covar"
download.file(url = url, destfile = "inst/examples/human/human.covar")

