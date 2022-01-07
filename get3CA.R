library(SummarizedExperiment)
library(tools)

input_dir <- "/pfs/download3CAData/"
out_dir <- "/pfs/out/"
# input_dir <- "~/Documents/pfs/download3CAData/"
# out_dir <- "~/Documents/pfs/get3CA/"

data <- read.csv("https://github.com/BHKLAB-Pachyderm/download3CAData/raw/main/download_links.csv")

getSummarizedExp <- function(data_dir, rds_name, files, cells){
  print(data_dir)
  if (file_ext(files[1]) == "mtx") {
    matrix <- Matrix::readMM(paste0(data_dir, "/", files[1]))
  }
  if (file_ext(files[1]) == "rds") {
    matrix <- readRDS(paste0(data_dir, "/", files[1]))
  }
  cols <- suppressWarnings(data.table::fread(paste0(data_dir, "/", cells[1])))
  if ("cell_name" %in% names(cols)) {
    return(SummarizedExperiment(assays = matrix,colData = cols$cell_name))
  }
  if ("Cell_names" %in% names(cols)) {
    return(SummarizedExperiment(assays = matrix,colData = cols$Cell_names))
  }
}

filenames <- c()

for(disease_name in data$disease) {
  print(disease_name)
  dir_name <- paste0(out_dir, disease_name)
  dir.create(dir_name)
  suppressWarnings(unzip(paste0(input_dir, disease_name, ".zip"), exdir=dir_name))
  
  dirs <- list.dirs(dir_name)
  if (length(dirs) > 1) {
    for (data_dir in dirs) {
      files <- list.files(path = data_dir, pattern = "exp", ignore.case = TRUE)
      cells <- list.files(path = data_dir, pattern = "cell", ignore.case = TRUE)
      if (length(files) >= 1) {
        names <- strsplit(data_dir,"/")
        rds_name <- paste0(disease_name, "_", names[[1]][length(names[[1]])])
        saveRDS(getSummarizedExp(data_dir, rds_name, files, cells), paste0(out_dir, "3CA_", rds_name, ".rds"))
        filenames <- c(filenames, paste0("3CA_", rds_name, ".rds"))
      }
    }
  }
  else {
    files <- list.files(path = dir_name, pattern = "exp", ignore.case = TRUE)
    cells <- list.files(path = dir_name, pattern = "cell", ignore.case = TRUE)
    saveRDS(getSummarizedExp(dir_name, disease_name, files, cells), paste0(out_dir, "3CA_", disease_name, ".rds"))
    filenames <- c(filenames, paste0("3CA_", disease_name, ".rds"))
  }
  unlink(dir_name, recursive=TRUE)
}

data_list <- data.frame(matrix(data=NA, ncol=0, nrow=length(filenames)))
data_list$filename <- filenames
data_list$doi <- NA
data_list$download_link <- NA
write.csv(data_list, paste0(out_dir, "data_list.csv"))

