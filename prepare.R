filename <- c("CiP_count.csv","CiP_metadata.csv")
for(file in filename) {
    download.file(
        paste0("https://github.com/marthazak/CiP/blob/main/data/", file),
        file)
}

