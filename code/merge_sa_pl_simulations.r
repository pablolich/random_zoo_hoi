#first load all files, in whatever format, and change them to be saved only with
#space, like currently.

for (n in c(1,2,3,4,5,6,7,8)){
  for (d in c(2,3,4,5,6,7)){
    f1 <-  paste0("n_", n, "_d_", d, ".csv")
    path = "../data/merged_files/sa_data2/"
    if (file.exists(paste0(path,f1))){
      dt <- read.table(paste0(path,f1), sep = ",", header = T) #load r
      colnames(dt) = c("n", "d", "nsol", "npos")
      #save in correct format
      write.table(dt, paste0("../data/merged_files/sa_data_mod/",f1), sep = " ", row.names = F)
    }
  }
}