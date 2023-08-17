library(RColorBrewer)

nsim = 1000
nsppmax = 15
dev.off()

for (nspp in seq(2,nsppmax)){
  #initialize matrix A
  Arand = matrix(0, nrow = nspp, ncol = nspp)
  Afeas = matrix(0, nrow = nspp, ncol = nspp)
  Anofeas = matrix(0, nrow = nspp, ncol = nspp)
  #initialize hs
  hrand = matrix(0, nrow = 1, ncol = nspp)
  hfeas = matrix(0, nrow = 1, ncol = nspp)
  hnofeas = matrix(0, nrow = 1, ncol = nspp)
  nfeas = 0
  nnofeas = 0
  for (sim in seq(0,nsim-1)){
    filename_A = paste("A_nspp_", nspp, "sim", sim,".csv", sep = "") #change to _sim_
    path_feas_A = paste("../data/feas/", filename_A, sep = "")
    path_nofeas_A = paste("../data/no_feas/", filename_A, sep = "") 
    filename_h = paste("h_nspp_", nspp, "sim", sim,".csv", sep = "")
    path_feas_h = paste("../data/feas/", filename_h, sep = "")
    path_nofeas_h = paste("../data/no_feas/", filename_h, sep = "") 
    if (file.exists(path_feas_A)){
      #load file
      A = read.csv(path_feas_A, header = F)
      Afeas = Afeas + A
      Arand = Arand + A
      h = read.csv(path_feas_h, header = F)
      hfeas = hfeas + h
      hrand = hrand + h
      nfeas = nfeas + 1
      
    }
    else if(file.exists(path_nofeas_A)){
      A = read.csv(path_nofeas_A, header = F)
      Anofeas = Anofeas + A
      Arand = Arand + A
      h = read.csv(path_nofeas_h, header = F)
      hnofeas = hnofeas + h
      hrand = hrand + h
      nnofeas = nnofeas + 1
    }
    else{
      next
    }
  }
  #average out
  #avoid dividing by 0
  if (nfeas == 0){
    Afeasmean = Afeas
    Anofeasmean = Anofeas/nnofeas
    Arandmean = Arand/nsim
    hfeasmean = hfeas
    hnofeasmean = hnofeas/nnofeas
    hrandmean = hrand/nsim
  }
  else if(nnofeas == 0){
    Afeasmean = Afeas/nfeas
    Anofeasmean = Anofeas
    Arandmean = Arand/nsim
    hfeasmean = hfeas/nfeas
    hnofeasmean = hnofeas
    hrandmean = hrand/nsim
  }
  else{
    Afeasmean = as.matrix(Afeas/nfeas)
    Anofeasmean = as.matrix(Anofeas/nnofeas)
    Arandmean = as.matrix(Arand/nsim)
    hfeasmean = as.matrix(hfeas/nfeas)
    hnofeasmean = as.matrix(hnofeas/nnofeas)
    hrandmean = as.matrix(hrand/nsim)
  }
  if (all(Afeasmean == 0) | all(Anofeasmean == 0)){
    next
  }
  else{
    par(mfrow = c(1, 3))
    image(Afeasmean, main = "Feasible", 
          col= colorRampPalette(brewer.pal(8, "Blues"))(25),
          asp = 1)
    image(Arandmean, main = "Random", 
          col= colorRampPalette(brewer.pal(8, "Blues"))(25),
          asp = 1)
    image(Anofeasmean,main = "Not feasible", 
            col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
          asp = 1)

    par(mfrow = c(3,1))
    image(hfeasmean, main = "Feasible",
          col= colorRampPalette(brewer.pal(8, "Blues"))(25))
    image(hrandmean, main = "Random",
          col= colorRampPalette(brewer.pal(8, "Blues"))(25))
    image(hnofeasmean, main = "Not feasible",
          col= colorRampPalette(brewer.pal(8, "Blues"))(25))
  }
}