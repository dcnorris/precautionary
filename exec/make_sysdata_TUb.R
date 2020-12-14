## Make the internal data {A, U, b} required for exact simulation of 3+3 design

T <- U <- b <- list()
for(D in 2:8){
  T_raw <- read.table(paste0("T", D, ".tab"))
  J <- nrow(T_raw)/2
  T[[D]] <- aperm(array(as.matrix(T_raw), dim = c(2, J, D)
                        , dimnames = list(c = 1:2, j = NULL, d = paste0("D",1:D)))
                  , perm = c("c","d","j"))
  U[[D]] <- cbind(apply(T[[D]], c("j","d"), sum, na.rm=TRUE),
                  apply(3 - T[[D]], c("j","d"), sum, na.rm=TRUE))
  dimnames(U[[D]]) <- list(j = NULL, pq = c(paste0("D", rep(1:D, 2))))
  b[[D]] <- apply(log(choose(3, T[[D]])), "j", sum, na.rm=TRUE)
}

usethis::use_data(T, U, b,
                  internal = TRUE,
                  overwrite = TRUE)
