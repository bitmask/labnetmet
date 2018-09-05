
library(devtools)
devtools::install_github("bitmask/labnetmet")
library(labnetmet)

# how different are the drosophilia networks from Froehlich2007?

L1 <- matrix( c(0, 1, 1, 1,  0, 0, 1, 0,  0, 1, 0, 0,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(L1) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

L2 <- matrix( c(0, 1, 1, 1,  0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(L2) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

L3 <- matrix( c(0, 0, 1, 1,  0, 0, 1, 0,  0, 1, 0, 0,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(L3) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

L4 <- matrix( c(0, 1, 0, 1,  0, 0, 1, 0,  0, 1, 0, 0,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(L4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

plot_dist(list(L1, L2, L3, L4), trans_dist, "~/grid-Dr.pdf")




# how different are the estrogen receptor networks from Froehlich2007?

#             1       2       3         4        5       6        7       8       9      10       11         12        13
sgenes <- c('ESR1', 'AKT2', 'HSPB8', 'CCNG2', 'FOXA1', 'STC2', 'XBP1', 'AKT1', 'DDR1', 'GPR30', 'BCL2', 'LOC120224', 'GDF15')

E1 <- matrix( c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 
                0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
               ), nrow=13, ncol=13, byrow=TRUE)
dimnames(E1) <- list(sgenes, sgenes)

E2 <- matrix( c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
               ), nrow=13, ncol=13, byrow=TRUE)
dimnames(E2) <- list(sgenes, sgenes)

E3 <- matrix( c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
               ), nrow=13, ncol=13, byrow=TRUE)
dimnames(E3) <- list(sgenes, sgenes)

E4 <- matrix( c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
               ), nrow=13, ncol=13, byrow=TRUE)
dimnames(E4) <- list(sgenes, sgenes)

plot_dist(list(E1, E2, E3, E4), trans_dist, "~/grid-E2.pdf")


# networks from the tests

a4 <- matrix( c(0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 0,  0, 1, 0, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(a4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

b4 <- matrix( c(0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 1, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(b4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

c4 <- matrix( c(0, 0, 0, 1,  0, 0, 0, 1,  0, 0, 0, 1,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(c4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

d4 <- matrix( c(0, 0, 0, 1,  0, 0, 0, 1,  1, 0, 0, 0,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(d4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

e4 <- matrix( c(0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(e4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

f4 <- matrix( c(0, 1, 1, 0,  0, 0, 1, 0,  0, 0, 0, 1,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(f4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

g4 <- matrix( c(0, 1, 1, 1,  0, 0, 1, 0,  0, 0, 0, 1,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(g4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

h4 <- matrix( c(0, 1, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(h4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

i4 <- matrix( c(0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(i4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

j4 <- matrix( c(0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(j4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

z4 <- matrix( c(0, 0, 0, 0,  1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(z4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))


plot_dist( list(a4, b4, c4, d4, e4, f4, g4, h4, i4, j4, z4) , trans_dist, "~/grid-tests.pdf")




