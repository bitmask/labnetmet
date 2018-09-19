
context("Simple sanity checks for distance measures")

test_that("dist to self is 0", {
    a2 <- matrix( c(0, 1, 0, 0), nrow=2, ncol=2, byrow=TRUE)
    dimnames(a2) <- list(c("A", "B"), c("A", "B"))

    expect_equal(hamming_dist(a2, a2), 0)
    expect_equal(subgraph_dist(a2, a2), 0)
    expect_equal(trans_dist(a2, a2), 0)

})


test_that("add inverse edge is difference", {
    # 2d examples
    a2 <- matrix( c(0, 1, 0, 0), nrow=2, ncol=2, byrow=TRUE)
    dimnames(a2) <- list(c("A", "B"), c("A", "B"))

    b2 <- matrix( c(0, 1, 1, 0), nrow=2, ncol=2, byrow=TRUE) # adds the inverse edge
    dimnames(b2) <- list(c("A", "B"), c("A", "B"))

    expect_equal(hamming_dist(a2, b2), 1)
    expect_equal(hamming_dist(b2, a2), 1)

    expect_equal(subgraph_dist(a2, b2), 1)
    expect_equal(subgraph_dist(b2, a2), 1)

    expect_equal(trans_dist(a2, b2), 1)
    expect_equal(trans_dist(b2, a2), 1)
})

test_that("3d tests", {

    # 3d examples

    # add one transitive edge
    a3 <- matrix( c(0, 1, 0,  0, 0, 1,  0, 0, 0), nrow=3, ncol=3, byrow=TRUE)
    dimnames(a3) <- list(c("A", "B", "C"), c("A", "B", "C"))

    b3 <- matrix( c(0, 1, 1,  0, 0, 1,  0, 0, 0), nrow=3, ncol=3, byrow=TRUE)
    dimnames(b3) <- list(c("A", "B", "C"), c("A", "B", "C"))

    expect_equal(hamming_dist(a3, b3), 1)
    expect_equal(hamming_dist(b3, a3), 1)

    expect_equal(subgraph_dist(a3, b3), 1/3)
    expect_equal(subgraph_dist(b3, a3), 1/3)

    expect_equal(trans_dist(a3, b3), 0)
    expect_equal(trans_dist(b3, a3), 0)

    # test same matricies but with mixed up rows
    aa3 <- matrix( c(0, 0, 1,  0, 1, 0,  0, 0, 0), nrow=3, ncol=3, byrow=TRUE)
    dimnames(aa3) <- list(c("B", "A", "C"), c("A", "B", "C"))

    expect_equal(hamming_dist(a3, aa3), 0)
    expect_equal(hamming_dist(aa3, a3), 0)

    expect_equal(subgraph_dist(a3, aa3), 0)
    expect_equal(subgraph_dist(aa3, a3), 0)

    expect_equal(trans_dist(a3, aa3), 0)
    expect_equal(trans_dist(aa3, a3), 0)

    # add the inverse edges
    ai3 <- matrix( c(0, 1, 0,  1, 0, 1,  0, 1, 0), nrow=3, ncol=3, byrow=TRUE)
    dimnames(ai3) <- list(c("A", "B", "C"), c("A", "B", "C"))

    expect_equal(hamming_dist(a3, ai3), 2)
    expect_equal(hamming_dist(ai3, a3), 2)

    expect_equal(subgraph_dist(a3, ai3), 2/3)
    expect_equal(subgraph_dist(ai3, a3), 2/3)

    expect_equal(trans_dist(a3, ai3), 3)
    expect_equal(trans_dist(ai3, a3), 3)

    # triangle inequality
    c3 <- matrix( c(0, 1, 1,  0, 0, 0,  0, 0, 0), nrow=3, ncol=3, byrow=TRUE)
    dimnames(c3) <- list(c("A", "B", "C"), c("A", "B", "C"))
    expect_true(trans_dist(a3, b3) + trans_dist(b3, c3) >= trans_dist(a3,c3))
    expect_true(trans_dist(a3, c3) + trans_dist(c3, b3) >= trans_dist(a3,b3))
    expect_true(trans_dist(c3, a3) + trans_dist(a3, b3) >= trans_dist(c3,b3))
})

test_that("4d tests", {

    # 4d examples

    # connects one level farther down the cascade  -- currently counted as 2 edits
    a4 <- matrix( c(0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 0,  0, 1, 0, 0), nrow=4, ncol=4, byrow=TRUE)
    dimnames(a4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

    b4 <- matrix( c(0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 1, 0), nrow=4, ncol=4, byrow=TRUE)
    dimnames(b4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

    expect_equal(hamming_dist(a4, b4), 2)
    expect_equal(hamming_dist(b4, a4), 2)

    expect_equal(subgraph_dist(a4, b4), 1/2)
    expect_equal(subgraph_dist(b4, a4), 1/2)

    expect_equal(trans_dist(a4, b4), 1)
    expect_equal(trans_dist(b4, a4), 1)

    # connects one layer higher
    c4 <- matrix( c(0, 0, 0, 1,  0, 0, 0, 1,  0, 0, 0, 1,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
    dimnames(c4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

    d4 <- matrix( c(0, 0, 0, 1,  0, 0, 0, 1,  1, 0, 0, 0,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
    dimnames(d4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

    expect_equal(hamming_dist(c4, d4), 2)
    expect_equal(hamming_dist(d4, c4), 2)

    expect_equal(subgraph_dist(c4, d4), 1/2)
    expect_equal(subgraph_dist(d4, c4), 1/2)

    expect_equal(trans_dist(c4, d4), 1)
    expect_equal(trans_dist(d4, c4), 1)

    # add two edges
    # triangle inequality for non-transitive equivalent measures
    e4 <- matrix( c(0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
    dimnames(e4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

    f4 <- matrix( c(0, 1, 1, 0,  0, 0, 1, 0,  0, 0, 0, 1,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
    dimnames(f4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

    expect_equal(hamming_dist(e4, f4), 1)

    expect_equal(subgraph_dist(e4, f4), 0.3)

    expect_equal(trans_dist(e4, f4), 0)

    g4 <- matrix( c(0, 1, 1, 1,  0, 0, 1, 0,  0, 0, 0, 1,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
    dimnames(g4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

    expect_equal(hamming_dist(e4, g4), 2)
    expect_equal(subgraph_dist(e4, g4), 1/2)
    expect_equal(trans_dist(e4, g4), 0)

    near1 <- subgraph_dist(e4,f4)
    near2 <- subgraph_dist(f4,g4)
    far <- subgraph_dist(e4,g4)
    expect_true(near1 + near2 >= far)
})

test_that("triangle inequality holds", {

    # triangle inequality for transitive equivalent measures
    h4 <- matrix( c(0, 1, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
    dimnames(h4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

    i4 <- matrix( c(0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
    dimnames(i4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

    j4 <- matrix( c(0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
    dimnames(j4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

    near1 <- trans_dist(h4,i4)
    near2 <- trans_dist(i4,j4)
    far <- trans_dist(h4,j4)
    expect_true(near1 + near2 >= far)

    # bigger example
    x5 <- matrix( c(0, 1, 0, 0, 0,  0, 0, 1, 0, 0,  0, 0, 0, 1, 0,  0, 0, 0, 0, 1,  0, 0, 0, 0, 0), nrow=5, ncol=5, byrow=TRUE)
    dimnames(x5) <- list(c("A", "B", "C", "D", "E"), c("A", "B", "C", "D", "E"))

    y5 <- matrix( c(0, 0, 0, 0, 1,  0, 0, 0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0, 0), nrow=5, ncol=5, byrow=TRUE)
    dimnames(y5) <- list(c("A", "B", "C", "D", "E"), c("A", "B", "C", "D", "E"))

    z5 <- matrix( c(0, 0, 0, 0, 0,  0, 0, 1, 0, 0,  0, 0, 0, 1, 0,  0, 0, 0, 0, 1,  0, 0, 0, 0, 0), nrow=5, ncol=5, byrow=TRUE)
    dimnames(z5) <- list(c("A", "B", "C", "D", "E"), c("A", "B", "C", "D", "E"))

    expect_true(trans_dist(x5, y5) + trans_dist(y5, z5) >= trans_dist(x5,z5))
    expect_true(trans_dist(x5, z5) + trans_dist(z5, y5) >= trans_dist(x5,y5))
    expect_true(trans_dist(z5, x5) + trans_dist(x5, y5) >= trans_dist(z5,y5))

})

test_that("inverse cascade is maximally different", {

    # inverse cascade is maximally different (no self edges)
    e4 <- matrix( c(0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
    dimnames(e4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))
    z4 <- matrix( c(0, 0, 0, 0,  1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0), nrow=4, ncol=4, byrow=TRUE)
    dimnames(z4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

    expect_equal(trans_dist(e4, z4), 12)
    expect_equal(trans_dist(z4, e4), 12)
})
