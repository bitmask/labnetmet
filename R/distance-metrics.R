
# prevent devtools::check() from throwing warnings about these names that are not even global variables
utils::globalVariables(c("variable", "distance", "graph_name", "graph_name_y"))

#' Find the Hamming distance between two graphs as sum of the absolute value differences between adjacency matricies.
#'
#' @param x An adjacency matrix
#' @param y An adjacency matrix
#' @return Number representing the distance between x and y
#' @export
hamming_dist <- function(x,y) {
    stopifnot(all.equal(sort(rownames(x)), sort(colnames(x))))
    stopifnot(all.equal(sort(rownames(y)), sort(colnames(y))))
    stopifnot(all.equal(nrow(x), nrow(y)))

    # first need to match the row, cols
    z <- y[rownames(x), colnames(x)]

    # count the differences 
    differences <- sum(abs(x - z))
    # TODO is not normalized yet
    return(differences)
}


#' Find the subgraph distance between two graphs as the number of subgraphs that differ between the networks.
#'
#' @param x An adjacency matrix
#' @param y An adjacency matrix
#' @return Number representing the distance between x and y
#' @export
subgraph_dist <- function(x,y) {
    # x and y are full matricies

    # find all subgraphs of x and y, then compare them

    stopifnot(all.equal(sort(rownames(x)), sort(colnames(x))))
    stopifnot(all.equal(sort(rownames(y)), sort(colnames(y))))
    stopifnot(all.equal(nrow(x), nrow(y)))

    min_subgraph_size <- 2
    max_subgraph_size <- max(nrow(y) - 1, min_subgraph_size)
    subgraph_sizes <- min_subgraph_size:max_subgraph_size

    subgraphs <- unlist(lapply(subgraph_sizes, function(size) {utils::combn(rownames(y), size, simplify=F)}), recursive=FALSE)

    ix <- igraph::graph_from_adjacency_matrix(x)
    iy <- igraph::graph_from_adjacency_matrix(y[rownames(x), colnames(x)])

    isomorph <- unlist(lapply(subgraphs, function(subgraph) {igraph::isomorphic(igraph::induced_subgraph(ix, subgraph), igraph::induced_subgraph(iy, subgraph))}))

    # normalize result to length of number of subgraphs
    return ((length(isomorph) - sum(isomorph)) / length(isomorph))

}

#' XOR for finiteness
#'
#' @param e A distance
#' @param f A distance
#' @return True if both are either finite or infinte, False if they are not
same_reachability <- function(e,f) {
    if ( (is.finite(e) && is.finite(f)) || (is.infinite(e) && is.infinite(f))) {
        return(TRUE)
    }
    return(FALSE)
}

#' Find the transitive distance between two graphs as the number of directed paths that differ between graphs x and y
#'
#' @param x An adjacency matrix
#' @param y An adjacency matrix
#' @return Number representing the distance between x and y
#' @export
trans_dist <- function(x,y) {
    # number of edits, considering pathways between genes -- a graph is equivalent to its transitive closure
    stopifnot(all.equal(sort(rownames(x)), sort(colnames(x))))
    stopifnot(all.equal(sort(rownames(y)), sort(colnames(y))))
    stopifnot(all.equal(nrow(x), nrow(y)))

    ix <- igraph::graph_from_adjacency_matrix(x)
    iy <- igraph::graph_from_adjacency_matrix(y[rownames(x), colnames(x)])

    dx <- igraph::distances(ix, mode="out")
    dy <- igraph::distances(iy, mode="out")

    score <- 0
    for (a in rownames(dx)) {
        for (b in rownames(dx)) {
            # is the reachability of a from b the same in both graphs?
            if(! same_reachability(dx[a,b], dy[a,b])) {
                score <- score + 1
            }
        }
    }
    # TODO normalize
    #score <- score / (length(rownames(dx))^2)
    return(score)
}

trans_dist_weighted <- function(x,y) {
    stopifnot(all.equal(sort(rownames(x)), sort(colnames(x))))
    stopifnot(all.equal(sort(rownames(y)), sort(colnames(y))))
    stopifnot(all.equal(nrow(x), nrow(y)))

    ix <- igraph::graph_from_adjacency_matrix(x, weighted=TRUE)
    iy <- igraph::graph_from_adjacency_matrix(y[rownames(x), colnames(x)], weighted=TRUE)

    dx <- igraph::distances(ix, mode="out")
    dy <- igraph::distances(iy, mode="out")

    score <- 0
    for (a in rownames(dx)) {
        for (b in rownames(dx)) {
            # is the reachability of a from b the same in both graphs?
            if(! same_reachability(dx[a,b], dy[a,b])) {
                score <- score + 1
            }
        }
    }
    # TODO normalize
    #score <- score / (length(rownames(dx))^2)
    return(score)

}

#' Generates all distances between graphs in provided list for a list based distance function
#'
#' @param l List of matricies representing graphs of interest
#' @param n Number of graphs to pass to the distance function
#' @param distance_function_list One of the distance metrics defined in this package, must be list based
#' @return Distance matricies
#' @export
generate_distances_list <- function(l, n, distance_function_list) {
    maxlength <- 10000
    l_ig <- lapply(l, igraph::graph_from_adjacency_matrix) # note these are igraphs here before the distance function is called
    comb <- combn(1:length(l), n, simplify=TRUE)
    d <- list()
    if (length(comb) > maxlength) {
        comb <- comb[,sample(1:ncol(comb), maxlength, replace=FALSE)]
    }
    for (i in 1:ncol(comb)) {
        d[[i]] <- distance_function_list(l_ig[comb[,i]])
    }
    return(d)
}

#' Generates all distances between graphs in provided list
#'
#' @param l List of matricies representing graphs of interest
#' @param distance_function One of the distance metrics defined in this package
#' @return Distance matricies
#' @export
generate_distances <- function(l, distance_function) {
    le <- length(l)
    dist <- matrix(0, ncol=le, nrow=le)
    for (i in 1:le) {
        for (j in 1:i) {
            d <- distance_function(l[[j]], l[[i]])
            dist[[i,j]] <- d
            dist[[j,i]] <- d
        }
    }
    return(dist)
}

#' Generates all distances between graphs in provided list, returns in long format
#'
#' @param l List of matricies representing graphs of interest
#' @param distance_function One of the distance metrics defined in this package
#' @return Distance matricies in long format
#' @export
generate_distances_long <- function(l, distance_function) {
    l_dists <- generate_distances(l, distance_function)
    colnames(l_dists) <- paste("graph", 1:nrow(l_dists), sep="")
    l_dists_t <- dplyr::tbl_df(l_dists)
    l_dists_t$graph_name <- paste("graph", 1:nrow(l_dists), sep="")
    l_dists_melt <- reshape2::melt(l_dists_t, id.vars="graph_name") %>% dplyr::rename(graph_name_y = variable, distance=value) %>% dplyr::filter(graph_name != graph_name_y)
    return(l_dists_melt)
}


#' Generates all distances between graphs returned in an mnem object
#'
#' @param res.list As returned from mnem
#' @param distance_function One of the distance metrics defined in this package
#' @return Vector of distance matricies
generate_distances_mnem <- function(res.list, distance_function) {
    dist <- list()
    for (k in 2:length(res.list)) {
        l <- length(res.list[[k]]$comp)
        all.equal(k,l)

        dist[[k]] <- matrix(0, ncol=k, nrow=k)
        for (i in 2:l) {
            for (j in 1:(i-1)) {
                print(paste("i=", i, " j=", j))
                d <- distance_function(res.list[[k]]$comp[[i]]$phi, res.list[[k]]$comp[[j]]$phi)
                dist[[k]][i,j] <- d
                dist[[k]][j,i] <- d
            }
        }
    }
    return(dist)
}

' Plots the networks and a heatmap of distances between them for distances calculated with a list based metric
#'
#' @param l A list of adjacency matricies representing the networks of interest
#' @param distance_function_list One of the distance metrics defined in this package
#' @param filename (optional) File to write plot out to
#' @return Plot object
#' @import magrittr
#' @export
plot_dist_list <- function(l, distance_function_list, filename="") {
    l_dists <- data.frame("idx"=integer(), "n"=integer(), "edges"=integer())
    edges <- unlist(lapply(l, sum))
    n <- rep(1, length(edges))
    idx <- 1:length(edges)
    l_dists <- rbind(l_dists, cbind(idx, cbind(n, edges)))

    r <- 2:length(l)
    for (i in r) {
        print(paste0("i=", i))
        edges <- unlist(generate_distances_list(l, i, distance_function_list))
        n <- rep(i, length(edges))
        idx <- 1:length(edges)
        l_dists <- rbind(l_dists, cbind(idx, cbind(n, edges)))
        if (sum(edges) == 0) {
            break # if all n-1 intersections are 0 then n and n+k will all be 0 also
        }
    }
    #ggplot2::ggplot(l_dists, aes(x=edges)) + geom_histogram() + facet_wrap(n ~ .)
    #ggplot2::ggplot(l_dists, aes(x=idx, y=edges)) + geom_point() + facet_wrap(n ~ .)
    ggplot2::ggplot(l_dists, ggplot2::aes(x=n, y=edges, group=n)) + ggplot2::geom_boxplot()
    if (filename != "") {
        ggplot2::ggsave(filename)
    }
}



#' Plots the networks and a heatmap of distances between them
#'
#' @param l A list of adjacency matricies representing the networks of interest
#' @param distance_function One of the distance metrics defined in this package
#' @param filename (optional) File to write plot out to
#' @return Plot object
#' @import magrittr
#' @export
plot_dist <- function(l, distance_function, filename="", draw_networks=NULL) {
    l_dists <- generate_distances(l, distance_function)
    generate_networks <- FALSE
    if (nrow(l_dists) <= 5) { # TODO: set default value for package
        generate_networks <- TRUE
    }
    if (! is.null(draw_networks)) {
        generate_networks <- draw_networks
    }
    if (generate_networks) {
        # generate network images
        images <- l
        for (i in 1:nrow(l_dists)) {
            net <- network::as.network(l[[i]], matrix.type="adjacency", directed=TRUE)
            palette <- RColorBrewer::brewer.pal(nrow(l[[1]]), "Paired")
            images[[i]] <- GGally::ggnet2(net, mode="circle", label=TRUE, alpha=0.7, size=3, color=palette, arrow.size=4, arrow.gap=.06, edge.color = "grey30", label.size=1.5, layout.exp = 0.3) + ggplot2::theme(aspect.ratio=1)
        }
        plot_networks <- do.call(gridExtra::grid.arrange, c(images, ncol = length(images)))
    }

    # turn matrix into long data format
    colnames(l_dists) <- paste("graph", 1:nrow(l_dists), sep="")
    l_dists_t <- dplyr::tbl_df(l_dists)
    l_dists_t$graph_name <- paste("graph", 1:nrow(l_dists), sep="")
    l_dists_melt <- reshape2::melt(l_dists_t, id.vars="graph_name") %>% dplyr::rename(graph_name_y = variable, distance=value) %>% dplyr::filter(graph_name != graph_name_y)

    # generate the heatmap
    if (generate_networks) {
        plot_heatmap <- ggplot2::ggplot(l_dists_melt, ggplot2::aes(x=graph_name, y=graph_name_y, label=distance, fill=distance)) + 
            ggplot2::geom_tile() + 
            ggplot2::geom_text() + 
            ggplot2::theme_minimal() + 
            ggplot2::scale_y_discrete(labels=rev(names(l)), limits=rev(levels(l_dists_melt$graph_name_y))) + 
            ggplot2::scale_x_discrete(labels=c()) + 
            ggplot2::scale_fill_distiller(guide=FALSE, palette="Blues", direction=1) +
            ggplot2::coord_fixed() +
            ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), axis.ticks=ggplot2::element_blank()) # axis.ticks isn't working??
    } else {
        plot_histogram <- ggplot2::ggplot(l_dists_melt, ggplot2::aes(x=distance)) + ggplot2::geom_bar() + ggplot2::theme_bw()
    }

    # is there a better way to save the output if a filename is provided?
    if (filename != "") {
        grDevices::pdf(filename)
    }
    if (generate_networks) {
        gridExtra::grid.arrange(plot_networks, plot_heatmap, nrow=2)
    } else {
        gridExtra::grid.arrange(plot_histogram, nrow=1)
    }
    if (filename != "") {
        grDevices::dev.off()
    }
    return(l_dists_melt)
}

#' Calculates the distances between k random networks of n nodes each
#'
#' @param distance_function One of the distance metrics defined in this package
#' @return A network as an adjacency matrix
#' @export
calculate_random_dist_list <- function(distance_function_list) {
    for (p in seq(0.2, 0.5, 0.1)) { #random probability
        n <- 10
        k <- 25
        filename <- paste0("~/projects/simulate/random_graph_distance/intersectdist_randomnet_n", n, "_k", k, "_p", p, ".pdf")
        l <- generate_random_networks(n, k, p, distance_function, filename=filename)
        plot_dist_list(l, distance_function_list, filename=filename)
    }
}

#' Calculates the distances between k scale free networks of n nodes each
#'
#' @param distance_function One of the distance metrics defined in this package
#' @return A network as an adjacency matrix
#' @export
calculate_barabasi_dist_list <- function(distance_function_list) {
    for (p in seq(1.2, 2.0, 0.1)) { #barabasi power
        n <- 10
        k <- 25
        filename <- paste0("~/projects/simulate/random_graph_distance/barabasinet_n", n, "_k", k, "_p", p, ".pdf")
        l <- generate_barabasi_networks(n, k, p, distance_function, filename=filename)
        plot_dist_list(l, distance_function_list, filename=filename)
    }
}


#' Calculates the distances between k random networks of n nodes each
#'
#' @param n The number of nodes desired
#' @param k The number of random networks
#' @param p The probability of generating an edge between two nodes
#' @param distance_function One of the distance metrics defined in this package
#' @return A network as an adjacency matrix
#' @export
generate_random_networks <- function(n, k, p, distance_function, filename="") {
    l <- list()
    for (i in 1:k) {
        ig <- igraph::erdos.renyi.game(n, p, type="gnp", directed=TRUE, loops=FALSE)
        #ig <- igraph::barabasi.game(n, power=p, directed=TRUE)
        m <- as.matrix(igraph::get.adjacency(ig))
        rownames(m) <- letters[1:n]
        colnames(m) <- rownames(m)
        l[[i]] <- m
    }
    return(l)
}

#' Calculates the distances between k random networks of n nodes each
#'
#' @param n The number of nodes desired
#' @param k The number of random networks
#' @param p The probability of generating an edge between two nodes
#' @param distance_function One of the distance metrics defined in this package
#' @return A network as an adjacency matrix
#' @export
generate_barabasi_networks <- function(n, k, p, distance_function, filename="") {
    l <- list()
    for (i in 1:k) {
        ig <- igraph::barabasi.game(n, power=p, directed=TRUE)
        m <- as.matrix(igraph::get.adjacency(ig))
        rownames(m) <- letters[1:n]
        colnames(m) <- rownames(m)
        l[[i]] <- m
    }
    return(l)
}

#' Find the largest subgraph shared by the arguments
#'
#' @param x 
#' @param y 
#' @param n
#' @return number of shared edges
#' @export
intersect_dist <- function(x, y) {
    ig <- igraph::intersection(igraph::graph_from_adjacency_matrix(x), igraph::graph_from_adjacency_matrix(y), byname=TRUE)
    return(igraph::gsize(ig))
}

#' Find the largest subgraph shared by the networks in the list
#'
#' @param l_ig A list of igraph objects
#' @return number of edges in the largest subgraph of all graphs
#' @export
intersect_dist_list <- function(l_ig) {
    ig <- l_ig[[1]]
    for (i in 2:length(l_ig)) {
        ig <- igraph::intersection(ig, l_ig[[i]], byname=TRUE)
    }
    return(igraph::gsize(ig))
}


