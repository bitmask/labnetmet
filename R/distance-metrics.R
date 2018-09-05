#' Find the Hamming distance between two graphs as sum of the absolute value differences between adjacency matricies.
#'
#' @param x An adjacency matrix
#' @param y An adjacency matrix
#' @return Number representing the distance between x and y
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
subgraph_dist <- function(x,y) {
    # x and y are full matricies

    # find all subgraphs of x and y, then compare them

    stopifnot(all.equal(sort(rownames(x)), sort(colnames(x))))
    stopifnot(all.equal(sort(rownames(y)), sort(colnames(y))))
    stopifnot(all.equal(nrow(x), nrow(y)))

    min_subgraph_size <- 2
    max_subgraph_size <- max(nrow(y) - 1, min_subgraph_size)
    subgraph_sizes <- min_subgraph_size:max_subgraph_size

    subgraphs <- unlist(lapply(subgraph_sizes, function(size) {combn(rownames(y), size, simplify=F)}), recursive=FALSE)

    ix <- igraph::graph_from_adjacency_matrix(x)
    iy <- igraph::graph_from_adjacency_matrix(y[rownames(x), colnames(x)])

    isomorphic <- unlist(lapply(subgraphs, function(subgraph) {isomorphic(induced_subgraph(ix, subgraph), induced_subgraph(iy, subgraph))}))

    # normalize result to length of number of subgraphs
    return ((length(isomorphic) - sum(isomorphic)) / length(isomorphic))

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

#' Generates all distances between graphs in provided list
#'
#' @param l List of matricies representing graphs of interest
#' @param distance_function One of the distance metrics defined in this package
#' @return Distance matricies
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

#' Plots the networks and a heatmap of distances between them
#'
#' @param l A list of adjacency matricies representing the networks of interest
#' @param distance_function One of the distance metrics defined in this package
#' @return Plot object
plot_dist <- function(l, distance_function, filename="") {
    l_dists <- generate_distances(l, distance_function)

    # generate network images
    images <- l
    for (i in 1:nrow(l_dists)) {
        net <- network::as.network(l[[i]], matrix.type="adjacency", directed=TRUE)
        # TODO: scale the sizes as a function of the number of nodes and networks
        images[[i]] <- GGally::ggnet2(net, mode="circle", label=TRUE, alpha=0, arrow.size=3, arrow.gap=.02, label.size=2)
    }
    plot_networks <- do.call(gridExtra::grid.arrange, c(images, ncol = length(images)))

    # turn matrix into long data format
    colnames(l_dists) <- letters[1:nrow(l_dists)]
    l_dists_t <- dplyr::tbl_df(l_dists)
    l_dists_t$graph_name <- letters[1:nrow(l_dists)]
    l_dists_melt <- reshape2::melt(l_dists_t, id.vars="graph_name") %>% dplyr::rename(graph_name_y = variable)

    # generate the heatmap
    plot_heatmap <- ggplot2::ggplot(l_dists_melt, ggplot2::aes(x=graph_name, y=graph_name_y, label=value, fill=value)) + ggplot2::geom_tile() + ggplot2::geom_text() + ggplot2::scale_y_discrete(limits=rev(levels(l_dists_melt$graph_name_y))) + ggplot2::theme_void() + ggplot2::coord_fixed() + ggplot2::scale_fill_distiller(guide=FALSE, palette="Blues", direction=1) 
    
    # there must be a better way to do this, but everything else seems broken
    if (filename != "") {
        pdf(filename)
        gridExtra::grid.arrange(plot_networks, plot_heatmap, nrow=2)
        dev.off()
    }
    gridExtra::grid.arrange(plot_networks, plot_heatmap, nrow=2)
}
