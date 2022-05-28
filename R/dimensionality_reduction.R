#' Calculating t-SNE coordinates for data
#'
#' This function calculates t-SNE coordinates on data. Random seed is set to 1. Returns full t-SNE object.
#' @param data (Required) Data in data frame format. Expects results from load_data
#' @param scale Default: TRUE. Option to standardise data prior to t-SNE calculations
#' @param center Default: FALSE. Option to center data prior to PCA within the Rtsne algorithm
#' @param perplexity Value for perplexity (Default: 30)
#' @param n_iter Number of iterations (Default: 1000)
#' @export
#' @examples
#' tsne.results <- run_TSNE(data)
#'
#' t-SNE coordinates can be extracted under Y:
#' tsne.coordinates <- tsne.results$Y
run_TSNE <- function ( data = NULL,
	scale = TRUE,
	center = FALSE,
	perplexity = 30,
	n_iter = 1000 ) {

	if ( is.null(data) ) {
		return ('Please provide data. Use load_data')
	}

	set.seed(1)

	if ( scale == TRUE ) {
		to.run <- scale(data)
	} else if ( scale == FALSE ) {
		to.run <- data
	} else {
		return('Please enter TRUE or FALSE for scale')
	}

	tsne.results <- Rtsne(to.run, perplexity = perplexity, max_iter = n_iter, pca_center = center, check_duplicates = FALSE)

	return ( tsne.results )

}


#' Calculating UMAP coordinates for data
#'
#' This function calculates UMAP coordinates on data. Random seed is set to 1. Returns full UMAP object.
#' @param data (Required) Data in data frame format. Expects results from load_data
#' @param scale Default: TRUE. Option to standardise data prior to UMAP calculations
#' @export
#' @examples
#' umap.results <- run_UMAP(data)
#'
#' UMAP coordinates can be extracted under layout:
#' umap.coordinates <- umap.results$layout
run_UMAP <- function ( data = NULL,
	scale = TRUE ) {

	if ( is.null(data) ) {
		return ('Please provide data. Use load_data')
	}

	set.seed(1)
	
	if ( scale == TRUE ) {
		to.run <- scale(data)
	} else if ( scale == FALSE ) {
		to.run <- data
	} else {
		return('Please enter TRUE or FALSE for scale')
	}

	umap.results <- umap(to.run)

	return ( umap.results )

}
