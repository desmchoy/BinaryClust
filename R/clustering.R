#' Performing binary classification
#'
#' This function performs binary classification on loaded data using a user-specified list of criteria. Returns a column in data frame class.
#' @param data (Required) Data in data frame format. Expects results from load_data
#' @param class.file (Required) Path to classification file in comma-separated format. Please see tutorial for details.
#' @param scale Standardisation of data (Default: FALSE)
#' @export
#' @examples
#' binary.results <- binary_class(data, class.file = 'cell_types.csv')
binary_class <- function ( data = NULL,
	class.file = NULL,
	scale = FALSE ) {

	if ( is.null(data) ) {
		return ('Please provide data. Use load_data')
	}

	if ( is.null(class.file) ) {
		return ('Please provide a classification file.')
	} else {

		data <- data[, colSums(data != 0, na.rm = TRUE) > 0]
		
		if ( scale == TRUE ) {
			data <- data.frame(scale(data))
		} else if ( scale == FALSE ) {
			data <- data
		} else {
			return ('Please choose either TRUE or FALSE for scale.')
		}

		kmeans.results <- data.frame(apply(data, 2, run_kmeans))
		colnames(kmeans.results) <- colnames(data)
		kmeans.limits <- data.frame(matrix(nrow = 0, ncol = 3))
		colnames(kmeans.limits) <- c('pos_min', 'Q1', 'Q3')
	
		for ( marker in colnames(kmeans.results) ) {
	
			if ( length(unique(kmeans.results[,marker])) == 1 & unique(kmeans.results[,marker])[1] == '-') {
				pos_min <- 999
				Q1 <- 0
				Q3 <- 0
			} else {
				pos_min <- as.numeric(min(data[kmeans.results[,marker] == '+',marker]))
				Q1 <- as.numeric(quantile(data[kmeans.results[,marker] == '+',marker], 0.25))
				Q3 <- as.numeric(quantile(data[kmeans.results[,marker] == '+',marker], 0.75))
			}
			kmeans.results[,marker] <- as.character(kmeans.results[,marker])
			kmeans.results[data[,marker] <= Q1 & kmeans.results[,marker] == '+',marker] <- 'L'
			kmeans.results[data[,marker] > Q1 & data[,marker] <= Q3,marker] <- 'I'
			kmeans.results[data[,marker] > Q3,marker] <- 'H'
			kmeans.limits[marker,] <- c(pos_min, Q1, Q3)
	
		}


			kmeans.results[,'Number'] <- rownames(kmeans.results)
			class.results <- data.frame(matrix(nrow = dim(data)[1], ncol = 1))
			colnames(class.results) <- c('Cell.Type')
			class.cat <- read.table(class.file, header = TRUE, sep = ',', row.names = 1, check.names = FALSE)
			class.cat <- apply(class.cat, c(1,2), function (x) { gsub(' ', '', x) })
			class.cat[class.cat == ''] <- 'A'
	
			for ( cell.type in rownames(class.cat) ) {
	
				temp.row <- class.cat[cell.type,, drop = FALSE]
				temp.row <- temp.row[,temp.row != 'A', drop = FALSE]
	
				temp.row.fine <- temp.row[,temp.row != '+', drop = FALSE]
				temp.row.all <- temp.row[,temp.row == '+', drop = FALSE]
	
				rows.required <- merge(kmeans.results, temp.row.fine)
				if ( length(colnames(temp.row.all)) == 1 ) {
					rows.required <- rows.required[rows.required[,colnames(temp.row.all)] %in% c('L', 'I', 'H'),][,'Number']
				} else {
					rows.required <- rows.required[apply(rows.required[,colnames(temp.row.all)], 1, function (x) { all(x %in% c('L', 'I', 'H')) }),][,'Number']
				}
				class.results[rows.required, 'Cell.Type'] <- cell.type
	
			}
	
			class.results[is.na(class.results[,'Cell.Type']),'Cell.Type'] <- 'Unclassified'
	
		return ( class.results )

	}

}




run_kmeans <- function ( data,
	cutoff.limit = 2 ) {

	k.means <- kmeans(data, 2, iter.max = 500)
	pos <- which(k.means$centers == max(k.means$centers))
	neg <- which(k.means$centers == min(k.means$centers))
	results <- k.means$cluster

	if ( min(k.means$centers) > cutoff.limit | max(k.means$centers) < cutoff.limit | abs(max(k.means$centers) - min(k.means$centers)) < 1 | abs(max(data) - max(k.means$centers)) > cutoff.limit ) {

		results[data > cutoff.limit] <- '+'
		results[data <= cutoff.limit] <- '-'

	} else {

		results[results == pos] <- '+'
		results[results == neg] <- '-'

	}

	return ( results )

}



#' Summarising binary classification results
#'
#' This function summarises binary classification results from binary_class. Returns a data frame of median expression values of all markers, abundance and percentage abundance for all user-defined cell types
#' @param data (Required) Data in data frame format. Expects results from load_data
#' @param binary.results (Required) Cell type results in data frame format. Expects results from binary_class
#' @export
#' @examples
#' binary.summary <- binary_summary(data, binary.results)
binary_summary <- function ( data = NULL,
	binary.results = NULL ) {

	if ( is.null(data) ) {
		return ('Please provide data. Use load_data')
	}

	if ( is.null(binary.results) ) {
		return ('Please provide binary classification results. Use binary_class')
	}

	frequencies <- data.frame(table(binary.results[,'Cell.Type']))
	colnames(frequencies) <- c('Cell.Type', 'Frequency')
	frequencies[,'Percentage'] <- frequencies[,'Frequency'] * 100 / sum(frequencies[,'Frequency'])

	medians <- data.frame(matrix(nrow = dim(frequencies)[1], ncol = dim(data)[2]))
	rownames(medians) <- frequencies[,1]
	colnames(medians) <- colnames(data)

	for ( i in frequencies[,1] ) {
		medians[i,] <- as.vector(apply(data[binary.results[,'Cell.Type'] == i,], 2, median))
	}

	medians[,'Cell.Type'] <- rownames(medians)

	to.return <- merge(medians, frequencies, by = 'Cell.Type')

	return ( to.return )

}



#' Clustering binary classified results
#'
#' This function clusters cells individually for each user-definied cell type. Returns a data frame.
#' @param data (Required) Data in data frame format. Expects results from load_data
#' @param binary.results (Required) Cell type results in data frame format. Expects results from binary_class
#' @param method Clustering method. k-means (Default) or FlowSOM
#' @param no_of_clusters Default: 40
#' @export
#' @examples
#' cluster.results <- cluster_subsets(data, binary.results)
#' cluster.results <- cluster_subsets(data, binary.results, method = 'FlowSOM')
cluster_subsets <- function ( data = NULL,
	binary.results = NULL,
	method = 'kmeans',
	no_of_clusters = 40 ) {
	
	if ( is.null(data) ) {
		return ('Please provide data. Use load_data')
	}

	if ( is.null(binary.results) ) {
		return ('Please provide binary classification results. Use binary_class')
	}


	binary.results[,'Cluster'] <- NA
	scaled.data <- scale(data)
	scaled.data <- data.frame(scaled.data)

	for ( cell.type in sort(unique(binary.results[,'Cell.Type'])) ) {

		if ( method == 'kmeans' ) {
			binary.results <- run_kmeans_cluster(scaled.data, binary.results, cell.type, no_of_clusters)
		} else if ( method == 'FlowSOM' ) {
			binary.results <- run_FlowSOM_cluster(scaled.data, binary.results, cell.type, no_of_clusters)
		} else {
			return ('Please choose a valid clustering method - kmeans or FlowSOM')
		}
			
	}

	cluster.results <- binary.results
#	cluster.results[,'Cell.Subtype'] <- ''

	return ( cluster.results )

}



run_kmeans_cluster <- function ( data,
	binary.results,
	cell.type,
	no_of_clusters ) {

	data.subset <- data[which(binary.results[,'Cell.Type'] == cell.type),]
	kmeans.nClus <- no_of_clusters

	if ( dim(data.subset)[1] == 0 ) {
		print(paste0('There are no cells in the subset ', cell.type))
	} else if ( dim(data.subset)[1] < kmeans.nClus ) {
		dummy.number <- 1
		for ( i in rownames(data.subset) ) {
			binary.results[i,'Cluster'] <- dummy.number
			dummy.number <- dummy.number + 1
		}
		print(paste0('There are less than ', kmeans.nClus, ' cells in the subset ', cell.type, '. No clustering is performed.'))
	} else {
		temp.kmeans <- kmeans(data.subset, kmeans.nClus, iter.max = 50)
		for ( i in names(temp.kmeans$cluster) ) {
			binary.results[as.numeric(i),'Cluster'] <- temp.kmeans$cluster[i]
		}
	}
	return ( binary.results )
}



run_FlowSOM_cluster <- function ( data,
	binary.results,
	cell.type,
	no_of_clusters ) {

	data.subset <- data[which(binary.results[,'Cell.Type'] == cell.type),]
	fSOM.xdim <- 10
	fSOM.ydim <- 10
	fSOM.nClus <- no_of_clusters

	if ( dim(data.subset)[1] == 0 ) {

		print(paste0('There are no cells in the subset ', cell.type))

	} else if ( dim(data.subset)[1] < 100 ) {

		print(paste0('Population is too small: ', dim(data.subset)[1], ', FlowSOM cannot be run for ', cell.type, '. K-means is run instead.'))
		binary.results <- run_kmeans_cluster(scaled.data, binary.results, cell.type, no_of_clusters)

	} else {

		flowframe <- new('flowFrame', exprs = as.matrix(data.subset))
		fSOM <- FlowSOM(flowframe,
			transform = FALSE,
			scale = FALSE,
			colsToUse = colnames(data.subset),
			xdim = fSOM.xdim,
			ydim = fSOM.ydim,
			nClus = fSOM.nClus)
		clusters <- GetClusters(fSOM)
		metaclusters <- GetMetaclusters(fSOM)
		binary.results[rownames(data.subset),'Cluster'] <- metaclusters
	}

	return ( binary.results )

}


#' Summarising clustering results
#'
#' This function summarises post-binary classification clustering results. Returns a data frame of median expressions and abundance for each cluster of each user-defined cell type.
#' @param data (Required) Data in data frame format. Expects results from load_data
#' @param cluster.results (Required) Cell type results in data frame format. Expects results from cluster_subsets
#' @export
#' @examples
#' cluster.summary <- cluster_summary(data, cluster.results)
cluster_summary <- function ( data = NULL,
	cluster.results = NULL	) {

	if ( is.null(data) ) {
		return ('Please provide data. Use load_data')
	}

	if ( is.null(cluster.results) ) {
		return ('Please provide cluster results. Use cluster_subset')
	}

	colnames(cluster.results) <- c('Cell.Type', 'Cluster')
	combined <- cbind(data, cluster.results[,c('Cell.Type', 'Cluster')])
	medians <- combined %>%	group_by(Cell.Type, Cluster) %>%
		summarise_all(median)
	medians <- as.data.frame(medians)

	percentage <- combined %>% group_by(Cell.Type, Cluster) %>%
		summarise(Frequency = n(), Percentage = (100 * n() / dim(combined)[1]))
	percentage <- as.data.frame(percentage)

	to.return <- merge(medians, percentage, by = c('Cell.Type', 'Cluster'))
	to.return['Cell.Subtype'] <- ''
	
	return ( to.return )
	
}



#' Updating clustering results with filled subtype column
#'
#' This function takes in the filled cell subtype column of a clustering summary and updates the associated all-cell clustering results. This function is intended for the UMAP plotting function. Returns a data frame.
#' @param cluster.results (Required) Cell type results in data frame format. Expects results from cluster_subsets
#' @param completed.summary (Required) Cluster summary in data frame format. Expects results from cluster_summary with a completed Cell.Subtype column
#' @export
#' @examples
#' updated.cluster.results <- update_cluster_results(cluster.results, completed.summary)
update_cluster_results <- function ( cluster.results = NULL,
	completed.summary = NULL) {

	if ( is.null(cluster.results) ) {
		return ('Please provide cluster results. Use cluster_subsets')
	}

	if ( is.null(completed.summary) ) {
		return ('Please provide cluster summary results. Use cluster_summary')
	}


	cluster.results[,'Number'] <- as.numeric(rownames(cluster.results))
	updated.results <- merge(cluster.results, completed.summary[,c('Cell.Type', 'Cluster', 'Cell.Subtype')], by = c('Cell.Type', 'Cluster'))
	updated.results <- updated.results[order(updated.results[,'Number']),]
	rownames(updated.results) <- NULL
	updated.results <- updated.results[,!(colnames(updated.results) %in% 'Number')]

	return ( updated.results )

}



