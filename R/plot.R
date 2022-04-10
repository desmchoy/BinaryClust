#' Calculating and plotting non-redundancy scores (NRS)
#'
#' This function calculates and plots non-redundancy scores of all the markers in data. Returns a graph.
#' @param data (Required) Data in data frame format. Expects results from load_data
#' @export
#' @examples
#' plot_NRS(data)
plot_NRS <- function ( data = NULL ) {

	if ( is.null(data) ) {
		return ('Please provide data (use load_data).')
	}


	calc_NRS <- function (x) {
		if ( dim(x)[1] < 3 ) {
			NRS <- data.frame(matrix(nrow=ncol(x), ncol=1))
			NRS[,1] <- NA
			return (NRS)
		} else {
			pc <- prcomp(x, center = TRUE, scale. = FALSE)
			NRS <- rowSums(abs(pc$rotation[, seq_len(3)]) * outer(rep(1, ncol(x)), pc$sdev[seq_len(3)]^2))
			NRS <- data.frame(NRS)
			return (NRS)
		}
	}

	NRS <- calc_NRS(data)
	NRS[,'Marker'] <- rownames(NRS)

	plot <- NRS %>% ggplot(aes(x=reorder(Marker, desc(NRS)), y=NRS)) +
		geom_col() +
		xlab('') +
		ylab('NRS') +
		theme_minimal() +
		theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			plot.title = element_text(face = 'bold', hjust = 0.5),
			axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5))

		
	return ( plot )

}



#' Plotting abundances based on binary classification results
#'
#' This function plots the abundanceresults of binary classification. Returns a graph.
#' @param binary.summary (Required) A data frame containing abundance results. Expects results from binary_summary
#' @param remove.unclass Option to remove unclassified cells from plot. (Default: FALSE)
#' @param xlabel Label for x-axis (Default: NULL)
#' @param ylabel Label for y-axis (Default: 'Percentage Abundance')
#' @param title Plot title (Default: 'Cell Abundances')
#' @export
#' @examples
#' plot_binary_abundances(binary.summary)
plot_binary_abundances <- function ( binary.summary = NULL,
	remove.unclass = FALSE,
	xlabel = '',
	ylabel = 'Percentage Abundance',
	title = 'Cell Abundances' ) {

	if ( is.null(binary.summary) ) {
		return ('Please provide binary classification summary (use binary_summary).')
	}


	if ( remove.unclass == FALSE ) {
		to.plot <- binary.summary
	} else if ( remove.unclass == TRUE ) {
		to.plot <- binary.summary[binary.summary[,'Cell.Type'] != 'Unclassified',]
	} else {
		return('Please choose TRUE or FALSE for remove.unclass')
	}

	plot <- to.plot %>% ggplot(aes(x=reorder(`Cell.Type`, desc(`Cell.Type`)), y=`Percentage`)) +
		xlab(xlabel) +
		ylab(ylabel) +
		ggtitle(title) +
		geom_col() +
		coord_flip() +
		theme_minimal() +
		theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			plot.title = element_text(face = 'bold', hjust = 0.5))

	return ( plot )

}



#' Plotting median expressions based on binary classification results
#'
#' This function plots the median expressions of all markers after binary classification. Returns a graph.
#' @param binary.summary (Required) A data frame containing abundance results. Expects results from binary_summary
#' @param remove.unclass Option to remove unclassified cells from plot. (Default: FALSE)
#' @param title Plot title (Default: 'Median Expressions')
#' @export
#' @examples
#' plot_binary_abundances(binary.summary)
plot_binary_median <- function ( binary.summary = NULL,
	remove.unclass = FALSE,
	title = 'Median Expressions' ) {

	if ( is.null(binary.summary) ) {
		return ('Please provide binary classification summary (use binary_summary).')
	}

	if ( remove.unclass == FALSE ) {
		to.plot <- binary.summary
	} else if ( remove.unclass == TRUE ) {
		to.plot <- binary.summary[binary.summary[,'Cell.Type'] != 'Unclassified',]
	}

	rownames(to.plot) <- to.plot[,'Cell.Type']
	to.plot <- subset(to.plot, select = -c(`Cell.Type`, Frequency, `Percentage`))

	medians.heatmap <- pheatmap(to.plot,
		cluster_rows = FALSE,
		cluster_cols = FALSE,
		cellwidth = 10,
		cellheight = 10,
		main = title)
	
	grid.newpage()
	plot <- grid.draw(medians.heatmap$gtable)

	return ( plot )

}



#' Plotting binary classification results on t-SNE coordinates
#'
#' This function plots t-SNE data coloured by binary classification results. Returns a graph.
#' @param tsne.results (Required) t-SNE object. Expects results from run_TSNE
#' @param binary.results (Required) Results from binary classification in data frame format. Expects results from binary_class
#' @param column Choose between plotting Cell.Type or Cell.Subtype (Default: Cell.Type) Cell.Subtype only availble for updated clustering results
#' @param title Title of the plot (Default: NULL)
#' @param facet Breaks plot into facets (Default: NULL)
#' @export
#' @examples
#' plot_TSNE(tsne.results, binary.results)
#' plot_TSNE(tsne.results, updated.cluster.results, column = 'Cell.Subtype', facet = TRUE)
plot_TSNE <- function ( tsne.results = NULL,
	binary.results = NULL,
	column = 'Cell.Type',
	title = '',
	facet = FALSE ) {

	if ( is.null(tsne.results) ) {
		return ('Please provide t-SNE object (use run_TSNE).')
	}

	if ( is.null(binary.results) ) {
		return ('Please provide binary classification results (use binary_class).')
	}

	if ( column == 'Cell.Type' ) {
		legend.label <- 'Cell Type'
		if ( column %in% colnames(binary.results) == FALSE ) {
			return ('Input file does not contain Cell.Type')
		}
	} else if ( column == 'Cell.Subtype' ) {
		legend.label <- 'Cell Subtype'
		if ( column %in% colnames(binary.results) == FALSE ) {
			return ('Input file does not contain Cell.Subtype')
		}
	} else {
		return ('Please choose a valid option for column - Cell.Type or Cell.Subtype')
	}


	to.plot <- cbind(tsne.results$Y, binary.results)
	colnames(to.plot)[c(1,2)] <- c('TSNE1', 'TSNE2')

	plot <- to.plot %>% ggplot(aes(x=TSNE1, y=TSNE2, color=.data[[column]])) +
		geom_point(size=0.05) +
		theme_minimal() +
		xlab('t-SNE 1') +
		ylab('t-SNE 2') +
		ggtitle(title) +
		scale_color_discrete(name = legend.label) +
		theme(legend.position = 'right',
			plot.title = element_text(hjust=0.5),
			panel.border = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.line = element_line(colour = "black")) +
		guides(colour = guide_legend(override.aes = list(size=2.5)))


	if ( facet == TRUE ) {
		plot <- plot + facet_wrap(vars(.data[[column]]))
	}

	return ( plot )

}



#' Plotting binary classification results on UMAP coordinates
#'
#' This function plots UMAP data coloured by binary classification results. Returns a graph.
#' @param umap.results (Required) UMAP object. Expects results from run_UMAP
#' @param binary.results (Required) Results from binary classification in data frame format. Expects results from binary_class
#' @param column Choose between plotting Cell.Type or Cell.Subtype (Default: Cell.Type) Cell.Subtype only availble for updated clustering results
#' @param title Title of the plot (Default: NULL)
#' @param facet Breaks plot into facets (Default: NULL)
#' @export
#' @examples
#' plot_UMAP(umap.results, binary.results)
#' plot_UMAP(umap.results, updated.cluster.results, column = 'Cell.Subtype', facet = TRUE)
plot_UMAP <- function ( umap.results = NULL,
	binary.results = NULL,
	column = 'Cell.Type',
	title = '',
	facet = FALSE ) {

	if ( is.null(umap.results) ) {
		return ('Please provide UMAP object (use run_UMAP).')
	}

	if ( is.null(binary.results) ) {
		return ('Please provide binary classification results (use binary_class).')
	}

	if ( column == 'Cell.Type' ) {
		legend.label <- 'Cell Type'
		if ( column %in% colnames(binary.results) == FALSE ) {
			return ('Input file does not contain Cell.Type')
		}
	} else if ( column == 'Cell.Subtype' ) {
		legend.label <- 'Cell Subtype'
		if ( column %in% colnames(binary.results) == FALSE ) {
			return ('Input file does not contain Cell.Subtype')
		}
	} else {
		return ('Please choose a valid option for column - Cell.Type or Cell.Subtype')
	}


	to.plot <- cbind(umap.results$layout, binary.results)
	colnames(to.plot)[c(1,2)] <- c('UMAP1', 'UMAP2')
	
	plot <- to.plot %>% ggplot(aes(x=UMAP1, y=UMAP2, color=.data[[column]])) +
		geom_point(size=0.05) +
		theme_minimal() +
		xlab('UMAP 1') +
		ylab('UMAP 2') +
		ggtitle(title) +
		scale_color_discrete(name = legend.label) +
		theme(legend.position = 'right',
			plot.title = element_text(hjust=0.5),
			panel.border = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.line = element_line(colour = "black")) +
		guides(colour = guide_legend(override.aes = list(size=2.5)))


	if ( facet == TRUE ) {
		plot <- plot + facet_wrap(vars(.data[[column]]))
	}

	return ( plot )

}



#' Colouring t-SNE results with marker expressions
#'
#' This function colours t-SNE results by marker expressions. Returns a graph.
#' @param tsne.results (Required) t-SNE object. Expects results from run_TSNE
#' @param data (Required) Data in data frame format. Expects results from load_data
#' @param marker (Required) Marker to be coloured
#' @param lower Lower bound of colour spectrum (Default: Lowest value of marker expression)
#' @param upper Upper bound of colour spectrum (Default: Highest value of marker expression)
#' @param palette ggplot2 colour spectrum (Default: 'Purples')
#' @param title Title of plot (Default: NULL)
#' @export
#' @examples
#' plot_TSNE_marker(tsne.results, data, marker = 'CD45')
plot_TSNE_marker <- function ( tsne.results = NULL,
	data = NULL,
	marker = NULL,
	lower = NULL,
	upper = NULL,
	palette = 'Purples',
	title = '' ) {

	if ( is.null(tsne.results) ) {
		return ('Please provide t-SNE object (use run_TSNE).')
	}

	if ( is.null(data) ) {
		return ('Please provide data (use load_data).')
	}

	if ( is.null(marker) ) {
		return ('Please specify a marker to plot.')
	}

	if ( is.null(lower) ) {
		lower.limit <- floor(min(data[,marker]))
	} else {
		lower.limit <- lower
	}

	if ( is.null(upper) ) {
		upper.limit <- ceiling(max(data[,marker]))
	} else {
		upper.limit <- upper
	}

	marker.name <- marker
#	marker <- gsub('^[^_]*_(*)', '\\1', marker)

	to.plot <- cbind(tsne.results$Y, data)
	colnames(to.plot)[c(1,2)] <- c('TSNE1', 'TSNE2')

	plot <- to.plot %>% ggplot(aes(x=TSNE1, y=TSNE2, color=.data[[marker]])) +
		geom_point(size=0.05) +
		theme_minimal() +
		xlab('t-SNE 1') +
		ylab('t-SNE 2') +
		ggtitle(title) +
		scale_colour_distiller(palette = palette, direction = 1, limits = c(lower.limit,upper.limit)) +
		theme(legend.position = 'right',
			plot.title = element_text(hjust=0.5),
			panel.border = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.line = element_line(colour = "black"))

	return ( plot )

}



#' Colouring UMAP results with marker expressions
#'
#' This function colours UMAP results by marker expressions. Returns a graph.
#' @param umap.results (Required) UMAP object. Expects results from run_UMAP
#' @param data (Required) Data in data frame format. Expects results from load_data
#' @param marker (Required) Marker to be coloured
#' @param lower Lower bound of colour spectrum (Default: Lowest value of marker expression)
#' @param upper Upper bound of colour spectrum (Default: Highest value of marker expression)
#' @param palette ggplot2 colour spectrum (Default: 'Purples')
#' @param title Title of plot (Default: NULL)
#' @export
#' @examples
#' plot_UMAP_marker(umap.results, data, marker = 'CD45')
plot_UMAP_marker <- function ( umap.results = NULL,
	data = NULL,
	marker = NULL,
	lower = NULL,
	upper = NULL,
	palette = 'Purples',
	title = '' ) {

	if ( is.null(umap.results) ) {
		return ('Please provide UMAP object (use run_UMAP).')
	}

	if ( is.null(data) ) {
		return ('Please provide data (use load_data).')
	}

	if ( is.null(marker) ) {
		return ('Please specify a marker to plot.')
	}

	if ( is.null(lower) ) {
		lower.limit <- floor(min(data[,marker]))
	} else {
		lower.limit <- lower
	}

	if ( is.null(upper) ) {
		upper.limit <- ceiling(max(data[,marker]))
	} else {
		upper.limit <- upper
	}

	marker.name <- marker
#	marker <- gsub('^[^_]*_(*)', '\\1', marker)

	to.plot <- cbind(umap.results$layout, data)
	colnames(to.plot)[c(1,2)] <- c('UMAP1', 'UMAP2')

	plot <- to.plot %>% ggplot(aes(x=UMAP1, y=UMAP2, color=.data[[marker]])) +
		geom_point(size=0.05) +
		theme_minimal() +
		xlab('UMAP 1') +
		ylab('UMAP 2') +
		ggtitle(title) +
		scale_colour_distiller(palette = palette, direction = 1, limits = c(lower.limit,upper.limit)) +
		theme(legend.position = 'right',
			plot.title = element_text(hjust=0.5),
			panel.border = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.line = element_line(colour = "black"))

	return ( plot )

}



#' Plotting heatmaps for post-binary classification clustering results
#'
#' This function plots a heatmap of clusters of a user-specified cell type. To be used after running cluster_subsets. Returns a graph.
#' @param data (Required) Data in data frame format. Expects results from load_data
#' @param cluster.results (Required) Clustering results in data frame format. Expects results from cluster_subsets
#' @param cell.type (Required) Cell subset to be plotted. This must be an entry from the binary classification specification 
#' @export
#' @examples
#' plot_cluster_heatmap(data, cluster.results, 'T Cells, CD4')
plot_cluster_heatmap <- function ( data = NULL,
	cluster.results = NULL,
	cell.type = NULL ) {

	if ( is.null(data) ) {
		return ('Please provide data. Use load_data')
	}

	if ( is.null(cluster.results) ) {
		return ('Please provide cluster results. Use cluster_subset')
	}


	if ( is.null(cell.type) ) {
		print('Please choose the cell type to plot.')
		return ()
	} else {

		colnames(cluster.results) <- c('Cell.Type', 'Cluster')
		to.plot <- cbind(data, cluster.results)
		medians <- to.plot %>%
			dplyr::filter(Cell.Type == cell.type) %>%
			select(-Cell.Type) %>%
			group_by(Cluster) %>%
			summarise_all(median) %>%
			select(-Cluster)
		medians <- as.matrix(medians)
		rownames(medians) <- seq(1, dim(medians)[1], 1)
		
		percentage <- to.plot %>%
			dplyr::filter(Cell.Type == cell.type) %>%
			select(-Cell.Type) %>%
			group_by(Cluster) %>%
			summarise(percentage=(100 * n() / dim(to.plot)[1]))

		myCol <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
		myBreaks <- seq(0, 3, length.out=100)

		colAnn <- HeatmapAnnotation("Cluster\nabundances (%)" = anno_barplot(
			as.matrix(percentage %>% select(percentage)),
			baseline = 0,
			border = FALSE,
			bar_width = 0.6,
			gp = gpar(fill = "black", fontsize = 16, fontface = "plain"),
			extend = 0.05,
			axis = TRUE,
			axis_param = list(gp = gpar(fontsize = 16, fontface = "plain")),
			width = NULL,
			height = unit(3, "cm")))

		hmap <- Heatmap(as.matrix(t(medians)),
			name = 'Expression',
			col = colorRamp2(myBreaks, myCol),

			heatmap_legend_param = list(
			color_bar = "continuous",
			legend_direction = "horizontal",
			legend_width = unit(8,"cm"),
			legend_height = unit(6,"cm"),
			title_position = "topcenter",
			title_gp = gpar(fontsize = 16, fontface = "bold"),
			labels_gp = gpar(fontsize = 16, fontface = "bold")),
	
			cluster_rows = TRUE,
			show_row_dend = TRUE,
			row_title = "Marker",
			row_title_side = "right",
			row_title_gp = gpar(fontsize = 16,  fontface="plain"),
			row_title_rot = 270,
			show_row_names = TRUE,
			row_names_gp = gpar(fontsize = 16, fontface="plain"),
			row_names_side = "right",
			row_dend_width = unit(25,"mm"),
	
			cluster_columns = TRUE,
			show_column_dend = TRUE,
			column_title = "Metacluster",
			column_title_side = "bottom",
			column_title_gp = gpar(fontsize = 16, fontface = "plain"),
			column_title_rot = 0,
			show_column_names = TRUE,
			column_names_gp = gpar(fontsize = 16, fontface = "plain"),
			column_dend_height = unit(25,"mm"),
	
			clustering_method_columns = "ward.D2",
			clustering_method_rows = "ward.D2",
	
			top_annotation = colAnn)

		plot <- draw(hmap, heatmap_legend_side="bottom", annotation_legend_side="top", row_sub_title_side="right")

		return ( plot )

	}

}
