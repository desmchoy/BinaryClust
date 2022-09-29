#' Loading BinaryClust results of multiple datasets
#'
#' This function concatenates multiple results of cluster_summary for differential analysis. Returns a data frame.
#' @param file.list (Required) A list of files in comma-separated format. Please see tutorial for details.
#' @export
#' @examples
#' multi.data <- load_multi_samples(file.list)
load_multi_samples <- function ( file.list = NULL ) {

	if ( is.null(file.list) ) {
		return ('Please provide a list of files.')
	}

	files <- read.csv(file.list, sep = ',', header = TRUE)
	conditions <- as.vector(unique(files[,'Condition']))
	if ( length(conditions) != 2 ) {
		return ('Please specify two (and only two) conditions.')
	}
	condition.1 <- conditions[1]
	condition.2 <- conditions[2]
	files.1 <- as.vector(unique(files[files[,'Condition'] == condition.1,'File']))
	files.2 <- as.vector(unique(files[files[,'Condition'] == condition.2,'File']))

	data.1 <- import_results(files, condition.1)
	data.2 <- import_results(files, condition.2)
	multi.data <- rbind(data.1, data.2)

	return ( multi.data )

}



import_results <- function ( files, condition ) {
	data <- lapply(as.vector(unique(files[files[,'Condition'] == condition,'File'])), read.csv, sep = ',', header = TRUE, row.names = 1)
	for ( i in c(1:dim(files[files[,'Condition'] == condition,])[1]) ) {
		data[[i]][,'Sample'] <- paste0(condition, ' Sample ', i)
		data[[i]][,'Condition'] <- condition
	}
	data <- do.call(rbind, data)
	return (data)
}



#' Plotting abundances of multiple datasets
#'
#' This function plots the abundance of multiple datasets. Returns a graph.
#' @param multi.data (Required) A data frame of clustering results. Expects results from load_multi_samples
#' @param column Choose between plotting Cell.Type or Cell.Subtype (Default: Cell.Type)
#' @param remove.unclass Option to remove unclassified cells from plot (Default: FALSE)
#' @param style Choose between two styles of plots, bar or dot (Default: bar)
#' @export
#' @examples
#' plot_diff_abundances(multi.data)
#' plot_diff_abundances(multi.data, column = 'Cell.Subtype', style = 'dot')
plot_diff_abundances <- function ( multi.data = NULL,
	column = 'Cell.Type',
	remove.unclass = FALSE,
	style = 'bar' ) {

	if ( is.null(multi.data) ) {
		return ('Please run load_multi_samples and provide the results.')
	}


	if ( remove.unclass == FALSE ) {
		multi.data <- multi.data
	} else if ( remove.unclass == TRUE ) {
		multi.data <- multi.data[multi.data[,'Cell.Type'] != 'Unclassified',]
	} else {
		return ('Please choose a valid option for remove.unclass - TRUE or FALSE')
	}

	if ( column == 'Cell.Type' ) {
		legend.label <- 'Cell Type'
	} else if ( column == 'Cell.Subtype' ) {
		legend.label <- 'Cell Subtype'
	} else {
		return ('Please choose a valid option for column - Cell.Type or Cell.Subtype')
	}


	if ( style == 'bar' ) {
		plot <- multi.data %>% group_by(Condition, Sample, .data[[column]]) %>%
			summarise(PertSum = sum(Percentage)) %>%
			ggplot(aes(x = Sample, y = PertSum, fill = .data[[column]])) +
			geom_col() +
			xlab('Sample') +
			ylab('Percentage Abundances') +
			labs(fill = legend.label) +
			theme_minimal() +
			theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.text.x = element_text(angle = 270, hjust = 1, vjust = 0.5),
				axis.line = element_line(colour = "black"),
				plot.title = element_text(face = "bold", hjust = 0.5)) +
			facet_grid(cols = vars(Condition), scales = 'free', space = 'free')
	} else if ( style == 'dot' ) {
		plot <- multi.data %>% group_by(Condition, Sample, .data[[column]]) %>%
			summarise(PertSum = sum(Percentage)) %>%
			ggplot(aes(x = Condition, y = PertSum, color = Condition)) +
			geom_jitter(width = 0.1) +
			ylab('Percentage Abundances') +
			theme_minimal() +
			theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				legend.position = 'none',
				axis.line = element_line(colour = "black"),
				plot.title = element_text(face = "bold", hjust = 0.5)) +
			facet_wrap(as.formula(paste("~", column)))
	} else {
		return ('Please choose a valid option for style - bar or dot')
	}


	return ( plot )

}

	

#' Performing differential analysis on multiple datasets
#'
#' This function performs differential analysis on datasets as imported by load_multiple_samples. Returns a list of three data frames: [[1]] percentage abundances of each cell (sub)type in each sample, [[2]] log2 fold changes of each cell (sub)type between two conditions and [[3]] p-values for the fold changes
#' @param multi.data (Required) A data frame of clustering results. Expects results from load_multi_samples
#' @param column Choose between plotting Cell.Type or Cell.Subtype (Default: Cell.Type)
#' @param remove.unclass Option to remove unclassified cells from plot (Default: FALSE)
#' @param transform Transformation method used for original analysis. 'default' assumes arcsinh transformation. 'none' indicates no transformation (Default: default)
#' @export
#' @examples
#' diff.results <- calculate_diff(multi.data)
#' diff.results <- calculate_diff(multi.data, transform = 'default', column = 'Cell.Subtype')
calculate_diff <- function ( multi.data = NULL,
	column = 'Cell.Type',
	remove.unclass = FALSE,
	transform = 'default' ) {

	if ( is.null(multi.data) ) {
		return ('Please run load_multi_samples and provide the results.')
	}

	conditions <- as.vector(unique(multi.data[,'Condition']))
	if ( length(conditions) != 2 ) {
		return ('Please specify two (and only two) conditions.')
	}
	condition.1 <- conditions[1]
	condition.2 <- conditions[2]


	if ( remove.unclass == FALSE ) {
		multi.data <- multi.data
	} else if ( remove.unclass == TRUE ) {
		multi.data <- multi.data[multi.data[,'Cell.Type'] != 'Unclassified',]
	} else {
		return ('Please choose a valid option for remove.unclass - TRUE or FALSE')
	}

	if ( column == 'Cell.Type' ) {
		neg <- 'Cell.Subtype'
		legend.label <- 'Cell Type'
	} else if ( column == 'Cell.Subtype' ) {
		neg <- 'Cell.Type'
		legend.label <- 'Cell Subtype'
	} else {
		return ('Please choose a valid option for column - Cell.Type or Cell.Subtype')
	}

	if ( transform == 'none' ) {
		to.process <- multi.data 
	} else if ( transform == 'default' ) {
		to.process <- multi.data %>% mutate_if(is.numeric, function (x) { 5 * sinh (x) })
	} else {
		return ('Please specify a valid transformation method - default or none')
	}


	abundances <- multi.data %>% group_by(Condition, Sample, .data[[column]]) %>%
		summarise(PertSum = sum(Percentage))
	abundances <- data.frame(abundances)


	log2.fold.change <- to.process %>%
		select(!c(Cluster, Frequency, Percentage, .data[[neg]])) %>%
		group_by(Condition, Sample, .data[[column]]) %>%
		summarise_all(mean) %>%
		ungroup() %>%
		select(!c(Sample)) %>%
		group_by(Condition, .data[[column]]) %>%
		summarise_all(mean) %>%
		group_by(.data[[column]]) %>%
		summarise_if(is.numeric, vars(log2(.x[Condition == condition.2] / .x[Condition == condition.1])))
	log2.fold.change <- data.frame(log2.fold.change)
	
	p.values.adj <- to.process %>%
		dplyr::filter(.data[[column]] %in% log2.fold.change[,column]) %>% 
		select(!c(Cluster, Frequency, Percentage, .data[[neg]])) %>%
		group_by(Condition, Sample, .data[[column]]) %>%
		summarise_all(mean) %>%
		ungroup() %>%
		select(!c(Sample)) %>%
		group_by(.data[[column]]) %>%
		summarise_if(is.numeric, vars(wilcox.test(.x[Condition == condition.1], .x[Condition == condition.2])$p.value)) 
	cell.column <- as.data.frame(p.values.adj %>% select(.data[[column]]))
	p.values.adj <- p.values.adj %>%
		summarise_if(is.numeric, vars(p.adjust(.x, method = 'fdr')))

	p.values.adj <- data.frame(p.values.adj)
	p.values.adj[is.na(p.values.adj)] <- 1
	p.values.adj <- cbind(cell.column, p.values.adj)

	
	to.return <- list(abundances, log2.fold.change, p.values.adj)

	return ( to.return )		

}



#' Plotting differential analysis results
#'
#' This function plots the results of differential analysis and returns a summary heatmap
#' @param diff.results (Required) A list of data frames from differential analysis. Expects results from calculate_diff
#' @param plot_p_values Option to display significance values on the heatmap (Default: TRUE)
#' @export
#' @examples
#' plot_diff_exp(diff.results)
plot_diff_exp <- function ( diff.results = NULL,
	plot_p_values = TRUE ) {


	if ( is.null(diff.results) ) {
		return ('Please run calculate_diff and provide the results.')
	}

	conditions <- as.vector(unique(diff.results[[1]][,'Condition']))
	condition.1 <- conditions[1]
	condition.2 <- conditions[2]
	column <- colnames(diff.results[[2]])[1]

	abundances <- diff.results[[1]] %>%
		dplyr::filter(.data[[column]] %in% diff.results[[2]][,column]) %>%
		group_by(Condition, .data[[column]]) %>%
		summarise(mean.freq = mean(PertSum))
	abundances.1 <- abundances %>%
		dplyr::filter(Condition == condition.1) %>%
		ungroup() %>%
		select(!c(Condition))
	abundances.1 <- data.frame(abundances.1)
	colnames(abundances.1)[2] <- 'mean.freq.1'

	abundances.2 <- abundances %>%
		dplyr::filter(Condition == condition.2) %>%
		ungroup() %>%
		select(!c(Condition))
	abundances.2 <- data.frame(abundances.2)
	colnames(abundances.2)[2] <- 'mean.freq.2'

	abundance.p <- diff.results[[1]] %>%
		dplyr::filter(.data[[column]] %in% diff.results[[2]][,column]) %>%
		group_by(.data[[column]]) %>%
		summarise(p.value = wilcox.test(PertSum[Condition == condition.1], PertSum[Condition == condition.2])$p.value) %>%
		mutate(p.value.adj = p.adjust(p.value, method = 'fdr'))
	abundance.p <- data.frame(abundance.p)

	is.na(diff.results[[2]]) <- sapply(diff.results[[2]], is.infinite)
	diff.results[[2]][is.na(diff.results[[2]])] <- 0

	summary.table <- merge(diff.results[[2]], abundances.1, by = column)
	summary.table <- merge(summary.table, abundances.2, by = column)
	summary.table <- merge(summary.table, abundance.p, by = column)
	rownames(summary.table) <- summary.table[,column]
	summary.table <- summary.table[,-1]


	dummy.background <- colorRamp2(c(0, 100), c("white", "white"))
	pch.function <- function (x) {
		if ( -log10(x) < 1.30 ) {
			pch = 6
		} else if ( -log10(x) >= 1.30 & -log10(x) < 2 ) {
			pch = 1
		} else if ( -log10(x) >= 2 & -log10(x) < 3 ) {
			pch = 3
		} else if ( -log10(x) >= 3 & -log10(x) < 4 ) {
			pch = 4
		} else if ( -log10(x) >= 4 ) {
			pch = 7
		}
		return (pch)
	}
#	p.value.pch <- sapply(summary.table[,'p.value.adj'], pch.function)
	p.value.pch <- sapply(summary.table[,'p.value'], pch.function)

	row.anno <- rowAnnotation('Condition 1' = anno_simple(summary.table[,'mean.freq.1'],
			pch = 20,
			pt_gp = gpar(col = "grey"),
			pt_size = unit(summary.table[,'mean.freq.1'] / 2, "mm"),
			col = dummy.background),
		'Condition 2' = anno_simple(summary.table[,'mean.freq.2'],
			pch = 20,
			pt_gp = gpar(col = "grey"),
			pt_size = unit(summary.table[,'mean.freq.2'] / 2, "mm"),
			col = dummy.background),
		'p-value' = anno_simple(summary.table[,'p.value.adj'],
			pch = p.value.pch,
			pt_size = unit(3, "mm"),
			col = dummy.background),
		' ' = anno_simple(rep(0, length(p.value.pch)),
			col = dummy.background),
		simple_anno_size_adjust = TRUE,
		width = unit(3, 'cm'),
		gap = unit(5, "mm"),
		annotation_name_side = "top")
	names(row.anno) <- c(condition.1, condition.2, 'p-value', '')

	if ( plot_p_values == FALSE ) {
		heatmap <- Heatmap(as.matrix(summary.table[,1:(dim(summary.table)[2]-4)]),
			name = 'Fold Change',
			row_names_side = 'left',
			cluster_rows = FALSE,
			cluster_columns = FALSE,
			right_annotation = row.anno,
			heatmap_legend_param = list(direction = "horizontal"),
			width = unit(30, "cm"), height = unit(20, "cm"))
	} else if ( plot_p_values == TRUE ) {
		heatmap <- Heatmap(as.matrix(summary.table[,1:(dim(summary.table)[2]-4)]),
			name = 'Fold Change',
			row_names_side = 'left',
			cluster_rows = FALSE,
			cluster_columns = FALSE,
			right_annotation = row.anno,
			heatmap_legend_param = list(direction = "horizontal"),
			width = unit(30, "cm"),
			height = unit(20, "cm"),
			cell_fun = function(j, i, x, y, width, height, fill) { if (as.matrix(diff.results[[3]][,-1])[i, j] < 0.05) { grid.text('S', x = x, y = y, gp = gpar(fontsize = 10)) } })
	} else {
		return ('Please enter either TRUE or FALSE for plot_p_values')
	}

	freq.legend <- Legend(title = 'Mean Abundance',
		labels = c('  1%', '  5%', '  10%'),
		type = 'points',
		pch = 19,
		size = unit(c(1, 5, 10) / 2, "mm"),
		legend_gp = gpar(col = "grey"),
		background = "white",
		grid_height = unit(0.8, "cm"))

	p.value.legend <- Legend(title = 'p-value',
		labels = c('  N.S.', '  < 0.05', '  < 0.01', '  < 0.001', '  < 0.0001'),
		type = 'points',
		pch = c(6, 1, 3, 4, 7),
		size = unit(3, "mm"),
		legend_gp = gpar(col = "black"),
			background = "white",
			grid_height = unit(0.8, "cm"))

	anno.legend <- packLegend(freq.legend, p.value.legend, row_gap = unit(1, "cm"))


	plot <- draw(heatmap, annotation_legend_list = anno.legend, heatmap_legend_side = 'bottom', annotation_legend_side = 'right')

	return ( plot )

}


