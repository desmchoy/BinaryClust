#' Examining an FCS File
#'
#' This function reads an FCS file and plots density plots for individual markers.
#' @param input.file (Required) Path to input file
#' @param transform 'arcsinh' or 'none'. 'arcsinh' performs an arcsinh transformation. 'none' performs no data trasformation
#' @param cofactor Cofactor for arcsinh transformation (Default: 5)
#' @export
#' @examples
#' examine_data('example.fcs')
examine_data <- function ( input.file = NULL,
	transform = 'arcsinh',
	cofactor = 5 ) {

	if ( is.null(input.file) ) {
		return ('Please provide a FCS file.')
	}

	flowframe <- read.FCS(input.file)
	raw.data <- data.frame(exprs(flowframe))
	desc <- flowframe@parameters@data$desc
	desc <- desc[grep('_', desc)]
	markers <- desc
	columns.to.use <- which(flowframe@parameters@data$desc %in% markers)
	useable.raw.data <- raw.data[columns.to.use]

	if ( transform == 'arcsinh' ) {

		transformed.data <- asinh(useable.raw.data / cofactor)

	} else if ( transform == 'none' ) {

		transformed.data <- useable.raw.data

	} else {

		print('Custom transformation currently unsupported. Transforming using default method.')
		transformed.data <- asinh(useable.raw.data / 5)

	}

	colnames(transformed.data) <- markers
	transformed.data[,'Cell'] <- c(1:dim(transformed.data)[1])
	to.return <- transformed.data %>%
		pivot_longer(-Cell, names_to = 'Marker', values_to = 'Expression') %>%
		ggplot(aes(x=Expression, fill=Marker)) +
			geom_density(fill = 'black') +
			xlab('Transformed Read Count') +
			ylab('Normalised Density') +
			ggtitle(paste0('Distributions of Marker Expressions')) +
			coord_cartesian(ylim=c(0,2)) +
			theme_minimal() +
			theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				plot.title = element_text(face = 'bold', hjust = 0.5)) +
			facet_wrap(vars(Marker))

		
	return ( to.return )

}



#' Printing the parameters of an FCS file
#'
#' This function reads an FCS file and prints out the parameter page of the flowCore flowFrame object. It is intended to aid feature selection.
#' @param input.file (Required) Path to input file
#' @export
#' @examples
#' print_parameters('example.fcs')
print_parameters <- function ( input.file = NULL ) {

	if ( is.null(input.file) ) {
		print('Please provide a FCS file.')
	} else {

		flowframe <- read.FCS(input.file)
		to.return <- flowframe@parameters@data

		return ( to.return )

	}

}



#' Loading an FCS File
#'
#' This function reads an FCS file and returns a data frame of read counts
#' @param input.file (Required) Path to input file
#' @param channels Default: NULL. Specify to subset data. Please provide a vector with full channel name (e.g. c('142Nd_CD19', '154Sm_CD3')). Channel names can be found using the print_parameters function.
#' @param transform 'arcsinh' or 'none'. 'arcsinh' performs an arcsinh transformation. 'none' performs no data trasformation
#' @param cofactor Cofactor for arcsinh transformation (Default: 5)
#' @export
#' @examples
#' data <- load_data('example.fcs')
load_data <- function ( input.file = NULL,
	channels = NULL,
	transform = 'arcsinh',
	cofactor = 5 ) {

	if ( is.null(input.file) ) {
		return ('Please provide a FCS file.')
	}

	flowframe <- read.FCS(input.file)
	raw.data <- data.frame(exprs(flowframe))
	desc <- flowframe@parameters@data$desc

	if ( is.null(channels) ) {

		desc <- desc[grep('_', desc)]
		markers <- desc[!grepl('[Ee]vent_[Ll]ength|[Bb]ead|DNA|[Ll]ive|[Dd]ead|[Vv]iability', desc)]
		markers.names <- sapply(markers, function (x) gsub('^[^_]*_(*)', '\\1', x))
						
	} else {
			
		markers <- desc[desc %in% channels]
		markers.names <- markers			

	}

	columns.to.use <- which(flowframe@parameters@data$desc %in% markers)
	useable.raw.data <- raw.data[columns.to.use]
	colnames(useable.raw.data) <- markers.names

	if ( transform == 'arcsinh' ) {

		transformed.data <- asinh(useable.raw.data / cofactor)

	} else if ( transform == 'none' ) {

		transformed.data <- useable.raw.data

	} else {

		print('Custom transformation currently unsupported. Transforming using default method.')
		transformed.data <- asinh(useable.raw.data / cofactor)

	}


	to.return <- transformed.data
	return ( to.return )
	
}
