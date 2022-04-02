#' Running BinaryClust with default settings
#'
#' This function performs all major steps of BinaryClust using default settings. Returns a list of data frames.
#' @param input.file (Required) Path to input file
#' @param class.file (Required) Path to classification file in comma-separated format. Please see tutorial for details.
#' @export
#' @examples
#' BinaryClust.results <- run_BinaryClust('example.fcs', class.file = 'cell_types.csv')
run_BinaryClust <- function ( input.file = NULL,
	class.file = NULL ) {

	if ( is.null(input.file) ) {
		return ('Please provide a FCS file.')
	}

	if ( is.null(class.file) ) {
		print('Please provide a classification file.')
	}

	
	data <- load_data(input.file)
	binary.results <- binary_class(data, class.file = class.file)
	binary.summary <- binary_summary(data, binary.results)
	cluster.results <- cluster_subsets(data, binary.results)
	cluster.summary <- cluster_summary(data, cluster.results)

	to.return <- list(data, binary.summary, cluster.results, cluster.summary)
	
	return ( to.return )

}



cluster_suggestion <- function ( data,
	cluster.summary,
	subtype.file = NULL ) {

	if ( is.null(subtype.file) ) {
		print('Please provide a subpopulation classification file.')
	} else {

#		cluster.summary['Cell.Subtype'] <- ''
		
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
	
	
		temp.data <- subset(cluster.summary, select = -c(Cell.Type, Cluster, Frequency, Percentage))
		for ( temp.marker in rownames(kmeans.limits) ) {
			temp.data[temp.data[,temp.marker] < kmeans.limits[temp.marker,'pos_min'],temp.marker] <- '-'
			temp.data[temp.data[,temp.marker] >= kmeans.limits[temp.marker,'pos_min'] & temp.data[,temp.marker] < kmeans.limits[temp.marker,'Q1'],temp.marker] <- 'L'
			temp.data[temp.data[,temp.marker] >= kmeans.limits[temp.marker,'Q1'] & temp.data[,temp.marker] < kmeans.limits[temp.marker,'Q3'],temp.marker] <- 'I'
			temp.data[temp.data[,temp.marker] >= kmeans.limits[temp.marker,'Q3'],temp.marker] <- 'H'
		}
		temp.data[,'Number'] <- rownames(temp.data)

		subclass.cat <- read.table(subtype.file, header = TRUE, sep = ',', check.names = FALSE)
		subclass.cat[,colnames(subclass.cat) %in% colnames(data)] <- lapply(subclass.cat[,colnames(subclass.cat) %in% colnames(data)], function (x) { gsub(' ', '', x) })
		subclass.cat[subclass.cat == ''] <- 'A'

		for ( temp.type in unique(cluster.summary[,'Cell.Type']) ) {

			if ( temp.type %in% subclass.cat[,'Cell Type'] ) {
	
				for ( cell.subtype in as.vector(subclass.cat[subclass.cat[,'Cell Type'] == temp.type,'Cell Subtype']) ) {

					temp.row <- subclass.cat[subclass.cat[,'Cell Subtype'] == cell.subtype,]
					rownames(temp.row) <- temp.row[,'Cell Subtype']
					temp.row <- subset(temp.row, select = -c(`Cell Type`, `Cell Subtype`))
					temp.row <- temp.row[,temp.row != 'A', drop = FALSE]
	
					temp.row.fine <- temp.row[,temp.row != '+', drop = FALSE]
					temp.row.all <- temp.row[,temp.row == '+', drop = FALSE]
	
					rows.required <- merge(temp.data[cluster.summary[,'Cell.Type'] == temp.type,], temp.row.fine)
					if ( length(colnames(temp.row.all)) == 1 ) {
						rows.required <- rows.required[rows.required[,colnames(temp.row.all)] %in% c('L', 'I', 'H'),][,'Number']
					} else {
						rows.required <- rows.required[apply(rows.required[,colnames(temp.row.all)], 1, function (x) { all(x %in% c('L', 'I', 'H')) }),][,'Number']
					}


					cluster.summary[rows.required, 'Cell.Subtype'] <- cell.subtype
	
				}
##### BROKEN BELOW


				cluster.summary[cluster.summary[,'Cell.Type'] == temp.type & cluster.summary[,'Cell.Subtype'] == '','Cell.Subtype'] <- paste0(temp.type, ', Unclassified')

			} else {
	
				cluster.summary[cluster.summary[,'Cell.Type'] == temp.type & cluster.summary[,'Cell.Subtype'] == '','Cell.Subtype'] <- temp.type

			}

##### EDIT BELOW


#			clustering.results[clustering.results[,'Cell Type'] == temp.type,'Cell Subtype'] <- temp.data[clustering.results[clustering.results[,'Cell Type'] == temp.type,'Cluster'],'Cell.Subtype']


		}

#		return ( clustering.results )
		return ( cluster.summary )
		
	
	}
}


