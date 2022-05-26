fix_negatives <- function(fix_data, column_to_fix) {
  #########################################################
  # Input: fix_data(dataframe) -> dataframe to fix        #
  #        column_to_fix(string) -> column to fix         #
  # Purpose: Assumes that the column to fix is some 'new' #
  #  column with an available 'total' column. Fixes all   #
  #  negative values in the given column.                 #
  # Output: (dataframe) res -> data formatted into form   #
  #                                                       #
  #########################################################
  
  # find coulmns that need to be fixed and assign fix ratios
  tofix <- filter(fix_data, !!as.symbol(column_to_fix)<0)
  tofix$ratio <- 1+dplyr::select(tofix, !!as.symbol(column_to_fix))/(dplyr::select(tofix, !!as.symbol(gsub("new", "total", column_to_fix, fixed = TRUE)))-dplyr::select(tofix, as.symbol(column_to_fix)))
  
  for( i in 1:nrow(tofix)){ #apply fix to all previous days
    # retrieve the column to fix
    colval <- tofix[i,]
    # extract data from the column
    startDate <- colval$date
    isoc <- colval$iso_code
    ra <- colval$ratio
    # find first and last index to fix
    s_ind <- head(which(fix_data$iso_code==isoc), 1)
    e_ind <- which(fix_data$iso_code==isoc & fix_data$date==startDate)
    # check the index
    # fix data
    fix_data[(s_ind:e_ind-1),column_to_fix] <- fix_data[(s_ind:e_ind-1),column_to_fix]*as.double(ra)
    fix_data[e_ind,column_to_fix] <- mean(fix_data[((e_ind-7):(e_ind-1)),column_to_fix])
  } 
  return(fix_data)
}