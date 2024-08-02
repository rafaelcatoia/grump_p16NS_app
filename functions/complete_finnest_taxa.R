#Explanation:
#Iterate Over Rows and Columns: Starting from the second column, check if the entry is empty. If it is, replace it with the value of the previous column concatenated with the current column name.
#Convert Back to Matrix: Once the operation is complete, convert the data frame back to a matrix.
#df_aux = df_test

complete_finnest_taxa <- function(df_aux) {
  col_names = colnames(df_aux)
  for (i in 1:(ncol(df_aux)-1)) {
    incomplete_idx = df_aux[,i+1] == ""
      df_aux[incomplete_idx,i+1] <- paste(col_names[i],df_aux[incomplete_idx,i],sep='_')
    }
  # Convert back to matrix
  return(df_aux)
}

# Example usage
#mat <- matrix(c("A", "", "C", "D", "E", "", "", "H", "I"), nrow = 3, byrow = TRUE)
#
#colnames(mat) <- c("X1", "X2", "X3")
#df_test <- data.frame(mat)
#complete_finnest_taxa(df_test)
