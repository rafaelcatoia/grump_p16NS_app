#Explanation:
#Iterate Over Rows and Columns: Starting from the second column, check if the entry is empty. If it is, replace it with the value of the previous column concatenated with the current column name.
#Convert Back to Matrix: Once the operation is complete, convert the data frame back to a matrix.

complete_finnest_taxa <- function(df_aux) {
  for (i in 1:nrow(df_aux)) {
    for (j in 2:ncol(df_aux)) { # Start from second column
      if (df_aux[i, j] == "") {
        df_aux[i, j] <- paste(colnames(df_aux)[j],df_aux[i, j-1], sep = "_")
      }
    }
  }
  
  # Convert back to matrix
  return(df_aux)
}

# Example usage
# mat <- matrix(c("A", "", "C", "D", "E", "", "", "H", "I"), nrow = 3, byrow = TRUE)
# 
# colnames(mat) <- c("X1", "X2", "X3")
# print(mat)
# print(complete_matrix(data.frame(mat)))
