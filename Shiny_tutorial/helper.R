# Generates a data frame of random numbers
make_df <- function(number_of_points){
  
  x <- rnorm(number_of_points, mean = 20, sd =  5)
  
  y <- rnorm(number_of_points, mean = 50, sd = 50)
  
  df <- data.frame(x,y)
  
  return(df)
}
















