make_table <- function(number_of_points){
    
    set.seed(5)
    x <- rnorm(number_of_points, 20, 5)
    
    set.seed(10)
    y <- rnorm(number_of_points, 50, 50)
    
    df <- data.frame(x,y)
}
