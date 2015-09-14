# df <- read.csv("genewise-count.csv", row.names = 1)
# 
# axis.names <- as.list(colnames(df))


make_plot <- function(number_of_points){
    library("ggplot2")
    set.seed(5)
    x <- rnorm(number_of_points, 20, 5)
    
    set.seed(10)
    y <- rnorm(number_of_points, 50,50)
    
    df <- data.frame(x,y)
    
    ggplot(data=df, aes(x=x, y=y))+ 
        geom_point()+
        ggtitle("Awesome Plot")
}