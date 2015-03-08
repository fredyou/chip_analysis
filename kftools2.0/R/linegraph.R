linegraph <-
function(data = "Book1.txt", xlab = "xaxis", ylab = "yaxis", legendTitle = "Samples") {
    require(ggplot2)
    require(reshape2)
    if (class(data) == "data.frame") {
        data = data
    } else if (class(data) == "character") {
        data = read.delim(data, quote = "", stringsAsFactors = F, check.names = F)
    }
    
    data = setNames(data.frame(t(data[, -1])), data[, 1])  #set colnames
    data$num = rownames(data)
    data.m = melt(data)
    
    # col=sample(colors()[-1], ncol(data), replace = FALSE)
    colours = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", 
        "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", 
        "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")
    if (dim(data)[1] > 26) {
        colours = colorRampPalette(c("red", "orange", "green", "purple", "blue"))(dim(data)[1])
    }
    # colours=data.frame(colour=colours[sample(nrow(colours),dim(data)[1]),],group=row.names(data))
    
    colours = data.frame(colour = sample(colours, dim(data)[1]), group = row.names(data))
    
    row.names(colours) = colours[, 2]
    data.m$colours = colours[data.m$num, 1]
    colours = data.frame(colour = colours[sample(nrow(colours), dim(data)[1]), ], group = row.names(data))
    # colours=data.frame(colour=sample(colours,dim(data)[1]),group=row.names(data))
    row.names(colours) = colours[, 2]
    data.m$colours = colours[data.m$num, 1]
    ggplot(data = data.m, aes(x = variable, y = value, group = num, colour = colours)) + geom_line() + geom_point() + 
        xlab(xlab) + ylab(ylab) + scale_color_discrete(name = legendTitle, breaks = data.m$colours, labels = data.m$num)
    ggsave(filename = "result_line.png", plot = last_plot(), height = par("din")[2] - 3)
    
}
