labelCol <-
function(x, sam_cols, colorCodes) {
    if (is.leaf(x)) {
        ## fetch label
        label <- attr(x, "label")
        code <- substr(label, 1, 1)
        ## set label color to red for A and B, to blue otherwise
        attr(x, "nodePar") <- list(lab.col = colorCodes[sam_cols[match(label, names(sam_cols))]])
        
        attr(x, "edgePar") <- list(col = colorCodes[sam_cols[match(label, names(sam_cols))]])
        # names(sam_cols) is the sample name
    }
    return(x)
}
