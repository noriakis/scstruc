
#' shdmat
shdmat <- function(netlist) {
    sapply(netlist, function(x) sapply(netlist, function(y) bnlearn::shd(x,y)))
}
