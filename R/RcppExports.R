# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

lossCpp <- function(xy, lambda, radius, ED, ThreeD, ToleranceofLoss, maximumStep, ToleranceofStepsize, proportional, ALPHA, Bool) {
    .Call('vennplot_loop_R', PACKAGE = 'vennplot', xy, lambda, radius, ED, ThreeD, ToleranceofLoss, maximumStep, ToleranceofStepsize, proportional, ALPHA, Bool)
}

transCpp <- function(xy, radius, radiusvec, radiusall) {
    .Call('vennplot_trans_R', PACKAGE = 'vennplot', xy, radius, radiusvec, radiusall)
}

allDisjointCpp <- function(xy1, xy2, radius1, radius2, delta) {
    .Call('vennplot_alldis_R', PACKAGE = 'vennplot', xy1, xy2, radius1, radius2, delta)
}

closeCpp <- function(xy1, xy2, radius1, radius2, delta, direc) {
    .Call('vennplot_close_R', PACKAGE = 'vennplot', xy1, xy2, radius1, radius2, delta, direc)
}

binaryIndexCpp <- function(M, xy, radius, k, yuan, xuan, num) {
    .Call('vennplot_binaryIndexCpp', PACKAGE = 'vennplot', M, xy, radius, k, yuan, xuan, num)
}

goThroughPixelCpp <- function(myList, m, num) {
    .Call('vennplot_goThroughPixelCpp', PACKAGE = 'vennplot', myList, m, num)
}

countCpp <- function(M, Me) {
    .Call('vennplot_countCpp', PACKAGE = 'vennplot', M, Me)
}

getRidofZeroCpp <- function(M) {
    .Call('vennplot_getRidofZeroCpp', PACKAGE = 'vennplot', M)
}

binaryIndexThreeDCpp <- function(myList, xy, radius, k, yuan, xuan, zuan, num) {
    .Call('vennplot_binaryIndexThreeDCpp', PACKAGE = 'vennplot', myList, xy, radius, k, yuan, xuan, zuan, num)
}

goThroughPixelThreeDCpp <- function(list, m, num) {
    .Call('vennplot_goThroughPixelThreeDCpp', PACKAGE = 'vennplot', list, m, num)
}

allConnectedCpp <- function(xy, radius, ThreeD) {
    .Call('vennplot_allConnectedCpp', PACKAGE = 'vennplot', xy, radius, ThreeD)
}

distanceCpp <- function(r1, r2, theta1, theta2, S, ThreeD) {
    .Call('vennplot_distanceCpp', PACKAGE = 'vennplot', r1, r2, theta1, theta2, S, ThreeD)
}

BoolScaleNMCpp <- function(proportional, value, LAMBDA, STRESS) {
    .Call('vennplot_BoolScaleNMCpp', PACKAGE = 'vennplot', proportional, value, LAMBDA, STRESS)
}

BoolScaleLCpp <- function(proportional, value, stress_n, stress) {
    .Call('vennplot_BoolScaleLCpp', PACKAGE = 'vennplot', proportional, value, stress_n, stress)
}

BoolDistanceCpp <- function(proportional, value, f1, f2, thetanew, theta) {
    .Call('vennplot_BoolDistanceCpp', PACKAGE = 'vennplot', proportional, value, f1, f2, thetanew, theta)
}

