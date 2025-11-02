package main

import (
	"gonum.org/v1/gonum/stat"
	//"math"
)

//gonum/stat package includes func to calculate correlation matrix, linear regression, mean, standard dev, R^2 (coefficient of determination)

//Pearson Formula: r = numSamples*(summation of gene X * gene Y expressionVal for each gene in samples ) - (summation of expression values for gene X)*(summation of expression values for gene Y) / sqrt((numSamples* summation X^2 - (summation gene X)^2))*((numSamples* summation Y^2 - (summation gene Y)^2))

// Another way to explain is covariance of 2 genes expression vectors divided by product of their standard deviations
// CalculatePearsonCoefficient takes a slice of floats for 2 genes x and y and returns the a float64 corresponding to its pearson correlation coefficient
func CalculatePearsonCoefficient(X, Y []float64) float64 {

	//Calculate the covariance
	covariance := stat.Covariance(X, Y, nil)

	//calculate the standard deviation
	stdDevX := stat.StdDev(X, nil)
	stdDevY := stat.StdDev(Y, nil)

	pearsonCoefficient := covariance / stdDevX * stdDevY

	return pearsonCoefficient
}
