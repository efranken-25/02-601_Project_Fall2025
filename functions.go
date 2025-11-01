package main

/*
import (
	"fmt"
	"gonum/stat" //will prob need for Pearson
	"math"
	"strconv"
)
*/

//gonum/stat package includes func to calculate correlation matrix, linear regression, mean, standard dev, R^2 (coefficient of determination)

//Pearson Formula: r = numSamples*(summation of gene X * gene Y expressionVal for each gene in samples ) - (summation of expression values for gene X)*(summation of expression values for gene Y) / sqrt((numSamples* summation X^2 - (summation gene X)^2))*((numSamples* summation Y^2 - (summation gene Y)^2))

//Another way to explain is covariance of 2 genes expression vectors divided by product of their standard deviations

//FilterLowExpression takes matrix, row, col in as input. if tpm levels are < 10, it will remove these genes. Returns new matrix without low expression data.
