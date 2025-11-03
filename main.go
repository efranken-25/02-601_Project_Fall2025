package main

import (
	"fmt"
	"log"
)

func main() {

	//Begin Analyzing Breast Cancer RNA seq gene count Data
	// Create a map of maps. Nested map mapping gene_name (string) to sampleID(string) to tpm value
	geneExpressionBreastMap, err1 := ReadGeneExpressionDirToGeneMap("/Users/bethany/go/src/ProgrammingProject/GeneExpressionData/BreastData")

	// log any errors that may arise from the parsing/ reading in the breast data files
	if err1 != nil {
		log.Fatal(err1)
	}

	fmt.Println("Breast Cancer Gene Expression Data read and map created.")

	//Begin Analyzing Ovarian Cancer RNA seq gene count Data
	//Create a map of maps. Nested map mapping gene_name (string) to sampleID(string) to tpm value
	geneExpressionOvarianMap, err2 := ReadGeneExpressionDirToGeneMap("/Users/bethany/go/src/ProgrammingProject/GeneExpressionData/OvarianData")

	// log any errors that may arise from the parsing/ reading in the ovarian data files
	if err2 != nil {
		log.Fatal(err2)
	}

	fmt.Println("Ovarian Cancer Gene Expression Data read and map created.")

	fmt.Println("Printing first line of geneExpressionBreastMap")
	//Check the format of the data file by printing the first line of our geneExpressionBreastMap
	for key, val := range geneExpressionBreastMap {
		fmt.Println(key, val)
		break
	}

	fmt.Println("Printing first line of geneExpressionOvarianMap")
	//Check the format of the data file by printing the first line of our geneExpressionOvarianMap
	for key, val := range geneExpressionOvarianMap {
		fmt.Println(key, val)
		break
	}

	// Store the sample names for the Breast Cancer Data
	sampleNames := GetSampleNames(geneExpressionBreastMap)

	// Use Mean-Based Filtering in order to normalize data while maintaining fixed number of samples across each gene
	filteredGeneExpressionBreastMap := MeanBasedFilter(geneExpressionBreastMap, sampleNames, 10.0)

	filteredGeneExpressionOvarianMap := MeanBasedFilter(geneExpressionOvarianMap, sampleNames, 10.0)

	/*stat.Correlation(x, y, weights []float64) float64
	Computes the Pearson correlation coefficient between two data slices, x and y.
	*/

}
