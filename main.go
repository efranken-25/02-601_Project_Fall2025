package main

import (
	"fmt"
	"log"
)

func main() {

	//Begin Analyzing Breast Cancer RNA seq gene count Data
	// Create a map of maps. Nested map mapping gene_name (string) to sampleID(string) to tpm value
	geneExpressionBreastMap, err1 := ReadGeneExpressionDirToGeneMap("/Project_Data/BreastData")

	// log any errors that may arise from the parsing/ reading in the breast data files
	if err1 != nil {
		log.Fatal(err1)
	}

	fmt.Println("Breast Cancer Gene Expression Data read and map created.")

	//Begin Analyzing Ovarian Cancer RNA seq gene count Data
	//Create a map of maps. Nested map mapping gene_name (string) to sampleID(string) to tpm value
	geneExpressionOvarianMap, err2 := ReadGeneExpressionDirToGeneMap("/Project_Data/OvarianData")

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

	fmt.Println("Gathering sample names for breast cancer data.")
	// Store the sample names for the Breast Cancer Data
	breastSampleNames := GetSampleNames(geneExpressionBreastMap)

	//Sort the breast sample names
	sortSampleNames(breastSampleNames)

	//Print the breast sample names to ensure they are sorted by numerical order
	fmt.Println(breastSampleNames)

	fmt.Println("Gathering sample names for ovarian cancer data.")
	// Store the sample names for the Ovarian Cancer Data
	ovarianSampleNames := GetSampleNames(geneExpressionOvarianMap)

	//Sort the ovarian sample names
	sortSampleNames(ovarianSampleNames)

	//Print the ovarian sample names to ensure they are sorted by numerical order
	fmt.Println(ovarianSampleNames)

	fmt.Println("Cleaning Data...Filtering out low expression genes...")
	// Use Mean-Based Filtering in order to normalize data while maintaining fixed number of samples across each gene
	filteredGeneExpressionBreastMap := MeanBasedFilter(geneExpressionBreastMap, breastSampleNames, 10.0)

	filteredGeneExpressionOvarianMap := MeanBasedFilter(geneExpressionOvarianMap, ovarianSampleNames, 10.0)

	fmt.Println("Gathering gene names for breast cancer data.")
	// Store the gene names for the Breast Cancer Data
	breastGeneNames := GetGeneNames(filteredGeneExpressionBreastMap)

	//Print the breast gene names to check the names and ensure they are sorted
	fmt.Println(breastGeneNames)

	fmt.Println("Gathering gene names for ovarian cancer data.")
	// Store the gene names for the Breast Cancer Data
	ovarianGeneNames := GetGeneNames(filteredGeneExpressionOvarianMap)

	//Print the breast gene names to check the names and ensure they are sorted
	fmt.Println(ovarianGeneNames)

	//Compute the Pearson Correlation
	fmt.Println("Computing Pearson Correlation between genes across all breast cancer samples.")

	PearsonCorrelationMatrix := ComputePearsonCorrelation(breastGeneNames, breastSampleNames, filteredGeneExpressionBreastMap)

	//temporary Checking the Pearson Correlation Matrix
	fmt.Println(PearsonCorrelationMatrix)
}
