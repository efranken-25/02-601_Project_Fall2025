package main

import (
	"math"
	"sort"
	"strconv"

	"golang.org/x/exp/stats"
	"gonum.org/v1/gonum/stat"
)

// Function TransformMapToSlice takes a map's outerkey (geneName) and the map its in as input
// and returns a slice corresponding to the float64 tpm values for each of the samples for that gene
func TransformMapToSlice(geneName string, geneExpressionBreastMap map[string]map[string]float64) []float64 {
	//initialize an empty slice to for a geneName, the float64 values are the tpm values for each of the samples
	geneNameslice := make([]float64, 0)

	//lookup that input geneName key in the map. The nested structure returns an inner map corresponding to
	//the inner keys being the sample names and values are the tpm values
	innerMap := geneExpressionBreastMap[geneName]

	//range through the inner map, and create the slice of tpm values
	for _, tpmVal := range innerMap {
		geneNameslice = append(geneNameslice, tpmVal)
	}

	//return the slice (pseudo vector) for one gene
	return geneNameslice

}

// Store all of the sample names from the Data
// Function GetSampleNames takes a map[string]map[string]float64 as input and grabs the inner string key corresponding to the
// sample names and returns it as a slice of sample names
func GetSampleNames(geneExpressionBreastMap map[string]map[string]float64) []string {

	samples := make([]string, 0)

	for _, innerMap := range geneExpressionBreastMap {
		for sample := range innerMap {
			// check if it's already in the list
			found := false
			for _, s := range samples {
				if s == sample {
					found = true
					break
				}
			}
			if !found {
				samples = append(samples, sample)
			}
		}
	}
	return samples //return slice of string names
}

// Store all of the gene names from the Data
// Function GetGeneNames takes a map[string]map[string]float64 as input and grabs the outer string key corresponding to the
// gene names and returns it as a sorted slice of gene names
func GetGeneNames(geneExpressionBreastMap map[string]map[string]float64) []string {

	//initialize list to store geneNames
	geneNames := make([]string, 0)

	//range over the nested map and extract the outer key corresponding to gene names and add it to my list for storage
	for g := range geneExpressionBreastMap {
		geneNames = append(geneNames, g)
	}

	sort.Strings(geneNames) //sort the list of gene names in alphabetical order

	return geneNames
}

// Use Mean-Based Filtering in order to normalize data while maintaining fixed number of samples across each gene
// Function MeanBasedFilter takes map[string]map[string]float64 as input and returns a filtered map[string]map[string]float64 as output
// Step 1: Find/confirm number of samples for each gene (21 samples)
// Step 2: For each gene, Calculate the mean of the tpm values across its samples and then,
// Step 3: If the mean of that gene's tpm values is <=10 then remove the gene from further analyses
func MeanBasedFilter(geneExpressionBreastMap map[string]map[string]float64, samples []string, meanThresh float64) map[string]map[string]float64 {

	//initialize a nested map for return corresponding to the input map filtered for only genes with a mean tpm greater than 10.0
	out := make(map[string]map[string]float64, len(geneExpressionBreastMap))

	//find the number of samples
	n := len(samples)

	//range through the input nested map
	for gene, row := range geneExpressionBreastMap {
		// Sum the tpm values over all samples for a given gene
		var sum float64 //declare the sum within the loop so that it resets to 0 with each loop for every gene
		for _, s := range samples {
			sum += row[s] //sum the tpm values for all the samples for a given gene
		}
		//calculate the mean by dividing the tpm value sum by number of samples
		mean := sum / float64(n)

		// Keep gene only if mean TPM > threshold(meanThresh)
		if mean > meanThresh {
			cpy := make(map[string]float64, len(row))
			for s, v := range row {
				cpy[s] = v
			}
			out[gene] = cpy
		}
	}
	return out
}

// Sorting by sample number in order helper function
// sortSampleNames take a slice of strings as input and sorts the strings in order based on their # following "BRCA"
func sortSampleNames(samples []string) {
	sort.Slice(samples, func(i, j int) bool {
		getNum := func(s string) int {
			//assume names like BRCA12. first 4 letters, then number
			n, _ := strconv.Atoi(s[4:])
			return n
		}
		return getNum(samples[i]) < getNum(samples[j])
	})
}

// Function ComputePearsonCorrelation takes as input a sorted list of gene names and sample names, as well as a nested filtered map
// of map[genenames]map[samplenames]TPMValues and it computes the Pearson Correlation between each pair of genes across samples
// it returns a 2D matrix of pearson coeff values ([][]float64) with row and column indices corresponding to indices of genes found in sortedGeneNames list
func ComputePearsonCorrelation(sortedGeneNames, sortedSampleNames []string, filteredGeneExpressionBreastMap map[string]map[string]float64) [][]float64 {

	//initialize a square matrix for return
	n := len(sortedGeneNames)
	CorrMatrix := make([][]float64, n)
	//make the columns for a gene x gene square matrix
	for i := range CorrMatrix {
		CorrMatrix[i] = make([]float64, n)
	}

	//fill in the diagonal as 1.0 when comparing same gene to itself
	for k := 0; k < n; k++ {
		CorrMatrix[k][k] = 1.0
	}

	//range through the sorted list of gene names to compare every gene with every other gene for the correlation
	for i := 0; i < len(sortedGeneNames); i++ {
		for j := i + 1; j < len(sortedGeneNames); j++ {
			g1 := sortedGeneNames[i] //store the first gene for comparison
			g2 := sortedGeneNames[j] //store the second gene for comparison

			//Now lookup each gene's inner map (row) containing its samples and TPM values
			row1 := filteredGeneExpressionBreastMap[g1]
			row2 := filteredGeneExpressionBreastMap[g2]

			var gene1, gene2 []float64 //initialize slice of floats for each gene's TPM values as input to the correlation function call

			//range through the sorted sample names in order
			for _, s := range sortedSampleNames {
				//lookup that sample for each gene and get their associated TPM value
				gene1TPM, ok1 := row1[s]
				gene2TPM, ok2 := row2[s]

				//ensure that both genes' TPMvalues exist
				if ok1 && ok2 {
					gene1 = append(gene1, gene1TPM)
					gene2 = append(gene2, gene2TPM)
				}
			}

			if len(gene1) >= 3 { //must have enough samples for the correlation
				//compute the pearson correlation between different genes across samples
				r := stat.Correlation(gene1, gene2, nil)
				// store r correlation value in 2D matrix. Row and Column indices correspond to gene names
				// found by the same indices stored in sortedGeneNames list for reference
				// Diagonal should read with 1's
				CorrMatrix[i][j] = r
				CorrMatrix[j][i] = r
			} else { //otherwise indicate as not computed
				CorrMatrix[i][j] = math.NaN()
				CorrMatrix[j][i] = math.NaN()
			}
		}
	}
	return CorrMatrix
}

/*
// WriteCorrelationCSV writes the correlation matrix to a CSV file.
// The first row/column contain gene names as headers.
func WriteCorrelationCSV(filename string, sortedGeneNames []string, corrMatrix [][]float64) error {
	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	w := csv.NewWriter(f)
	defer w.Flush()

	n := len(sortedGeneNames)

	// Optional: sanity check
	if len(corrMatrix) != n {
		return fmt.Errorf("corrMatrix dimension (%d) does not match geneNames length (%d)", len(corrMatrix), n)
	}

	// ---- Write header row: "", gene1, gene2, ... ----
	header := make([]string, n+1)
	header[0] = "" // top-left corner empty
	for i, name := range sortedGeneNames {
		header[i+1] = name
	}
	if err := w.Write(header); err != nil {
		return err
	}

	// ---- Write each row: geneName, corr values... ----
	for i, row := range corrMatrix {
		if len(row) != n {
			return fmt.Errorf("row %d length (%d) does not match geneNames length (%d)", i, len(row), n)
		}

		record := make([]string, n+1)
		record[0] = sortedGeneNames[i]

		for j, v := range row {
			if math.IsNaN(v) {
				// choose how you want to represent missing correlations:
				record[j+1] = "" // or "NaN"
			} else {
				// format with, say, 4 decimal places
				record[j+1] = strconv.FormatFloat(v, 'f', 4, 64)
			}
		}

		if err := w.Write(record); err != nil {
			return err
		}
	}

	return nil
}*/

// calculate the quantile and choose the 95% as the cutoff corr value threshold
// need to go thru the matrix, append all the values into a growing slice
// sort the slice
// call the quantile operation from the stats package to compute the 90 95 99%
func TransformMatrixToSlice(BreastPearsonCorrelationMatrix [][]float64) []float64 {
	quantileSlice := make([]float64, 0)

	numRows := len(BreastPearsonCorrelationMatrix)
	numCols := len(BreastPearsonCorrelationMatrix[0])

	//only check upper triangle of the distance matrix for efficiency
	//and to ignore the diagonal "1"s
	for i := 0; i < numRows; i++ {
		for j := i + 1; j < numCols; j++ {
			value := BreastPearsonCorrelationMatrix[i][j]
			quantileSlice = append(quantileSlice, value)
		}
	}
	return quantileSlice
}

// Function SortCorrVals takes as input a slice of float64 values,
// it sorts the values in increasing numerical order and it returns the
// sorted []float64
func SortCorrVals(quantileSlice []float64) []float64 {
	sort.Float64s(quantileSlice)
	return quantileSlice
}

// Function ComputeQuantile takes a sorted slice of float64 values corresponding to the
// correlation values from the corr matrix, and it computes the 90th, 95th, and 99th quantile
// to be analyzed for determining the correlation threshold for edge determination
func ComputeQuantile(quantileSlice []float64) []float64 {

	quantileProbability := []float64{0.90, 0.95, 0.99}

	computedQuantiles := stats.Quantiles(quantileSlice, quantileProbability...)

	return computedQuantiles
} //prints: [0.5576559638595919 0.6474333018882147 0.779461567805639]

// Function BuildGraph takes a correlation matrix as input along with the threshold and sorted list of gene names
// it ranges over the values and determines
// if there is an edge between correlated genes based on the computed threshold and builds a graph network.
// If the absolute value of the correlation between 2 genes is >= 0.65, then there is an edge between those
// two genes (nodes)
func BuildGraph(BreastPearsonCorrelationMatrix [][]float64, breastGeneNames []string, threshold float64) GraphNetwork {

	numGenes := len(BreastPearsonCorrelationMatrix)
	numCols := len(BreastPearsonCorrelationMatrix[0])

	//input the nodes representing each gene in the breastGraph
	breastGraph := make(GraphNetwork, numGenes)
	for i := 0; i < numGenes; i++ {
		breastGraph[i] = &Node{
			ID:       i,
			GeneName: breastGeneNames[i],
			Edges:    []*Edge{},
		}
	}

	//range over the breast correlation matrix upper triangle
	for i := 0; i < numGenes; i++ {
		for j := i + 1; j < numCols; j++ {
			corr := BreastPearsonCorrelationMatrix[i][j]
			weight := math.Abs(corr)
			//check if the absolute value of the correlation is above or equal to the threshold to determine if there is an edge between them
			if weight >= threshold {
				//then each gene is a node and
				//there exists an edge between the 2 genes
				u := breastGraph[i]
				v := breastGraph[j]

				//create the edges
				e1 := &Edge{To: v, Weight: corr}
				e2 := &Edge{To: u, Weight: corr}

				//give the nodes those edges
				u.Edges = append(u.Edges, e1)
				v.Edges = append(v.Edges, e2)
			}

		}
	}

	return breastGraph
}

// Identify clustering using levain
