package main

import (
	"02-601_Project_Fall2025/louvain"
	"context"
	"log"
	"math"
	"sort"
	"strconv"

	fisher "github.com/glycerine/golang-fisher-exact"
	"golang.org/x/exp/stats"
	"gonum.org/v1/gonum/stat"
)

// Function TransformMapToSlice takes a map's outerkey (geneName) and the map its in as input
// and returns a slice corresponding to the float64 tpm values for each of the samples for that gene
func TransformMapToSlice(geneName string, geneExpressionBreastMap ExpressionMap) []float64 {
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
func GetSampleNames(geneExpressionBreastMap ExpressionMap) []string {

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
func GetGeneNames(geneExpressionBreastMap ExpressionMap) []string {

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
func MeanBasedFilter(geneExpressionBreastMap ExpressionMap, samples []string, meanThresh float64) ExpressionMap {

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
func ComputePearsonCorrelation(sortedGeneNames, sortedSampleNames []string, filteredGeneExpressionBreastMap ExpressionMap) CorrelationMatrix {

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

// Plan
// calculate the quantile and choose the 95% as the cutoff corr value threshold
// need to go through the matrix, append all the values into a growing slice
// sort the slice
// call the quantile operation from the stats package to compute the 90 95 99% quantile probabilities

// Function TransformMatrixToSlice takes a symmetrical square matrix as input and flattens
// it into a single slice of float64 values for return.
func TransformMatrixToSlice(BreastPearsonCorrelationMatrix CorrelationMatrix) []float64 {
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
// it ranges over the values and determines if there is an edge between correlated genes based on the computed
// threshold and builds an undirected, weighted graph network.
// If the absolute value of the correlation between 2 genes is >= threshold, then there is an edge between those
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

// Identify clustering using Louvain algorithm from package by building Louvain graph type from the correlation matrix and threshold inputs.
/*
func BuildLouvainGraph is a helper function that converts the correlation matrix into a pointer to Louvain graph type
that is required as input for the Louvain algorithm.

This function also uses 2 calls from the package library:
	func NewGraph(numNodes int) *Graph
	func (g *Graph) AddEdge(u, v int, weight float64) error

This package requires a weighted undirected graph with this structure:
type Graph struct {
	NumNodes    int
	Adjacency   [][]int       // adjacency[i] = list of neighbors of node i
	Weights     [][]float64   // weights[i][j] = weight of edge from node i to neighbor adjacency[i][j]
	Degrees     []float64     // degrees[i] = weighted degree of node i
	TotalWeight float64       // sum of all edge weights
}
*/
func BuildLouvainGraph(BreastPearsonCorrelationMatrix [][]float64, threshold float64) *louvain.Graph {
	numGenes := len(BreastPearsonCorrelationMatrix) //equivalent to numNodes
	graph := louvain.NewGraph(numGenes)

	//range over the expression matrix
	for i := 0; i < numGenes; i++ {
		for j := i + 1; j < numGenes; j++ {
			//grabs the absolute value of all the edge weights (corr between genes)
			//then each gene is a node and
			//there exists an edge between the 2 genes
			weight := math.Abs(BreastPearsonCorrelationMatrix[i][j])
			if weight >= threshold {
				graph.AddEdge(i, j, weight)
			}
		}
	}
	return graph //pointer to a graph
}

// Adapted from (github.com/gilchrisn/graph-clustering-service/pkg2/louvain)
// RunLouvainOnMatrix: runs Louvain algorithm to mathematically determine communities (clusters of nodes) in a graph network to
// determine which genes are more strongly correlated to one another than to other genes.
// Function takes as input a pointer to louvain graph type
// outputs a map of nodeIDs to their community(cluster) ID and a measure of their modularity
func RunLouvainOnMatrix(graph *louvain.Graph) (map[int]int, float64) {

	// Configure Louvain
	config := louvain.NewConfig()
	config.Set("algorithm.max_levels", 5)
	config.Set("algorithm.min_modularity_gain", 1e-5)

	// Run clustering
	ctx := context.Background()
	result, err := louvain.Run(graph, config, ctx)
	if err != nil {
		log.Fatal(err)
	}

	// FinalCommunities: map[nodeID]communityID
	return result.FinalCommunities, result.Modularity
}

// Function CountCommunities counts the number of unique communities(clusters) produced from the Louvain algorithm
// takes as input: map of nodeIDs to their community(cluster) ID
// output: count of number of communities
func CountCommunities(communityMap map[int]int) int {
	count := 0

	//default to fault
	seen := make(map[int]bool)

	//range over the map and extract unique communities
	for _, community := range communityMap {
		if seen[community] != true {
			seen[community] = true
			count++ //count of unique communities
		}
	}

	return count
}

// Analyzing the graph networks

// Function CalculateNumEdges takes a graph as input and calculates the total number of edges in the graph
func CalculateNumEdges(BreastGraph GraphNetwork) int {
	count := 0
	for _, node := range BreastGraph {
		fromID := node.ID

		for _, edge := range node.Edges {
			toID := edge.To.ID

			// Only count each undirected edge once
			// count the "smaller -> larger" direction
			if fromID >= toID {
				continue
			}

			if edge.Weight != 0 {
				count++
			}
		}
	}

	return count
}

// Function ComputeAverageDegree takes as input total number of edges and total number of nodes in a graph,
// computes the average degree of a graph network with Average Degree = (2 * Total Edges) / Total Nodes
// and outputs the average degree of a graph
func ComputeAverageDegree(numEdges, numNodes int) float64 {
	return 2.0 * float64(numEdges) / float64(numNodes)
}

// Function ComputeEdgeDensity computes the edge density of a graph to quantify how connected a graph is.
// It takes the total number of edges and total number of nodes in the graph as input, computes the formula
// ((2*numEdges)/(numNodes*(numNodes-1))) and returns the measure of edge density.
func ComputeEdgeDensity(numEdges, numNodes int) float64 {
	return 2.0 * float64(numEdges) / (float64(numNodes) * float64(numNodes-1))
}

// Function EdgeStats computes the number of positively and negatively correlated edges in the graph input
func EdgeStats(graph GraphNetwork) (pos, neg int) {
	for _, node := range graph {
		fromID := node.ID

		for _, edge := range node.Edges {
			toID := edge.To.ID

			// Only count each undirected edge once
			// count the "smaller -> larger" direction
			if fromID >= toID {
				continue
			}

			if edge.Weight > 0 {
				pos++
			} else if edge.Weight < 0 {
				neg++
			}
		}
	}
	return pos, neg
}

// Function ComputeModuleSizes takes a map[nodeID]communityID and returns
// a map[communityID]moduleSize to analyze module size
func ComputeModuleSizes(clusterMap map[int]int) map[int]int {
	moduleSizes := make(map[int]int)

	//range over the cluster map, and fill in new map
	for _, communityID := range clusterMap {
		moduleSizes[communityID]++ //track count size for each module
	}
	return moduleSizes //return map
}

// Function LargestModule takes as input a map[communityID]moduleSize and returns the communityID with the
// largest size as well as its size
// If moduleSizes is empty, it returns (-1, 0).
func LargestModule(moduleSizes map[int]int) (int, int) {
	//guard against an empty map
	if len(moduleSizes) == 0 {
		return -1, 0
	}
	//initialize maximum community and maximum size
	maxCommunity := -1
	maxSize := -1

	//range over the map to find the largest module and its size
	for communityID, moduleSize := range moduleSizes {
		if moduleSize > maxSize {
			maxSize = moduleSize
			maxCommunity = communityID
		}
	}

	return maxCommunity, maxSize
}

// *ChatGPT generated*
// SummarizeModuleSizes computes basic statistics on module sizes.
// Inputs:
//   - moduleSizes: map[communityID]moduleSize
//   - totalGenes: total number of nodes in the graph (for % calculations)
//
// Outputs:
//   - minSize, maxSize, meanSize, medianSize: module size stats
//   - numModules: total number of modules
//   - numSmall: number of modules with size < 5
//   - largestModuleID: ID of the largest module
//   - largestPct: % of genes in the largest module
func SummarizeModuleSizes(
	moduleSizes map[int]int,
	totalGenes int,
) (
	minSize, maxSize, meanSize, medianSize float64,
	numModules, numSmall int,
	largestModuleID int,
	largestPct float64,
) {
	numModules = len(moduleSizes)
	if numModules == 0 {
		return
	}

	sizes := make([]int, 0, numModules)

	sum := 0
	first := true
	var minSizeInt, maxSizeInt, largestSize int

	for commID, size := range moduleSizes {
		sizes = append(sizes, size)
		sum += size

		if first {
			minSizeInt = size
			maxSizeInt = size
			largestModuleID = commID
			largestSize = size
			first = false
		} else {
			if size < minSizeInt {
				minSizeInt = size
			}
			if size > maxSizeInt {
				maxSizeInt = size
				largestModuleID = commID
				largestSize = size
			}
		}

		if size < 5 {
			numSmall++
		}
	}

	minSize = float64(minSizeInt)
	maxSize = float64(maxSizeInt)
	meanSize = float64(sum) / float64(numModules)

	// median
	sort.Ints(sizes)
	if numModules%2 == 1 {
		medianSize = float64(sizes[numModules/2])
	} else {
		mid := numModules / 2
		medianSize = float64(sizes[mid-1]+sizes[mid]) / 2.0
	}

	if totalGenes > 0 {
		largestPct = 100.0 * float64(largestSize) / float64(totalGenes)
	}

	return
}

// Function InvertMapWithGeneNames converts a map[nodeID]communityID into
// map[communityID][]geneName, using the provided geneNames slice
// where geneNames[nodeID] gives the gene name.
func InvertMapWithGeneNames(clusterMap map[int]int, geneNames []string) map[int][]string {
	//initialize inverted map to store community ID to gene names mapping
	moduleMap := make(map[int][]string)

	//range over the input cluster map
	for nodeID, communityID := range clusterMap {
		//get the gene name using the nodeID as index into geneNames slice
		gene := geneNames[nodeID]
		//build the list of gene names for each module key in the map
		moduleMap[communityID] = append(moduleMap[communityID], gene)
	}

	return moduleMap
}

// Function CountOverlap takes as input 2 lists of genes corresponding to 2 different modules from different cancer types
// and returns the number of overlapping genes (intersection count) between the two lists
func CountOverlap(a, b []string) int {
	// Build set for A
	setA := make(map[string]bool, len(a))
	//range over list of genes and put each gene into the map as a key and set its value to true
	for _, x := range a {
		setA[x] = true
	}

	// Count intersection |A ∩ B|
	intersection := 0
	for _, y := range b {
		if setA[y] { //if that gene already exists in the map then it is an overlap
			intersection++ //increase the count
		}
	}
	//return overlapping count
	return intersection
}

// Function Jaccard takes as input 2 lists of genes corresponding to 2 different modules from different cancer types
// and computes the Jaccard index by (Number of overlapping genes) / (Total number of unique genes in both modules)
// and returns the float64 index
func Jaccard(a, b []string) float64 {
	// Build set for A
	setA := make(map[string]bool, len(a))
	for _, x := range a {
		setA[x] = true
	}

	// Count intersection |A ∩ B|
	intersection := 0
	for _, y := range b {
		if setA[y] {
			intersection++
		}
	}

	// Compute union: |A| + |B| - |A ∩ B|
	union := len(a) + len(b) - intersection
	if union == 0 {
		return 0.0
	}
	// J = |A ∩ B| / |A ∪ B|
	return float64(intersection) / float64(union)
}

// Function JaccardMatrix takes 2 map[nodeID]communityID for each cancer type as input, and constructs a matrix for return corresponding to
// the Jaccard index between each pair of modules across each cancer type. It also returns lists of the breast and ovarian module IDs
func JaccardMatrix(breastClusterMap, ovarianClusterMap map[int]int, breastGeneNames, ovarianGeneNames []string) ([][]float64, []int, []int) {
	// Invert maps to map each module to a slice of gene names
	breastModuleMap := InvertMapWithGeneNames(breastClusterMap, breastGeneNames)
	ovarianModuleMap := InvertMapWithGeneNames(ovarianClusterMap, ovarianGeneNames)

	// Collect and sort module IDs to fix matrix ordering
	breastIDs := make([]int, 0, len(breastModuleMap))
	for id := range breastModuleMap {
		breastIDs = append(breastIDs, id)
	}
	sort.Ints(breastIDs) //sort

	ovarianIDs := make([]int, 0, len(ovarianModuleMap))
	for id := range ovarianModuleMap {
		ovarianIDs = append(ovarianIDs, id)
	}
	sort.Ints(ovarianIDs)

	// Initialize Jaccard matrix: rows = breast modules, cols = ovarian modules
	jaccardMatrix := make([][]float64, len(breastIDs))
	for i := range jaccardMatrix {
		jaccardMatrix[i] = make([]float64, len(ovarianIDs))
	}

	// Fill matrix with Jaccard(breast module, ovarian module)
	for i, breastID := range breastIDs {
		for j, ovarianID := range ovarianIDs {
			genesB := breastModuleMap[breastID]
			genesO := ovarianModuleMap[ovarianID]
			jaccardMatrix[i][j] = Jaccard(genesB, genesO)
		}
	}

	// Return matrix with jaccard indeces and the ID order for rows (breast mods) and columns (ovarian mods)
	return jaccardMatrix, breastIDs, ovarianIDs
}

// Function BestMatchesByRow takes a Jaccard matrix, row module IDs, and column module IDs as input
// and finds the best matching column module for each row module based on maximum Jaccard index.
// Returns a slice of ModuleMatch structs containing source module ID, target module ID, and Jaccard value.
func BestMatchesByRow(jaccardMatrix [][]float64, rowModuleIDs, colModuleIDs []int) []ModuleMatch {
	//initialize the list of ModuleMatch structs for return containing the best matches
	bestMatches := make([]ModuleMatch, 0, len(rowModuleIDs))

	//range over the row modules (breast cancer modules)
	for i, rowID := range rowModuleIDs {
		maxJ := -1.0    //initialize with -1.0 because minimum real jaccard index is 0
		bestColID := -1 //actual modules exist from 0

		//range over the column modules (ovarian cancer modules)
		for j, colID := range colModuleIDs {
			//find the maximum jaccard index pair corresponding to the best matched ovarian module for every breast module
			jacc := jaccardMatrix[i][j]
			if jacc > maxJ {
				maxJ = jacc
				bestColID = colID
			}
		}

		//build the struct for return
		bestMatches = append(bestMatches, ModuleMatch{
			SourceModuleID: rowID,
			TargetModuleID: bestColID,
			JaccardValue:   maxJ,
		})
	}

	return bestMatches
}

// Function BestMatchesByColumn takes a Jaccard matrix, row module IDs, and column module IDs as input
// and finds the best matching row module for each column module based on maximum Jaccard index.
// Returns a slice of ModuleMatch structs containing source module ID, target module ID, and Jaccard value.
func BestMatchesByColumn(jaccMat [][]float64, rowModuleIDs, colModuleIDs []int) []ModuleMatch {
	//initialize the list of ModuleMatch structs for return containing the best matches
	bestMatches := make([]ModuleMatch, 0, len(colModuleIDs))

	//range over the column modules (ovarian cancer modules)
	for j, colID := range colModuleIDs {
		maxJ := -1.0    //initialize with -1.0 because minimum real jaccard index is 0
		bestRowID := -1 //actual modules exist from 0

		//range over the row modules (breast cancer modules)
		for i, rowID := range rowModuleIDs {
			//find the maximum jaccard index pair corresponding to the best matched breast module for every ovarian module
			jacc := jaccMat[i][j]
			if jacc > maxJ {
				maxJ = jacc
				bestRowID = rowID
			}
		}

		//build the struct for return
		bestMatches = append(bestMatches, ModuleMatch{
			SourceModuleID: colID,
			TargetModuleID: bestRowID,
			JaccardValue:   maxJ,
		})
	}

	return bestMatches
}

// Function TransformJaccardMatrixToSlice takes a nonsymmetrical matrix as input and flattens
// it into a single slice of float64 values for return.
func TransformJaccardMatrixToSlice(JaccardMatrix [][]float64) []float64 {
	flattenedMatrix := make([]float64, 0)

	numRows := len(JaccardMatrix)
	numCols := len(JaccardMatrix[0])

	//range over the matrix
	for i := 0; i < numRows; i++ {
		for j := 0; j < numCols; j++ {
			value := JaccardMatrix[i][j]
			//only keep nonzero Jaccard values
			if value > 0 {
				flattenedMatrix = append(flattenedMatrix, value)
			}
		}
	}
	return flattenedMatrix
}

// Function ComputeJaccardQuantiles takes a matrix of Jaccard index values between cross-cancer modules
// and computes the 90th, 95th, and 99th percentile values.
// We define strongly overlapping module pairs as those with Jaccard index above the 95th percentile
// of all non-zero cross-cancer Jaccard values in the matrix.
func ComputeJaccardQuantiles(JaccardMatrix [][]float64) []float64 {
	flattenedMatrix := TransformJaccardMatrixToSlice(JaccardMatrix) //flattens matrix
	jaccardValueList := SortCorrVals(flattenedMatrix)               //sorts the jaccard values
	jaccardQuantiles := ComputeQuantile(jaccardValueList)           //computes 90, 95, 99th quantile probabilities from nonzero jaccard values
	return jaccardQuantiles
}

// Function ExtractPairsAboveJaccardThreshold takes as input the Jaccard matrix,
// the lists of breast and ovarian module IDs, and a Jaccard threshold.
// It extracts all cross-cancer module pairs whose Jaccard index is greater than
// or equal to the threshold and returns a slice of these module matches.
func ExtractPairsAboveJaccardThreshold(jaccardMatrix [][]float64, breastModuleIDs, ovarianModuleIDs []int, threshold float64) []ModuleMatch {
	topMatches := make([]ModuleMatch, 0)
	// range over all breast × ovarian module pairs
	for i := 0; i < len(jaccardMatrix); i++ {
		for j := 0; j < len(jaccardMatrix[0]); j++ {

			jValue := jaccardMatrix[i][j]

			// keep only pairs above threshold
			if jValue >= threshold {
				topMatches = append(topMatches, ModuleMatch{
					SourceModuleID: breastModuleIDs[i],
					TargetModuleID: ovarianModuleIDs[j],
					JaccardValue:   jValue,
				})
			}
		}
	}
	return topMatches
}

// Function TotalNumGenes takes as input lists of breast and ovarian gene names in the dataset and returns a count of
// the number of genes in the union of both sets.
func TotalNumGenes(breastGeneNames, ovarianGeneNames []string) int {
	//initialize map of gene names
	geneUniverse := make(map[string]bool)

	//range over each list and fill in the map (no repeated keys)
	for _, g := range breastGeneNames {
		geneUniverse[g] = true
	}
	for _, g := range ovarianGeneNames {
		geneUniverse[g] = true
	}

	return len(geneUniverse)
}

// *ChatGPT generated*
// BuildOverlapTables takes the significant module pairs and builds full
// contingency tables + metadata for each pair.
func BuildOverlapTables(
	modulePairs []ModuleMatch,
	breastModuleSizes, ovarianModuleSizes map[int]int,
	breastModuleMap, ovarianModuleMap map[int][]string,
	totalGenes int,
) []ModuleOverlapStats {

	overlapStats := make([]ModuleOverlapStats, 0, len(modulePairs))

	for _, p := range modulePairs {
		// module sizes
		sizeA := breastModuleSizes[p.SourceModuleID]
		sizeB := ovarianModuleSizes[p.TargetModuleID]

		// gene lists
		breastGenes := breastModuleMap[p.SourceModuleID]
		ovarianGenes := ovarianModuleMap[p.TargetModuleID]

		// compute overlap |A ∩ B|
		overlap := CountOverlap(breastGenes, ovarianGenes)

		// 2×2 table counts
		a := overlap
		b := sizeB - overlap
		c := sizeA - overlap
		d := totalGenes - (a + b + c)
		if d < 0 {
			d = 0 // guard, should not happen if totalGenes is correct
		}

		overlapStats = append(overlapStats, ModuleOverlapStats{
			BreastModuleID:  p.SourceModuleID,
			OvarianModuleID: p.TargetModuleID,
			SizeA:           sizeA,
			SizeB:           sizeB,
			Overlap:         overlap,
			A:               a,
			B:               b,
			C:               c,
			D:               d,
			TotalGenes:      totalGenes,
			Jaccard:         p.JaccardValue,
		})
	}

	return overlapStats
}

// *ChatGPT generated*
// ComputeOverlapPValues runs Fisher's exact test (right-tailed) for each 2x2 table
// and then applies Benjamini–Hochberg FDR correction.
//
// Right-tailed p-value here means:
//
//	P[X >= observed overlap | independence]
//
// which tests for enrichment of shared genes between modules.
func ComputeOverlapPValues(stats []ModuleOverlapStats) []ModuleOverlapResult {
	n := len(stats)
	results := make([]ModuleOverlapResult, n)
	pvals := make([]float64, n)

	// 1) raw Fisher p-values (right-tailed; enrichment of overlap)
	for i, s := range stats {
		// Table:
		//   n11 = A (in A and in B)
		//   n12 = B (not in A, in B)
		//   n21 = C (in A, not in B)
		//   n22 = D (not in A, not in B)
		_, _, rightp, _ := fisher.FisherExactTest(s.A, s.B, s.C, s.D)

		pvals[i] = rightp
		results[i] = ModuleOverlapResult{
			ModuleOverlapStats: s,
			PValue:             rightp,
		}
	}

	// 2) Benjamini–Hochberg FDR across all tested pairs
	qvals := benjaminiHochberg(pvals)
	for i := range results {
		results[i].QValue = qvals[i]
	}

	return results
}

// *ChatGPT generated*
// benjaminiHochberg applies BH FDR to a slice of p-values.
func benjaminiHochberg(p []float64) []float64 {

	n := len(p)
	arr := make([]idxP, n)
	for i, v := range p {
		arr[i] = idxP{idx: i, p: v}
	}

	// sort by p ascending
	sort.Slice(arr, func(i, j int) bool { return arr[i].p < arr[j].p })

	q := make([]float64, n)
	prev := 1.0

	// walk from largest p to smallest
	for k := n - 1; k >= 0; k-- {
		rank := float64(k + 1)
		val := arr[k].p * float64(n) / rank
		if val > 1.0 {
			val = 1.0
		}
		if val > prev {
			val = prev
		}
		prev = val
		q[arr[k].idx] = val
	}

	return q
}

// Function LocalClusteringCoeff takes a graph network as input, computes the clustering coefficient of each node (gene) with the formula:
// Ci = 2ni/ki(ki-1), where ni is the number of observed links connecting the ki neighbors of node i and ki(ki-1)/2
// is the total number of possible links. The function returns a slice of local clustering coefficients for the
// graph network.
func LocalClusteringCoeff(graph GraphNetwork) []float64 {
	clusteringCoeffs := make([]float64, len(graph))

	//range over the graph network (slice of pointer nodes)
	for i, node := range graph {

		// grab the count of neighbors from the number of edges each node has to other nodes
		k := len(node.Edges)
		//k cannot be 0 or 1 or the formula becomes undefined
		if k < 2 {
			clusteringCoeffs[i] = math.NaN() //set the coeff for these to NaN so we skip them
			continue
		}

		// build a set of neighbor IDs for quick lookup
		neighborSet := make(map[int]bool)
		for _, e := range node.Edges {
			neighborSet[e.To.ID] = true
		}

		// count edges among neighbors
		n := 0
		for _, edge := range node.Edges {
			u := edge.To.ID
			// check u's edges
			for _, edge2 := range graph[u].Edges {
				v := edge2.To.ID
				if u < v && neighborSet[v] {
					n++
				}
			}
		}

		// Then, compute the clustering coefficient from formula
		clusteringCoeffs[i] = (2.0 * float64(n)) / (float64(k) * float64(k-1))
	}

	return clusteringCoeffs
}

// Function GlobalClusteringCoeff takes as input a list of local clustering coefficients and computes
// the mean of these values for all nodes with at least two neighbors. This outputs a global measure
// of how “clustered” the whole network is.
func GlobalClusteringCoeff(clusteringCoeffs []float64) float64 {
	sum := 0.0
	nonNaNCount := 0.0

	//find number of nonNaN clustering coeffs
	for i := 0; i < len(clusteringCoeffs); i++ {
		if !math.IsNaN(clusteringCoeffs[i]) {
			nonNaNCount++
			//sum nonNaN clustering coeffs
			sum = sum + clusteringCoeffs[i]
		}
	}
	//compute the average
	return sum / float64(nonNaNCount)
}
