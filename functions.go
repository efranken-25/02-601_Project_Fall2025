package main

import (
	"02-601_Project_Fall2025/louvain"
	"context"
	"log"
	"math"
	"math/rand"
	"sort"
	"strconv"
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
				r := calculatePearsonCorrelation(gene1, gene2)
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

// Custom Pearson correlation calculation function
// Function calculatePearsonCorrelation computes the Pearson correlation coefficient between two slices of float64 values using the formula:
// r = Σ[(xi - x̄)(yi - ȳ)] / √[Σ(xi - x̄)²] × √[Σ(yi - ȳ)²]
// where x̄ and ȳ are the means of x and y respectively
// The function takes two slices of float64 values as input and returns the correlation coefficient as float64
func calculatePearsonCorrelation(x, y []float64) float64 {
	n := len(x)
	// Check for invalid inputs: mismatched lengths or empty slices
	if n != len(y) || n == 0 {
		return math.NaN()
	}

	// Calculate means: x̄ = Σxi/n and ȳ = Σyi/n
	sumX, sumY := 0.0, 0.0
	for i := 0; i < n; i++ {
		sumX += x[i]
		sumY += y[i]
	}
	meanX := sumX / float64(n)
	meanY := sumY / float64(n)

	// Calculate the numerator Σ[(xi - x̄)(yi - ȳ)] and denominators Σ(xi - x̄)² and Σ(yi - ȳ)²
	numerator := 0.0           // Σ[(xi - x̄)(yi - ȳ)] - covariance sum
	sumSqX, sumSqY := 0.0, 0.0 // Σ(xi - x̄)² and Σ(yi - ȳ)² - variance sums

	for i := 0; i < n; i++ {
		dx := x[i] - meanX   // deviation from mean for x
		dy := y[i] - meanY   // deviation from mean for y
		numerator += dx * dy // accumulate covariance
		sumSqX += dx * dx    // accumulate x variance
		sumSqY += dy * dy    // accumulate y variance
	}

	// Calculate correlation coefficient: r = numerator / √(sumSqX × sumSqY)
	denominator := math.Sqrt(sumSqX * sumSqY)
	if denominator == 0 {
		return math.NaN() // Handle division by zero (perfect constant values)
	}

	return numerator / denominator
}

// Function TransformMatrixToSlice takes a symmetrical square matrix as input and flattens
// it into a single slice of float64 values for return.
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

	// call the CalculateQuantile operation to compute the 90 95 99% quantile probabilities
	quantileProbability := []float64{0.90, 0.95, 0.99}

	computedQuantiles := make([]float64, len(quantileProbability))
	//compute the 90, 95, and 99 quantile probability
	for i, p := range quantileProbability {
		computedQuantiles[i] = CalculateQuantile(quantileSlice, p)
	}

	return computedQuantiles
}

// Function CalculateQuantile computes the quantile (output) for a given probability p (input) from a sorted slice (input)
func CalculateQuantile(sortedData []float64, p float64) float64 {
	n := len(sortedData)
	//check that the data slice is not empty or length 1, if it is then return 0.0 probability
	if n == 0 {
		return 0.0
	}
	if n == 1 {
		return sortedData[0]
	}

	// Calculate the position: (n-1) * p
	pos := float64(n-1) * p

	// Get the integer and fractional parts
	lower := int(pos)
	fraction := pos - float64(lower)

	// Handle edge cases
	if lower >= n-1 {
		return sortedData[n-1]
	}
	if lower < 0 {
		return sortedData[0]
	}

	// Linear interpolation between the two closest values
	return sortedData[lower] + fraction*(sortedData[lower+1]-sortedData[lower])
}

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

// Function CalculateNumEdges takes a graph network as input and calculates the total number of edges in the graph for return
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

// Function ComputeModuleEdgeCounts takes a graph network and a clusterMap[nodeID]communityID as input
// and returns a map[communityID]numEdges counting only edges whose
// both endpoints lie in the same community
func ComputeModuleEdgeCounts(graph GraphNetwork, clusterMap map[int]int) map[int]int {
	//initialize the map of communities to its number of edges
	moduleEdges := make(map[int]int)

	//range over the graph network
	for _, node := range graph {
		//get the source node ID and its community assignment
		fromID := node.ID
		commFrom := clusterMap[fromID]

		//range over all edges from this node
		for _, edge := range node.Edges {
			//get the target node ID and its community assignment
			toID := edge.To.ID
			commTo := clusterMap[toID]

			// Only count edges fully inside one community (intra-community edges)
			if commFrom != commTo {
				continue //skip inter-community edges
			}

			// For undirected graph: only count smaller to larger once to avoid double counting
			if fromID >= toID {
				continue
			}

			//only count edges with non-zero weight
			if edge.Weight != 0 {
				moduleEdges[commFrom]++ //increment edge count for this community
			}
		}
	}

	return moduleEdges
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

// Function ComputeModuleDensities takes a map of community to node count and map of community to edge count as inputs
// and returns a map of community to its edge density
func ComputeModuleDensities(moduleSizes map[int]int, moduleEdges map[int]int) map[int]float64 {

	//initialize the map for storing community densities for return
	moduleDensities := make(map[int]float64, len(moduleSizes))

	//range over the module sizes map to compute density for each community
	for communityID, numNodes := range moduleSizes {
		//get the number of edges for this community from the module edges map
		numEdges := moduleEdges[communityID]

		//only compute density if the module has more than one node and has edges
		if numNodes > 1 && numEdges > 0 {
			//compute edge density using the existing ComputeEdgeDensity function
			moduleDensities[communityID] = ComputeEdgeDensity(numEdges, numNodes)
		} else {
			// For singleton modules (1 node) or edgeless modules, density is 0
			moduleDensities[communityID] = 0.0
		}
	}

	return moduleDensities
}

// Function EdgeStats takes a graph network as input and computes the number of positively and negatively
// correlated edges in the graph (output)
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
// SummarizeModuleSizes computes comprehensive statistical measures on module (community) sizes using standard descriptive statistics.
// This function analyzes the distribution of cluster sizes to understand network modularity structure.
// Inputs:
//   - moduleSizes: map[communityID]moduleSize - mapping from each community ID to its node count
//   - totalGenes: total number of nodes in the graph (for percentage calculations)
//
// Outputs:
//   - minSize, maxSize: minimum and maximum module sizes in the network
//   - meanSize: arithmetic mean μ = (Σxi)/n where xi is size of module i, n is number of modules
//   - medianSize: middle value when module sizes are sorted; for even n: median = (x[n/2-1] + x[n/2])/2
//   - numModules: total count of distinct modules/communities
//   - numSmall: count of modules with size < 5 (small module threshold)
//   - largestModuleID: community ID of the module with maximum size
//   - largestPct: percentage of total genes contained in largest module = (largestSize/totalGenes) × 100
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
	// Handle edge case: empty module map returns all zero values
	if numModules == 0 {
		return
	}

	// Initialize slice to store all module sizes for median calculation
	sizes := make([]int, 0, numModules)

	// Initialize variables for statistical calculations
	sum := 0                                    // Σxi for mean calculation
	first := true                               // flag for first iteration initialization
	var minSizeInt, maxSizeInt, largestSize int // tracking min/max and largest module

	// Single pass through moduleSizes map to collect statistics
	for commID, size := range moduleSizes {
		sizes = append(sizes, size) // store for median calculation
		sum += size                 // accumulate sum for mean: Σxi

		// Initialize or update min/max tracking on first iteration
		if first {
			minSizeInt = size
			maxSizeInt = size
			largestModuleID = commID // track ID of largest module
			largestSize = size
			first = false
		} else {
			// Update minimum module size
			if size < minSizeInt {
				minSizeInt = size
			}
			// Update maximum module size and largest module tracking
			if size > maxSizeInt {
				maxSizeInt = size
				largestModuleID = commID // update ID of largest module
				largestSize = size
			}
		}

		// Count small modules (size < 5 threshold)
		if size < 5 {
			numSmall++
		}
	}

	// Convert integer statistics to float64 for output
	minSize = float64(minSizeInt)
	maxSize = float64(maxSizeInt)
	// Calculate arithmetic mean: μ = Σxi/n
	meanSize = float64(sum) / float64(numModules)

	// Calculate median: middle value of sorted distribution
	sort.Ints(sizes) // sort module sizes in ascending order
	if numModules%2 == 1 {
		// Odd number of modules: median = middle element
		medianSize = float64(sizes[numModules/2])
	} else {
		// Even number of modules: median = average of two middle elements
		mid := numModules / 2
		medianSize = float64(sizes[mid-1]+sizes[mid]) / 2.0
	}

	// Calculate percentage of genes in largest module
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

// Function MapIntValuesToFloatSlice converts a map[int]int (input)
// into a []float64 of just the values (output)
func MapIntValuesToFloatSlice(m map[int]int) []float64 {
	out := make([]float64, 0, len(m))
	//range over the map keys and extract the values to append to growing slice
	for _, v := range m {
		out = append(out, float64(v))
	}
	return out
}

// Function MapFloatValuesToSlice converts a map[int]float64 (input)
// into a []float64 of just the values (output)
func MapFloatValuesToSlice(m map[int]float64) []float64 {
	out := make([]float64, 0, len(m))
	//range over the map keys and extract the values to append to growing slice
	for _, v := range m {
		out = append(out, v)
	}
	return out
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

// Function FilterNaNs takes a slice of floats as input and returns a new slice
// with all NaN values removed.
func FilterNaNs(vals []float64) []float64 {
	out := make([]float64, 0, len(vals))
	//range through the input
	for _, v := range vals {
		//if it is not NaN then go ahead and input it into the out slice
		if !math.IsNaN(v) {
			out = append(out, v)
		}
	}
	return out
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

// Function ComputeDegrees takes a graph network structure as input and computes for each node in the graph its
// degree = its number of connected edges. It returns a slice containing a list of these values (# of degrees) for each node in the graph.
func ComputeDegrees(graph GraphNetwork) []float64 {
	degree := make([]float64, len(graph))
	//range over the graph network
	for i := range graph {
		//grab the count of its neighbors
		//and input the node's degree count in the slice storing all degrees for the graph network
		degree[i] = float64(len(graph[i].Edges))
	}
	return degree
}

// Function MeanStd computes the mean and standard deviation (output) of a slice of float64 values (input)
func MeanStd(vals []float64) (mean float64, sd float64) {
	//find the number of total values
	n := float64(len(vals))
	sum := 0.0
	//sum the values
	for _, v := range vals {
		sum = sum + v
	}
	//divide by the number of values to compute the mean
	mean = sum / n

	// compute variance = Sum of each (datapoint - mean)^2
	varSum := 0.0

	for _, v := range vals {
		varSum += (v - mean) * (v - mean)
	}
	//population standard deviation = variance sum / number of datapoints
	sd = math.Sqrt(varSum / n)
	//return mean and standard deviation
	return mean, sd
}

// Function RandomGraphGenerator creates a random graph using the Erdős-Rényi G(n,p) model.
// For each possible pair of nodes (genes), an edge is included independently with probability pVal.
// The algorithm examines all (n choose 2) = n(n-1)/2 possible edges and includes each with fixed probability p.
// This generates an undirected, unweighted random graph following the G(n,p) distribution.
// Inputs: geneNames (slice of gene names for node labels), pVal (edge creation probability)
// Outputs: random GraphNetwork with randomly connected genes where each edge has weight 1.0
func RandomGraphGenerator(geneNames []string, pVal float64) GraphNetwork {
	n := len(geneNames)

	// Initialize nodes for the random graph
	graph := make(GraphNetwork, n)
	// Create new pointer nodes for every gene in this network
	for i := 0; i < n; i++ {
		graph[i] = &Node{
			ID:       i,
			GeneName: geneNames[i],
			Edges:    []*Edge{},
		}
	}

	// Iterate over all possible unordered pairs of nodes (i,j)
	// For each pair, include an edge with probability pVal
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if rand.Float64() < pVal { // edge created with probability pVal
				u := graph[i]
				v := graph[j]

				// Create undirected edge: add edge in both directions with weight 1.0
				e1 := &Edge{To: v, Weight: 1.0}
				e2 := &Edge{To: u, Weight: 1.0}
				//add the edge to that node's edge list
				u.Edges = append(u.Edges, e1)
				v.Edges = append(v.Edges, e2)
			}
		}
	}
	//return the random graph
	return graph
}

// Function calculateECDF calculates the Empirical Cumulative Distribution Function (ECDF) value
// for a given sample (distribution) at a specific point x using the formula: ECDF(x) = (number of values ≤ x) / n
func calculateECDF(sample []float64, x float64) float64 {
	count := 0
	//range over the sample and count the number of values that are <= x
	for _, value := range sample {
		if value <= x {
			count++
		}
	}
	//divide the count by total sample size for the Empirical CDF
	return float64(count) / float64(len(sample))
}

// Function KSTwoSampleStatistic calculates the Kolmogorov-Smirnov D statistic for two samples
// using the formula: D = max|F₁(x) - F₂(x)| where F₁ and F₂ are empirical CDFs
func KSTwoSampleStatistic(sample1, sample2 []float64) float64 {

	length1 := len(sample1)
	length2 := len(sample2)
	combinedLength := length1 + length2

	// Combine both lists
	combined := make([]float64, 0, combinedLength)
	combined = append(combined, sample1...)
	combined = append(combined, sample2...)

	//Sort the combined list
	sort.Float64s(combined)

	// Remove duplicates from combined slice to only check unique points
	uniqueCombined := make([]float64, 0, len(combined))
	//add the first value to the unique list
	uniqueCombined = append(uniqueCombined, combined[0])
	//range through the combined list and if the value is not equal to the previous one,
	//then it is not a duplicate and you can add it to unique list
	for i := 1; i < len(combined); i++ {
		if combined[i] != combined[i-1] {
			uniqueCombined = append(uniqueCombined, combined[i])
		}
	}

	// Find maximum absolute difference between empirical CDFs
	maxD := 0.0
	//calculate the empirical CDF for every value in the unique list
	for _, x := range uniqueCombined {
		ecdf1 := calculateECDF(sample1, x)
		ecdf2 := calculateECDF(sample2, x)
		//calculate the absolute value of the difference between the CDFs
		d := math.Abs(ecdf1 - ecdf2)
		//find the maximum difference
		if d > maxD {
			maxD = d
		}
	}

	//maximum difference between empirical CDFs = Dstatistic
	return maxD
}

// Function KSTest performs the two-sample Kolmogorov-Smirnov test to determine if two datasets
// come from the same underlying distribution using the KS test statistic:
// D = max|F₁(x) - F₂(x)| where F₁ and F₂ are the empirical cumulative distribution functions
// The p-value is computed using the infinite series formula: p = 2*∑[k=1 to ∞] (-1)^(k-1)*e^(-2k² * λ²)
// where λ = √(neff) * D and neff = (n₁ × n₂)/(n₁ + n₂)
// Inputs: dataset1, dataset2 - two slices of float64 values to compare
// Outputs: dStatistic (KS test statistic), pValue (probability of observing this difference under null hypothesis)
func KSTest(dataset1, dataset2 []float64) (float64, float64) {

	// Calculate KS test statistic
	dStatistic := KSTwoSampleStatistic(dataset1, dataset2)

	// Calculate the effective sample size
	m := float64(len(dataset1)*len(dataset2)) / float64(len(dataset1)+len(dataset2))

	//Calculate lamda: λ=sqrt(neff​)*dStatistic
	lambda := math.Sqrt(m) * dStatistic

	//compute pvalue from formula: p=2*∑​(−1)^(k−1)*e^(−2k^2 * λ^2)
	sum := 0.0
	for k := 1; k <= 100; k++ {
		term := math.Pow(-1, float64(k-1)) * math.Exp(-2.0*float64(k*k)*lambda*lambda)
		sum += term
	}

	pValue := 2.0 * sum

	// Ensure p-value is positive
	if pValue < 0 {
		pValue = 0.0
	}

	//return the D statistic and p value
	return dStatistic, pValue
}
