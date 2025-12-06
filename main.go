package main

import (
	"fmt"
	"log"
	"sort"
)

func main() {

	//Begin Analyzing Breast Cancer RNA seq gene count Data
	// Create a map of maps. Nested map mapping gene_name (string) to sampleID(string) to tpm value
	geneExpressionBreastMap, err1 := ReadGeneExpressionDirToGeneMap("GeneExpressionData/BreastData")

	// log any errors that may arise from the parsing/ reading in the breast data files
	if err1 != nil {
		log.Fatal(err1)
	}

	fmt.Println("Breast Cancer Gene Expression Data read and map created.")

	//Begin Analyzing Ovarian Cancer RNA seq gene count Data
	//Create a map of maps. Nested map mapping gene_name (string) to sampleID(string) to tpm value
	geneExpressionOvarianMap, err2 := ReadGeneExpressionDirToGeneMap("GeneExpressionData/OvarianData")

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

	//Print the first few breast sample names to ensure they are sorted by numerical order
	fmt.Println(breastSampleNames[:5])

	fmt.Println("Gathering sample names for ovarian cancer data.")
	// Store the sample names for the Ovarian Cancer Data
	ovarianSampleNames := GetSampleNames(geneExpressionOvarianMap)

	//Sort the ovarian sample names
	sortSampleNames(ovarianSampleNames)

	//Print the ovarian sample names to ensure they are sorted by numerical order
	fmt.Println(ovarianSampleNames[:5])

	fmt.Println("Cleaning Data...Filtering out low expression genes...")
	// Use Mean-Based Filtering in order to normalize data while maintaining fixed number of samples across each gene
	filteredGeneExpressionBreastMap := MeanBasedFilter(geneExpressionBreastMap, breastSampleNames, 10.0)

	filteredGeneExpressionOvarianMap := MeanBasedFilter(geneExpressionOvarianMap, ovarianSampleNames, 10.0)

	fmt.Println("Gathering gene names for breast cancer data.")
	// Store the gene names for the Breast Cancer Data (sorted)
	breastGeneNames := GetGeneNames(filteredGeneExpressionBreastMap)

	//Print the first few breast gene names to check the names and ensure they are sorted
	fmt.Println(breastGeneNames[:5])

	fmt.Println("Gathering gene names for ovarian cancer data.")
	// Store the gene names for the Breast Cancer Data
	ovarianGeneNames := GetGeneNames(filteredGeneExpressionOvarianMap)

	//Print the first few breast gene names to check the names and ensure they are sorted
	fmt.Println(ovarianGeneNames[:5])

	//Compute the Pearson Correlation on breast cancer data
	fmt.Println("Computing Pearson Correlation between genes across all breast cancer samples.")

	BreastPearsonCorrelationMatrix := ComputePearsonCorrelation(breastGeneNames, breastSampleNames, filteredGeneExpressionBreastMap)

	//Compute the Pearson Correlation on ovarian cancer data
	fmt.Println("Computing Pearson Correlation between genes across all ovarian cancer samples.")

	OvarianPearsonCorrelationMatrix := ComputePearsonCorrelation(ovarianGeneNames, ovarianSampleNames, filteredGeneExpressionOvarianMap)

	fmt.Println("Transforming the breast correlation matrix...")
	breastQuantileSlice := TransformMatrixToSlice(BreastPearsonCorrelationMatrix)

	fmt.Println("Transforming the ovarian correlation matrix...")
	ovarianQuantileSlice := TransformMatrixToSlice(OvarianPearsonCorrelationMatrix)

	sortedBreastQuantileSlice := SortCorrVals(breastQuantileSlice)

	sortedOvarianQuantileSlice := SortCorrVals(ovarianQuantileSlice)

	fmt.Println("Computing the breast quantile probabilities...")
	computedBreastQuantiles := ComputeQuantile(sortedBreastQuantileSlice)

	fmt.Println("Computing the ovarian quantile probabilities...")
	computedOvarianQuantiles := ComputeQuantile(sortedOvarianQuantileSlice)

	fmt.Println("Breast Cancer Quantile Probabilities:", computedBreastQuantiles)
	fmt.Println("Ovarian Cancer Quantile Probabilities:", computedOvarianQuantiles)

	//The 95% quantile probability chosen as the cutoff correlation value threshold
	breastThreshold := computedBreastQuantiles[1]
	ovarianThreshold := computedOvarianQuantiles[1]

	fmt.Println("Building the Graph Networks...")

	//Build weighted, undirected cancer graph network with edge weights corresponding to raw computed correlations between genes (nodes)
	BreastGraph := BuildGraph(BreastPearsonCorrelationMatrix, breastGeneNames, breastThreshold)
	OvarianGraph := BuildGraph(OvarianPearsonCorrelationMatrix, ovarianGeneNames, ovarianThreshold)

	////////////////////////////////////////////////////////////////

	//Analyzing Graph Properties (Network Level per cancer type)
	//Breast Cancer
	numBreastNodes := len(BreastGraph)
	numBreastEdges := CalculateNumEdges(BreastGraph)
	breastGraphDegree := ComputeAverageDegree(numBreastEdges, numBreastNodes)
	breastEdgeDensity := ComputeEdgeDensity(numBreastEdges, numBreastNodes)
	posBreastEdges, negBreastEdges := EdgeStats(BreastGraph)
	fmt.Printf("\n===== Breast Cancer Graph Properties =====\n")
	fmt.Printf("Nodes (N):                  %d\n", numBreastNodes)
	fmt.Printf("Edges (E):                  %d\n", numBreastEdges)
	fmt.Printf("Average Degree:             %.3f\n", breastGraphDegree)
	fmt.Printf("Edge Density:               %.6f\n", breastEdgeDensity)
	fmt.Printf("Positive Edges:             %d (%.3f%%)\n",
		posBreastEdges, 100*float64(posBreastEdges)/float64(numBreastEdges))
	fmt.Printf("Negative Edges:             %d (%.3f%%)\n",
		negBreastEdges, 100*float64(negBreastEdges)/float64(numBreastEdges))

	//Ovarian Cancer
	numOvarianNodes := len(OvarianGraph)
	numOvarianEdges := CalculateNumEdges(OvarianGraph)
	ovarianGraphDegree := ComputeAverageDegree(numOvarianEdges, numOvarianNodes)
	ovarianEdgeDensity := ComputeEdgeDensity(numOvarianEdges, numOvarianNodes)
	posOvarianEdges, negOvarianEdges := EdgeStats(OvarianGraph)
	fmt.Printf("\n===== Ovarian Cancer Graph Properties =====\n")
	fmt.Printf("Nodes (N):                  %d\n", numOvarianNodes)
	fmt.Printf("Edges (E):                  %d\n", numOvarianEdges)
	fmt.Printf("Average Degree:             %.3f\n", ovarianGraphDegree)
	fmt.Printf("Edge Density:               %.6f\n", ovarianEdgeDensity)
	fmt.Printf("Positive Edges:             %d (%.3f%%)\n",
		posOvarianEdges, 100*float64(posOvarianEdges)/float64(numOvarianEdges))
	fmt.Printf("Negative Edges:             %d (%.3f%%)\n",
		negOvarianEdges, 100*float64(negOvarianEdges)/float64(numOvarianEdges))

	//Build Louvain Graph type of cancer data to run the Louvain Algorithm
	breastLouvainGraph := BuildLouvainGraph(BreastPearsonCorrelationMatrix, breastThreshold)
	ovarianLouvainGraph := BuildLouvainGraph(OvarianPearsonCorrelationMatrix, ovarianThreshold)

	fmt.Println("Running Louvain algorithm and assigning gene community clusters...")
	//Run Louvain Algorithm from package to determine communities assignments (clusters of nodes) for further analysis and visualization
	breastClusterMap, breastModularity := RunLouvainOnMatrix(breastLouvainGraph)
	ovarianClusterMap, ovarianModularity := RunLouvainOnMatrix(ovarianLouvainGraph)

	////////////////////////////////////////////////////////////////
	//Module level Analyses per cancer type

	// Breast Cancer
	breastModuleSizes := ComputeModuleSizes(breastClusterMap)
	minB, maxB, meanB, medianB,
		numModsB, numSmallB,
		largestModB, largestPctB :=
		SummarizeModuleSizes(breastModuleSizes, numBreastNodes)

	fmt.Printf("\n===== Breast Cancer Module Properties =====\n")
	fmt.Printf("Number of modules:                       %d\n", numModsB)
	fmt.Printf("Module size (min / median / mean / max): %.0f / %.0f / %.1f / %.0f genes\n",
		minB, medianB, meanB, maxB)
	fmt.Printf("Modules with < 5 genes:                  %d\n", numSmallB)
	fmt.Printf("Largest module:                          ID %d with %.0f genes (%.2f%% of genes)\n",
		largestModB, maxB, largestPctB)

	// Ovarian Cancer
	ovarianModuleSizes := ComputeModuleSizes(ovarianClusterMap)
	minO, maxO, meanO, medianO,
		numModsO, numSmallO,
		largestModO, largestPctO :=
		SummarizeModuleSizes(ovarianModuleSizes, numOvarianNodes)

	fmt.Printf("\n===== Ovarian Cancer Module Properties =====\n")
	fmt.Printf("Number of modules:                       %d\n", numModsO)
	fmt.Printf("Module size (min / median / mean / max): %.0f / %.0f / %.1f / %.0f genes\n",
		minO, medianO, meanO, maxO)
	fmt.Printf("Modules with < 5 genes:                  %d\n", numSmallO)
	fmt.Printf("Largest module:                          ID %d with %.0f genes (%.2f%% of genes)\n",
		largestModO, maxO, largestPctO)

	// Access results for modularity and communities(clusters)
	fmt.Printf("Final breast cancer graph modularity: %.4f\n", breastModularity)
	fmt.Printf("Number of breast cancer communities: %d\n", CountCommunities(breastClusterMap))
	fmt.Printf("Final ovarian cancer modularity: %.4f\n", ovarianModularity)
	fmt.Printf("Number of ovarian cancer communities: %d\n", CountCommunities(ovarianClusterMap))

	////////////////////////////////////////////////////////////////
	//Analyses between Cancer type

	fmt.Println("Conducting analyses between cancer types...")
	fmt.Println("\n===== Jaccard Similarity: Breast vs Ovarian Modules =====")

	// Compute Jaccard matrix between all breast and ovarian modules
	jaccardMatrix, breastModuleIDs, ovarianModuleIDs := JaccardMatrix(breastClusterMap, ovarianClusterMap, breastGeneNames, ovarianGeneNames)

	// Find which ovarian module is most similar to each breast module
	breastToOvarianMatches := BestMatchesByRow(jaccardMatrix, breastModuleIDs, ovarianModuleIDs)

	fmt.Println("\n--- Best Ovarian Match for Each Breast Module ---")
	for _, match := range breastToOvarianMatches {
		fmt.Printf("Breast Module %d -> Ovarian Module %d : Jaccard = %.4f\n",
			match.SourceModuleID, match.TargetModuleID, match.JaccardValue)
	}

	// Find which breast module is most similar to each ovarian module
	ovarianToBreastMatches := BestMatchesByColumn(jaccardMatrix, breastModuleIDs, ovarianModuleIDs)

	fmt.Println("\n--- Best Breast Match for Each Ovarian Module ---")
	for _, match := range ovarianToBreastMatches {
		fmt.Printf("Ovarian Module %d -> Breast Module %d : Jaccard = %.4f\n",
			match.SourceModuleID, match.TargetModuleID, match.JaccardValue)
	}

	//compute quantiles from the jaccard matrix to determine a threshold for significant module similarity
	jaccardQuantiles := ComputeJaccardQuantiles(jaccardMatrix)
	jaccardThreshold := jaccardQuantiles[1] // 95th percentile used as threshold for significance

	// Extract all module pairs with Jaccard index >= jaccardThreshold to identify highly similar cross-cancer modules
	modulePairs := ExtractPairsAboveJaccardThreshold(jaccardMatrix, breastModuleIDs, ovarianModuleIDs, jaccardThreshold)
	fmt.Printf("\n===== Jaccard Threshold (95th percentile): %.4f =====\n", jaccardThreshold)

	fmt.Println("\n===== Cross-Cancer Module Pairs Above Threshold =====")

	// total gene universe size
	totalGenes := TotalNumGenes(breastGeneNames, ovarianGeneNames)

	// invert cluster maps to get lists of gene names for each module ID
	breastModuleMap := InvertMapWithGeneNames(breastClusterMap, breastGeneNames)
	ovarianModuleMap := InvertMapWithGeneNames(ovarianClusterMap, ovarianGeneNames)

	// build full overlap + contingency info for each significant pair
	overlapStats := BuildOverlapTables(
		modulePairs,
		breastModuleSizes, ovarianModuleSizes,
		breastModuleMap, ovarianModuleMap,
		totalGenes,
	)

	// compute Fisher p-values (right-tailed; enrichment of overlap) + BH FDR
	overlapResults := ComputeOverlapPValues(overlapStats)

	// sort by q-value (or p-value) for nicer printing
	sort.Slice(overlapResults, func(i, j int) bool {
		return overlapResults[i].QValue < overlapResults[j].QValue
	})

	// print descriptive info about each significant cross-cancer module pair
	for _, r := range overlapResults {
		fmt.Printf(
			"Breast Module %d (|A|=%d) ↔ Ovarian Module %d (|B|=%d): Overlap=%d, Jaccard=%.4f, p=%.3e, q=%.3e\n",
			r.BreastModuleID, r.SizeA,
			r.OvarianModuleID, r.SizeB,
			r.Overlap, r.Jaccard,
			r.PValue, r.QValue,
		)
	}

	////////////////////////////////////////////////////////////////
	//Computing Clustering Coefficients for each Graph Network
	fmt.Println("Computing Clustering Coefficients for each Graph Network...")
	breastLocalC := LocalClusteringCoeff(BreastGraph)
	ovarianLocalC := LocalClusteringCoeff(OvarianGraph)

	breastGlobalC := GlobalClusteringCoeff(breastLocalC)
	ovarianGlobalC := GlobalClusteringCoeff(ovarianLocalC)
	fmt.Println("\n===== Global Clustering Coefficients =====")
	fmt.Printf("Breast Cancer global clustering coefficient:  %.4f\n", breastGlobalC)
	fmt.Printf("Ovarian Cancer global clustering coefficient: %.4f\n", ovarianGlobalC)

	////////////////////////////////////////////////////////////////

	//Identifying Hub Genes

	///////////////////////////////////////////////////////////////
	//Write gene network and community assignments to CSV files
	fmt.Println("Writing breast cancer gene network and community assignments to CSV files...")

	err3 := WriteCommunitiesCSV("ShinyApp/breast_nodes_communities.csv", breastGeneNames, breastClusterMap)

	if err3 != nil {
		panic(err3)
	}

	err4 := WriteEdgesCSV("ShinyApp/breast_edges.csv", BreastGraph)

	if err4 != nil {
		panic(err4)
	}

	fmt.Println("Writing ovarian cancer gene network and community assignments to CSV files...")

	err5 := WriteCommunitiesCSV("ShinyApp/ovarian_nodes_communities.csv", ovarianGeneNames, ovarianClusterMap)

	if err5 != nil {
		panic(err5)
	}

	err6 := WriteEdgesCSV("ShinyApp/ovarian_edges.csv", OvarianGraph)

	if err6 != nil {
		panic(err6)
	}
}
