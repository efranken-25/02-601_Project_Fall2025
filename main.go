package main

import (
	"fmt"
	"log"
	"math/rand"
	"os"
	"sort"
	"strconv"
)

const (
	//USER: Label your datasets here.
	dataset1Label = "Breast Cancer"
	dataset2Label = "Ovarian Cancer"

	dataset1Path = "GeneExpressionData/Dataset1"
	dataset2Path = "GeneExpressionData/Dataset2"

	outputDir = "ShinyApp"
)

func main() {
	// Parse command-line arguments - simple approach
	arg := 1 // default to both datasets
	if len(os.Args) > 1 {
		if argValue, err := strconv.Atoi(os.Args[1]); err == nil {
			arg = argValue
		}
	}

	//Begin Analyzing Dataset1 RNA seq gene count Data
	// Create a map of maps. Nested map mapping gene_name (string) to sampleID(string) to tpm value

	geneExpressionBreastMap, err1 := ReadGeneExpressionDirToGeneMap(dataset1Path)

	// log any errors that may arise from the parsing/ reading in the breast data files
	if err1 != nil {
		log.Fatal(err1)
	}

	fmt.Printf("%s Gene Expression Data read and map created.\n", dataset1Label)

	/*fmt.Println("Printing first line of geneExpressionBreastMap")
	  //Check the format of the data file by printing the first line of our geneExpressionBreastMap
	  for key, val := range geneExpressionBreastMap {
	      fmt.Println(key, val)
	      break
	  }*/

	fmt.Printf("Gathering sample names for %s data.\n", dataset1Label)

	// Store the sample names for the Breast Cancer Data
	breastSampleNames := GetSampleNames(geneExpressionBreastMap)

	//Sort the breast sample names
	sortSampleNames(breastSampleNames)

	//Print the first few breast sample names to ensure they are sorted by numerical order
	//fmt.Println(breastSampleNames[:5])

	fmt.Println("Cleaning Data...Filtering out low expression genes...")
	// Use Mean-Based Filtering in order to normalize data while maintaining fixed number of samples across each gene
	filteredGeneExpressionBreastMap := MeanBasedFilter(geneExpressionBreastMap, breastSampleNames, 10.0)

	fmt.Printf("Gathering gene names for %s data.\n", dataset1Label)

	// Store the gene names for the Breast Cancer Data (sorted)
	breastGeneNames := GetGeneNames(filteredGeneExpressionBreastMap)

	//Print the first few breast gene names to check the names and ensure they are sorted
	//fmt.Println(breastGeneNames[:5])

	//Compute the Pearson Correlation on dataset1 data
	fmt.Printf("Computing Pearson Correlation between genes across all %s samples.\n", dataset1Label)

	//Compute the Pearson Correlation on dataset1 data
	//fmt.Printf("Computing Pearson Correlation between genes across all %s samples.\n", dataset3Label)
	BreastPearsonCorrelationMatrix := ComputePearsonCorrelation(breastGeneNames, breastSampleNames, filteredGeneExpressionBreastMap)

	fmt.Printf("Transforming the %s correlation matrix...\n", dataset1Label)
	breastQuantileSlice := TransformMatrixToSlice(BreastPearsonCorrelationMatrix)

	//sort the correlation values
	sortedBreastQuantileSlice := SortCorrVals(breastQuantileSlice)

	//The 95% quantile probability chosen as the cutoff correlation value threshold

	fmt.Printf("Computing the %s quantile probabilities...\n", dataset1Label)
	computedBreastQuantiles := ComputeQuantile(sortedBreastQuantileSlice)
	fmt.Printf("%s Quantile Probabilities: %v\n", dataset1Label, computedBreastQuantiles)
	breastThreshold := computedBreastQuantiles[1]

	fmt.Println("Building the Graph Network...")

	//Build weighted, undirected cancer graph network with edge weights corresponding to raw computed correlations between genes (nodes)
	BreastGraph := BuildGraph(BreastPearsonCorrelationMatrix, breastGeneNames, breastThreshold)

	/*
	   //Begin Analyzing Dataset3 RNA seq gene count Data
	   // Create a map of maps. Nested map mapping gene_name (string) to sampleID(string) to tpm value
	   geneExpressionLungMap, err1 := ReadGeneExpressionDirToGeneMap(dataset3Path)

	   // log any errors that may arise from the parsing/ reading in the breast data files
	   if err1 != nil {
	       log.Fatal(err1)
	   }

	   fmt.Printf("%s Gene Expression Data read and map created.\n", dataset3Label)

	   fmt.Printf("Gathering sample names for %s data.\n", dataset3Label)

	   // Store the sample names for the Lung Cancer Data
	   lungSampleNames := GetSampleNames(geneExpressionLungMap)

	   //Sort the lung sample names
	   sortSampleNames(lungSampleNames)

	   fmt.Println("Cleaning Data...Filtering out low expression genes...")
	   // Use Mean-Based Filtering in order to normalize data while maintaining fixed number of samples across each gene
	   filteredGeneExpressionLungMap := MeanBasedFilter(geneExpressionLungMap, lungSampleNames, 20.0)

	   fmt.Printf("Gathering gene names for %s data.\n", dataset3Label)

	   // Store the gene names for the Lung Cancer Data (sorted)
	   lungGeneNames := GetGeneNames(filteredGeneExpressionLungMap)

	   //Compute the Pearson Correlation on dataset1 data
	   fmt.Printf("Computing Pearson Correlation between genes across all %s samples.\n", dataset3Label)
	   LungPearsonCorrelationMatrix := ComputePearsonCorrelation(lungGeneNames, lungSampleNames, filteredGeneExpressionLungMap)

	   fmt.Printf("Transforming the %s correlation matrix...\n", dataset3Label)
	   lungQuantileSlice := TransformMatrixToSlice(LungPearsonCorrelationMatrix)

	   //sort the correlation values
	   sortedLungQuantileSlice := SortCorrVals(lungQuantileSlice)

	   //The 95% quantile probability chosen as the cutoff correlation value threshold (high value)
	   fmt.Printf("Computing the %s quantile probabilities...\n", dataset3Label)
	   computedLungQuantiles := ComputeQuantile(sortedLungQuantileSlice)
	   fmt.Printf("%s Quantile Probabilities: %v\n", dataset3Label, computedLungQuantiles)
	   LungThreshold := computedLungQuantiles[2]

	   fmt.Println("Building the Graph Network...")

	   //Build weighted, undirected cancer graph network with edge weights corresponding to raw computed correlations between genes (nodes)
	   LungGraph := BuildGraph(LungPearsonCorrelationMatrix, lungGeneNames, LungThreshold)

	*/
	////////////////////////////////////////////////////////////////
	//Analyzing Graph Properties (Network Level per cancer type)

	//Dataset1 (Breast Cancer)
	numBreastNodes := len(BreastGraph)
	numBreastEdges := CalculateNumEdges(BreastGraph)
	breastGraphDegree := ComputeAverageDegree(numBreastEdges, numBreastNodes)
	breastEdgeDensity := ComputeEdgeDensity(numBreastEdges, numBreastNodes)
	posBreastEdges, negBreastEdges := EdgeStats(BreastGraph)
	fmt.Printf("\n===== %s Graph Properties =====\n", dataset1Label)
	fmt.Printf("Nodes (N):                  %d\n", numBreastNodes)
	fmt.Printf("Edges (E):                  %d\n", numBreastEdges)
	fmt.Printf("Average Degree:             %.3f\n", breastGraphDegree)
	fmt.Printf("Edge Density:               %.6f\n", breastEdgeDensity)
	fmt.Printf("Positive Edges:             %d (%.3f%%)\n",
		posBreastEdges, 100*float64(posBreastEdges)/float64(numBreastEdges))
	fmt.Printf("Negative Edges:             %d (%.3f%%)\n",
		negBreastEdges, 100*float64(negBreastEdges)/float64(numBreastEdges))
	// Compute node degrees for each graph (for node-level distributions)
	breastDegrees := ComputeDegrees(BreastGraph)

	//Build Louvain Graph type of cancer data to run the Louvain Algorithm
	breastLouvainGraph := BuildLouvainGraph(BreastPearsonCorrelationMatrix, breastThreshold)
	fmt.Println("Running Louvain algorithm and assigning gene community clusters...")
	//Run Louvain Algorithm from package to determine communities assignments (clusters of nodes) for further analysis and visualization
	breastClusterMap, breastModularity := RunLouvainOnMatrix(breastLouvainGraph)

	/*
	   //Dataset3 (Lung Cancer)
	   numLungNodes := len(LungGraph)
	   numLungEdges := CalculateNumEdges(LungGraph)
	   lungGraphDegree := ComputeAverageDegree(numLungEdges, numLungNodes)
	   lungEdgeDensity := ComputeEdgeDensity(numLungEdges, numLungNodes)
	   posLungEdges, negLungEdges := EdgeStats(LungGraph)
	   fmt.Printf("\n===== %s Graph Properties =====\n", dataset3Label)
	   fmt.Printf("Nodes (N):                  %d\n", numLungNodes)
	   fmt.Printf("Edges (E):                  %d\n", numLungEdges)
	   fmt.Printf("Average Degree:             %.3f\n", lungGraphDegree)
	   fmt.Printf("Edge Density:               %.6f\n", lungEdgeDensity)
	   fmt.Printf("Positive Edges:             %d (%.3f%%)\n",
	       posLungEdges, 100*float64(posLungEdges)/float64(numLungEdges))
	   fmt.Printf("Negative Edges:             %d (%.3f%%)\n",
	       negLungEdges, 100*float64(negLungEdges)/float64(numLungEdges))
	   // Compute node degrees for each graph (for node-level distributions)
	   LungDegrees := ComputeDegrees(LungGraph)

	   //Build Louvain Graph type of cancer data to run the Louvain Algorithm
	   LungLouvainGraph := BuildLouvainGraph(LungPearsonCorrelationMatrix, LungThreshold)
	   fmt.Println("Running Louvain algorithm and assigning gene community clusters...")
	   //Run Louvain Algorithm from package to determine communities assignments (clusters of nodes) for further analysis and visualization
	   lungClusterMap, lungModularity := RunLouvainOnMatrix(LungLouvainGraph)
	*/
	////////////////////////////////////////////////////////////////
	//Module level Analyses per cancer type

	// Breast Cancer
	breastModuleSizes := ComputeModuleSizes(breastClusterMap)
	minB, maxB, meanB, medianB,
		numModsB, numSmallB,
		largestModB, largestPctB :=
		SummarizeModuleSizes(breastModuleSizes, numBreastNodes)

	fmt.Printf("\n===== %s Module Properties =====\n", dataset1Label)
	fmt.Printf("Number of modules:                       %d\n", numModsB)
	fmt.Printf("Module size (min / median / mean / max): %.0f / %.0f / %.1f / %.0f genes\n",
		minB, medianB, meanB, maxB)
	fmt.Printf("Modules with < 5 genes:                  %d\n", numSmallB)
	fmt.Printf("Largest module:                          ID %d with %.0f genes (%.2f%% of genes)\n",
		largestModB, maxB, largestPctB)
	// Access results for modularity
	fmt.Printf("%s graph modularity: %.4f\n", dataset1Label, breastModularity)

	// Per-community edge counts and densities
	breastModuleEdges := ComputeModuleEdgeCounts(BreastGraph, breastClusterMap)
	breastModuleDensities := ComputeModuleDensities(breastModuleSizes, breastModuleEdges)

	// per-community node count, edge count, density table
	fmt.Println("Module Structural Summary:")
	fmt.Println("ModuleID | NumNodes | NumEdges | Density")
	keys := make([]int, 0, len(breastModuleSizes))
	for k := range breastModuleSizes {
		keys = append(keys, k)
	}
	sort.Ints(keys)

	for _, commID := range keys {
		nNodes := breastModuleSizes[commID]
		nEdges := breastModuleEdges[commID]
		dens := breastModuleDensities[commID]
		fmt.Printf("%11d | %8d | %8d | %.4f\n", commID, nNodes, nEdges, dens)
	}

	/*
	   // Lung Cancer
	   LungModuleSizes := ComputeModuleSizes(lungClusterMap)
	   minL, maxL, meanL, medianL,
	       numModsL, numSmallL,
	       largestModL, largestPctL :=
	       SummarizeModuleSizes(LungModuleSizes, numLungNodes)

	   fmt.Printf("\n===== %s Module Properties =====\n", dataset3Label)
	   fmt.Printf("Number of modules:                       %d\n", numModsL)
	   fmt.Printf("Module size (min / median / mean / max): %.0f / %.0f / %.1f / %.0f genes\n",
	       minL, medianL, meanL, maxL)
	   fmt.Printf("Modules with < 5 genes:                  %d\n", numSmallL)
	   fmt.Printf("Largest module:                          ID %d with %.0f genes (%.2f%% of genes)\n",
	       largestModL, maxL, largestPctL)
	   // Access results for modularity
	   fmt.Printf("%s graph modularity: %.4f\n", dataset3Label, lungModularity)

	   // Per-community edge counts and densities
	   LungModuleEdges := ComputeModuleEdgeCounts(LungGraph, lungClusterMap)
	   LungModuleDensities := ComputeModuleDensities(LungModuleSizes, LungModuleEdges)

	   // per-community node count, edge count, density table
	   fmt.Println("Module Structural Summary:")
	   fmt.Println("ModuleID | NumNodes | NumEdges | Density")
	   keys := make([]int, 0, len(LungModuleSizes))
	   for k := range LungModuleSizes {
	       keys = append(keys, k)
	   }
	   sort.Ints(keys)

	   for _, commID := range keys {
	       nNodes := LungModuleSizes[commID]
	       nEdges := LungModuleEdges[commID]
	       dens := LungModuleDensities[commID]
	       fmt.Printf("%11d | %8d | %8d | %.4f\n", commID, nNodes, nEdges, dens)
	   }
	*/
	////////////////////////////////////////////////////////////////
	//Computing Clustering Coefficients for Graph Network
	fmt.Println("Computing Clustering Coefficients for Graph Network...")
	breastLocalC := LocalClusteringCoeff(BreastGraph)
	//LungLocalC := LocalClusteringCoeff(LungGraph)

	breastGlobalC := GlobalClusteringCoeff(breastLocalC)
	//LungGlobalC := GlobalClusteringCoeff(LungLocalC)
	fmt.Println("\n===== Global Clustering Coefficient =====")
	fmt.Printf("%s global clustering coefficient:  %.4f\n", dataset1Label, breastGlobalC)
	//fmt.Printf("%s global clustering coefficient:  %.4f\n", dataset3Label, LungGlobalC)

	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
	if arg == 2 {
		//Analyzing Dataset2

		//Begin Analyzing Dataset2 RNA seq gene count Data
		//Create a map of maps. Nested map mapping gene_name (string) to sampleID(string) to tpm value
		geneExpressionOvarianMap, err2 := ReadGeneExpressionDirToGeneMap(dataset2Path)

		// log any errors that may arise from the parsing/ reading in the ovarian data files
		if err2 != nil {
			log.Fatal(err2)
		}

		fmt.Printf("%s Gene Expression Data read and map created.\n", dataset2Label)

		/*fmt.Println("Printing first line of geneExpressionOvarianMap")
		  //Check the format of the data file by printing the first line of our geneExpressionOvarianMap
		  for key, val := range geneExpressionOvarianMap {
		      fmt.Println(key, val)
		      break
		  }*/

		fmt.Printf("Gathering sample names for %s data.\n", dataset2Label)
		// Store the sample names for the Ovarian Cancer Data
		ovarianSampleNames := GetSampleNames(geneExpressionOvarianMap)

		//Sort the ovarian sample names
		sortSampleNames(ovarianSampleNames)

		//Print the ovarian sample names to ensure they are sorted by numerical order
		//fmt.Println(ovarianSampleNames[:5])

		fmt.Println("Cleaning Data...Filtering out low expression genes...")
		// Use Mean-Based Filtering in order to normalize data while maintaining fixed number of samples across each gene
		filteredGeneExpressionOvarianMap := MeanBasedFilter(geneExpressionOvarianMap, ovarianSampleNames, 10.0)

		fmt.Printf("Gathering gene names for %s data.\n", dataset2Label)
		// Store the gene names for the Breast Cancer Data
		ovarianGeneNames := GetGeneNames(filteredGeneExpressionOvarianMap)

		//Print the first few breast gene names to check the names and ensure they are sorted
		//fmt.Println(ovarianGeneNames[:5])

		//Compute the Pearson Correlation on dataset2 data
		fmt.Printf("Computing Pearson Correlation between genes across all %s samples.\n", dataset2Label)

		OvarianPearsonCorrelationMatrix := ComputePearsonCorrelation(ovarianGeneNames, ovarianSampleNames, filteredGeneExpressionOvarianMap)

		fmt.Printf("Transforming the %s correlation matrix...\n", dataset2Label)
		ovarianQuantileSlice := TransformMatrixToSlice(OvarianPearsonCorrelationMatrix)

		//sort correlation values
		sortedOvarianQuantileSlice := SortCorrVals(ovarianQuantileSlice)

		fmt.Printf("Computing the %s quantile probabilities...\n", dataset2Label)
		computedOvarianQuantiles := ComputeQuantile(sortedOvarianQuantileSlice)

		fmt.Printf("%s Quantile Probabilities: %v\n", dataset2Label, computedOvarianQuantiles)

		//The 95% quantile probability chosen as the cutoff correlation value threshold

		ovarianThreshold := computedOvarianQuantiles[1]

		fmt.Println("Building the Graph Networks...")

		//Build weighted, undirected cancer graph network with edge weights corresponding to raw computed correlations between genes (nodes)

		OvarianGraph := BuildGraph(OvarianPearsonCorrelationMatrix, ovarianGeneNames, ovarianThreshold)

		////////////////////////////////////////////////////////////////
		//Analyzing Graph Properties (Network Level per cancer type)
		// Dataset2 (Ovarian Cancer)
		numOvarianNodes := len(OvarianGraph)
		numOvarianEdges := CalculateNumEdges(OvarianGraph)
		ovarianGraphDegree := ComputeAverageDegree(numOvarianEdges, numOvarianNodes)
		ovarianEdgeDensity := ComputeEdgeDensity(numOvarianEdges, numOvarianNodes)
		posOvarianEdges, negOvarianEdges := EdgeStats(OvarianGraph)
		fmt.Printf("\n===== %s Graph Properties =====\n", dataset2Label)
		fmt.Printf("Nodes (N):                  %d\n", numOvarianNodes)
		fmt.Printf("Edges (E):                  %d\n", numOvarianEdges)
		fmt.Printf("Average Degree:             %.3f\n", ovarianGraphDegree)
		fmt.Printf("Edge Density:               %.6f\n", ovarianEdgeDensity)
		fmt.Printf("Positive Edges:             %d (%.3f%%)\n",
			posOvarianEdges, 100*float64(posOvarianEdges)/float64(numOvarianEdges))
		fmt.Printf("Negative Edges:             %d (%.3f%%)\n",
			negOvarianEdges, 100*float64(negOvarianEdges)/float64(numOvarianEdges))
		// Compute node degrees for each graph (for node-level distributions)
		ovarianDegrees := ComputeDegrees(OvarianGraph)

		//Build Louvain Graph type of cancer data to run the Louvain Algorithm
		ovarianLouvainGraph := BuildLouvainGraph(OvarianPearsonCorrelationMatrix, ovarianThreshold)

		fmt.Println("Running Louvain algorithm and assigning gene community clusters...")
		//Run Louvain Algorithm from package to determine communities assignments (clusters of nodes) for further analysis and visualization

		ovarianClusterMap, ovarianModularity := RunLouvainOnMatrix(ovarianLouvainGraph)

		////////////////////////////////////////////////////////////////
		// Module level Analyses per cancer type
		// Ovarian Cancer
		ovarianModuleSizes := ComputeModuleSizes(ovarianClusterMap)
		minO, maxO, meanO, medianO,
			numModsO, numSmallO,
			largestModO, largestPctO :=
			SummarizeModuleSizes(ovarianModuleSizes, numOvarianNodes)

		fmt.Printf("\n===== %s Module Properties =====\n", dataset2Label)
		fmt.Printf("Number of modules:                       %d\n", numModsO)
		fmt.Printf("Module size (min / median / mean / max): %.0f / %.0f / %.1f / %.0f genes\n",
			minO, medianO, meanO, maxO)
		fmt.Printf("Modules with < 5 genes:                  %d\n", numSmallO)
		fmt.Printf("Largest module:                          ID %d with %.0f genes (%.2f%% of genes)\n",
			largestModO, maxO, largestPctO)
		// Access results for modularity
		fmt.Printf("%s graph modularity: %.4f\n", dataset2Label, ovarianModularity)

		// Per-community edge counts and densities (Ovarian)
		ovarianModuleEdges := ComputeModuleEdgeCounts(OvarianGraph, ovarianClusterMap)
		ovarianModuleDensities := ComputeModuleDensities(ovarianModuleSizes, ovarianModuleEdges)

		// per-community node count, edge count, density table
		fmt.Println("Module Structural Summary:")
		fmt.Println("ModuleID | NumNodes | NumEdges | Density")
		keys = make([]int, 0, len(ovarianModuleSizes))
		for k := range ovarianModuleSizes {
			keys = append(keys, k)
		}
		sort.Ints(keys)

		for _, commID := range keys {
			nNodes := ovarianModuleSizes[commID]
			nEdges := ovarianModuleEdges[commID]
			dens := ovarianModuleDensities[commID]
			fmt.Printf("%11d | %8d | %8d | %.4f\n", commID, nNodes, nEdges, dens)
		}

		////////////////////////////////////////////////////////////////
		//Computing Clustering Coefficients for each Graph Network
		fmt.Println("Computing Clustering Coefficients for Graph Network...")
		ovarianLocalC := LocalClusteringCoeff(OvarianGraph)
		ovarianGlobalC := GlobalClusteringCoeff(ovarianLocalC)

		fmt.Println("\n===== Global Clustering Coefficient =====")
		fmt.Printf("%s global clustering coefficient: %.4f\n", dataset2Label, ovarianGlobalC)

		////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////
		// Network Topology Analysis: Real vs Random Graph Comparisons

		fmt.Println("\nNetwork Topology Analysis: Analyzing whether gene co-expression networks exhibit non-random topological properties...")

		//Setting the seed
		rand.Seed(42)

		// Generate random graphs using Erdős-Rényi G(n,p) model for both cancer types
		fmt.Println("Generating breast cancer random graph using Erdős-Rényi G(n,p) model...")
		breastRandomGraph := RandomGraphGenerator(breastGeneNames, 0.3)

		fmt.Println("Generating ovarian cancer random graph using Erdős-Rényi G(n,p) model...")
		ovarianRandomGraph := RandomGraphGenerator(ovarianGeneNames, 0.3)

		// Compute degree distributions for random graphs
		fmt.Println("Computing degree distributions for random graph networks...")
		breastRandomDegrees := ComputeDegrees(breastRandomGraph)
		ovarianRandomDegrees := ComputeDegrees(ovarianRandomGraph)

		// Statistical comparison: Real vs Random degree distributions using Kolmogorov-Smirnov test
		fmt.Printf("\n===== %s: Real vs Random Network Comparison =====\n", dataset1Label)
		breastDStatistic, breastPValue := KSTest(breastDegrees, breastRandomDegrees)
		fmt.Printf("Kolmogorov-Smirnov test comparing real vs random degree distributions:\n")
		fmt.Printf("D-statistic: %.4f\n", breastDStatistic)
		// Handle very small p-values for scientific reporting
		if breastPValue < 0.001 {
			fmt.Printf("P-value: < 0.001 (highly significant)\n")
		} else {
			fmt.Printf("P-value: %.3f\n", breastPValue)
		}

		// Interpret statistical significance for breast cancer
		if breastPValue < 0.001 {
			fmt.Printf("Result: Significantly different degree distributions (p < 0.001)\n")
			fmt.Printf("The real %s network has a significantly different structure than the random model.\n", dataset1Label)
		} else {
			fmt.Printf("Result: No significant difference in degree distributions (p ≥ 0.001)\n")
			fmt.Printf("The real %s network structure is not significantly different from the random model.\n", dataset1Label)
		}

		fmt.Printf("\n===== %s: Real vs Random Network Comparison =====\n", dataset2Label)
		ovarianDStatistic, ovarianPValue := KSTest(ovarianDegrees, ovarianRandomDegrees)
		fmt.Printf("Kolmogorov-Smirnov test comparing real vs random degree distributions:\n")
		fmt.Printf("D-statistic: %.4f\n", ovarianDStatistic)
		// Handle very small p-values for scientific reporting
		if ovarianPValue < 0.001 {
			fmt.Printf("P-value: < 0.001 (highly significant)\n")
		} else {
			fmt.Printf("P-value: %.3f\n", ovarianPValue)
		}

		// Interpret statistical significance for ovarian cancer
		if ovarianPValue < 0.001 {
			fmt.Printf("Result: Significantly different degree distributions (p < 0.001)\n")
			fmt.Printf("The real %s network has a significantly different structure than the random model.\n", dataset2Label)
		} else {
			fmt.Printf("Result: No significant difference in degree distributions (p ≥ 0.001)\n")
			fmt.Printf("The real %s network structure is not significantly different from the random model.\n", dataset2Label)
		}

		// Cross-Cancer Type Analysis: Breast vs Ovarian Degree Distribution Comparison
		fmt.Printf("\n===== Degree Distribution Cross-Cancer Analysis: %s vs %s =====\n", dataset1Label, dataset2Label)

		// Statistical comparison: Breast vs Ovarian degree distributions using Kolmogorov-Smirnov test
		crossCancerDStatistic, crossCancerPValue := KSTest(breastDegrees, ovarianDegrees)
		fmt.Printf("Kolmogorov-Smirnov test comparing %s vs %s degree distributions:\n", dataset1Label, dataset2Label)
		fmt.Printf("D-statistic: %.4f\n", crossCancerDStatistic)
		// Handle very small p-values for scientific reporting
		if crossCancerPValue < 0.001 {
			fmt.Printf("P-value: < 0.001 (highly significant)\n")
		} else {
			fmt.Printf("P-value: %.3f\n", crossCancerPValue)
		}

		// Interpret statistical significance for cross-cancer comparison
		if crossCancerPValue < 0.001 {
			fmt.Printf("Result: Significantly different degree distributions between cancer types (p < 0.001)\n")
			fmt.Printf("The %s and %s networks have significantly different topological structures.\n", dataset1Label, dataset2Label)
		} else {
			fmt.Printf("Result: No significant difference in degree distributions between cancer types (p ≥ 0.001)\n")
			fmt.Printf("The %s and %s networks have similar topological structures.\n", dataset1Label, dataset2Label)
		}

		////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////
		// Local Clustering Coefficient Comparison between Cancer Types
		fmt.Printf("\n===== Clustering Coefficient Cross-Cancer Analysis: %s vs %s =====\n", dataset1Label, dataset2Label)

		// Filter out NaN values from clustering coefficients before statistical comparison
		breastClusteringFiltered := FilterNaNs(breastLocalC)
		ovarianClusteringFiltered := FilterNaNs(ovarianLocalC)

		// Statistical comparison: Breast vs Ovarian clustering coefficient distributions using Kolmogorov-Smirnov test
		clusteringDStatistic, clusteringPValue := KSTest(breastClusteringFiltered, ovarianClusteringFiltered)
		fmt.Printf("Kolmogorov-Smirnov test comparing %s vs %s clustering coefficient distributions:\n", dataset1Label, dataset2Label)
		fmt.Printf("D-statistic: %.4f\n", clusteringDStatistic)
		// Handle very small p-values for scientific reporting
		if clusteringPValue < 0.001 {
			fmt.Printf("P-value: < 0.001 (highly significant)\n")
		} else {
			fmt.Printf("P-value: %.3f\n", clusteringPValue)
		}

		// Interpret statistical significance for clustering coefficient comparison
		if clusteringPValue < 0.001 {
			fmt.Printf("Result: Significantly different clustering coefficient distributions between cancer types (p < 0.001)\n")
			fmt.Printf("The %s and %s networks have significantly different local clustering patterns.\n", dataset1Label, dataset2Label)
		} else {
			fmt.Printf("Result: No significant difference in clustering coefficient distributions between cancer types (p ≥ 0.001)\n")
			fmt.Printf("The %s and %s networks have similar local clustering patterns.\n", dataset1Label, dataset2Label)
		}

		////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////
		// Data Export: Write Network Analysis Results to CSV Files
		// Export network structures, community assignments, and statistical measures
		// for visualization and further analysis in R/Shiny applications

		fmt.Printf("Writing %s gene network information to CSV files...\n", dataset2Label)

		err3 := WriteCommunitiesCSV(outputDir+"/dataset2_nodes_communities.csv", ovarianGeneNames, ovarianClusterMap)

		if err3 != nil {
			panic(err3)
		}

		err4 := WriteEdgesCSV(outputDir+"/dataset2_edges.csv", OvarianGraph)

		if err4 != nil {
			panic(err4)
		}

		// Write per-community stats for Ovarian Cancer to CSV (for Shiny + R distributions)
		errOCStats := WriteCommunityStatsCSV(outputDir+"/dataset2_community_stats.csv", ovarianModuleSizes, ovarianModuleEdges, ovarianModuleDensities)

		if errOCStats != nil {
			panic(errOCStats)
		}

		// Write node-level stats to CSV for Ovarian Cancer
		errONode := WriteNodeStatsCSV(outputDir+"/dataset2_node_stats.csv", ovarianGeneNames, ovarianClusterMap, ovarianDegrees, ovarianLocalC)

		if errONode != nil {
			panic(errONode)
		}

		/////////////// Largest Module Heatmap ////////////
		// grabs largest module and saves data in .csv file
		SaveLargestModule(OvarianGraph, ovarianClusterMap, "largest_module_ovarian.csv")

		///save Correlation Matrix in .csv/////////
		errCorrMatrixO := SaveMatrixAsCSV("ExtraData/OvarianPearsonCorrelationMatrix.csv",
			OvarianPearsonCorrelationMatrix, ovarianGeneNames, ovarianGeneNames)

		if errCorrMatrixO != nil {
			fmt.Println("Error saving ovarian Pearson matrix:", errCorrMatrixO)
		} else {
			fmt.Println("Ovarian Pearson correlation matrix saved.")
		}

		///////// Save OVARIAN module within 400-550 genes ///////////
		inputCSVO := "ShinyApp/dataset2_nodes_communities.csv"
		outputDirO := "ExtraData"
		minSize := 400
		maxSize := 550

		if errModuleInRangeO := SaveModulesInSizeRange(inputCSVO, outputDirO, minSize, maxSize); errModuleInRangeO != nil {
			log.Fatalf("Error: %v", errModuleInRangeO)
		}

		///get top 200 genes///
		errSTopNodesO := SaveTopNGenes(OvarianGraph, 100, "top_100_genes_ovarian.csv")
		if errSTopNodesO != nil {
			log.Fatalf("Error: %v", errSTopNodesO)
		}

		// Call the subroutine to save edge counts for top 200 genes
		errCEdgesO := SaveEdgeCounts(OvarianGraph, 100, "top100_edge_counts_ovarian.csv")
		if errCEdgesO != nil {
			fmt.Println("Error:", errCEdgesO)
		}
	}

	fmt.Printf("Writing %s gene network information to CSV files...\n", dataset1Label)

	err5 := WriteCommunitiesCSV(outputDir+"/dataset1_nodes_communities.csv", breastGeneNames, breastClusterMap)

	if err5 != nil {
		panic(err5)
	}

	err6 := WriteEdgesCSV(outputDir+"/dataset1_edges.csv", BreastGraph)

	if err6 != nil {
		panic(err6)
	}

	// Write per-community stats for Breast Cancer to CSV (for Shiny + R distributions)
	errBCStats := WriteCommunityStatsCSV(outputDir+"/dataset1_community_stats.csv", breastModuleSizes, breastModuleEdges, breastModuleDensities)

	if errBCStats != nil {
		panic(errBCStats)
	}

	// Write node-level stats to CSV for Breast Cancer
	errBNode := WriteNodeStatsCSV(outputDir+"/dataset1_node_stats.csv", breastGeneNames, breastClusterMap, breastDegrees, breastLocalC)

	if errBNode != nil {
		panic(errBNode)
	}
	fmt.Printf("Writing %s gene network information to CSV files...\n", dataset1Label)

	// Write per-community stats for Breast Cancer to CSV (for Shiny + R distributions)
	errBCStats = WriteCommunityStatsCSV(outputDir+"/dataset1_community_stats.csv", breastModuleSizes, breastModuleEdges, breastModuleDensities)

	// Write per-community stats for Lung Cancer to CSV (for Shiny + R distributions)
	//errBCStats := WriteCommunityStatsCSV(outputDir+"/dataset3_community_stats.csv", LungModuleSizes, LungModuleEdges, LungModuleDensities)

	if errBCStats != nil {
		panic(errBCStats)
	}

	// Write node-level stats to CSV for Breast Cancer
	errBNode = WriteNodeStatsCSV(outputDir+"/dataset1_node_stats.csv", breastGeneNames, breastClusterMap, breastDegrees, breastLocalC)

	// Write node-level stats to CSV for Breast Cancer
	//errBNode := WriteNodeStatsCSV(outputDir+"/dataset3_node_stats.csv", lungGeneNames, lungClusterMap, LungDegrees, LungLocalC)

	if errBNode != nil {
		panic(errBNode)
	}

	/////////////// Largest Module Heatmap ////////////
	// grabs largest module and saves data in .csv file
	SaveLargestModule(BreastGraph, breastClusterMap, "largest_module.csv")

	/*
		///save Correlation Matrix in .csv/////////
		//breast (dataset1)
		errCorrMatrixB := SaveMatrixAsCSV("BreastPearsonCorrelationMatrix.csv",
			BreastPearsonCorrelationMatrix, breastGeneNames, breastGeneNames)

		if errCorrMatrixB != nil {
			fmt.Println("Error saving Pearson matrix:", errCorrMatrixB)
		} else {
			fmt.Println("Pearson correlation matrix saved.")
		}
	*/

	//////// Save module within 400-550 genes ///////////
	// Breast - Read node–community CSV
	inputCSV := "ShinyApp/dataset1_nodes_communities.csv"
	outputDir := "ExtraData"
	minSize := 350
	maxSize := 550

	if errModuleInRange := SaveModulesInSizeRange(inputCSV, outputDir, minSize, maxSize); errModuleInRange != nil {
		log.Fatalf("Error: %v", errModuleInRange)
	}

	//return top 200 genes based on degree
	errSTopNodes := SaveTopNGenes(BreastGraph, 100, "top_100_genes_breast.csv")
	if errSTopNodes != nil {
		log.Fatalf("Error: %v", errSTopNodes)
	}

	// Call the subroutine to save edge counts for top 200 genes
	errCEdges := SaveEdgeCounts(BreastGraph, 100, "top100_edge_counts_breast.csv")
	if errCEdges != nil {
		fmt.Println("Error:", errCEdges)
	}

}
