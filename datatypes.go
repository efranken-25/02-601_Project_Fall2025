//Noemi Banda, Emma Franken, Beth Vazquez Smith
//Programming for Scientists Project
//October 31, 2025

package main

type ExpressionMap map[string]map[string]float64

type CorrelationMatrix [][]float64

// Graph network is a slice of pointers to nodes
type GraphNetwork []*Node

// Node is an object that represents a node of a graph network.
type Node struct {
	ID       int
	GeneName string
	Edges    []*Edge
}

// Edge is an object representing an edge between nodes of a graph network
type Edge struct {
	To     *Node
	Weight float64 // correlation
}

type ModuleMatch struct {
	SourceModuleID int     // module that you are analyzing
	TargetModuleID int     // the module it matches best in the other cancer type
	JaccardValue   float64 // the max Jaccard value for this pair
}

// ModuleOverlapStats stores all stats for one breast–ovarian module pair
type ModuleOverlapStats struct {
	BreastModuleID  int
	OvarianModuleID int

	SizeA   int // |A| = size of breast module
	SizeB   int // |B| = size of ovarian module
	Overlap int // |A ∩ B|

	// 2×2 contingency table counts
	A int // in A and in B
	B int // not in A, but in B
	C int // in A, but not in B
	D int // not in A, not in B

	TotalGenes int     // N
	Jaccard    float64 // Jaccard index for this pair
}

// Holds final results with p- and q-values
type ModuleOverlapResult struct {
	ModuleOverlapStats

	PValue float64
	QValue float64
}

/*
type idxP struct {
	idx int
	p   float64
}
*/
