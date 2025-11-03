//Noemi Banda, Emma Franken, Beth Vazquez Smith
//Programming for Scientists Project
//October 31, 2025

package main

type ExpressionMap map[string]map[string]float64 //upgma was [][]float64

type ExpressionMatrix [][]float64

// using Network instead of tree
// Network is slice of pointers to nodes
type Network []*Node

type Edge struct {
	Nodes *Node   //node this edge connects to
	Corr  float64 //weight for weighted matrix
}

type Node struct {
	GeneName string
	Edges    []*Edge //all edges from this node and to others
}
