//Noemi Banda, Emma Franken, Beth Vazquez Smith
//Programming for Scientists Project
//October 31, 2025

package main

// map of strings to a map of strings of floats
type ExpressionMap map[string]map[string]float64

// 2D slice of float64
type CorrelationMatrix [][]float64

// Tree is a slice of pointers to nodes
type GraphNetwork []*Node

// Node is an object that represents a node of a tree.
type Node struct {
	ID       int
	GeneName string
	Edges    []*Edge
}

type Edge struct {
	To     *Node
	Weight float64 // correlation
}
