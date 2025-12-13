// Gene Co-Expression Network Analysis Pipeline
// Programming for Scientists (Group 2)
// Date: December 12, 2025

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
