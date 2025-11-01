//Noemi Banda, Emma Franken, Beth Vazquez Smith
//Programming for Scientists Project
//October 31, 2025

package main

type ExpressionMatrix [][]float64 //or would this be map[string]float64

// using Network instead of tree
// Network is slice of pointers to nodes
type Network []*Node

type Node struct {
	ExpressionVal  *float64  //strength of co-expression
	Weights        []float64 //need weight for weighted matrix
	GeneName       string
	Child1, Child2 *Node //need to figure out a way to not have finite num of children (i think its Children []*Node)
}
