//Noemi Banda, Emma Franken, Beth Vazquez Smith
//Programming for Scientists Project
//October 31, 2025

package main

type ExpressionMatrix map[string]float64 //upgma was [][]float64

// using Network instead of tree
// Network is slice of pointers to nodes
type Network []*Node

type Node struct {
	ExpressionVal *float64  //strength of co-expression
	Weights       []float64 //need weight for weighted matrix
	GeneName      string
	Children      []*Node //need to figure out a way to not have finite num of children (i think its Children []*Node)
}
