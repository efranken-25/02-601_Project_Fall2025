//This file has all of the helper functions for the testing code.

package main

import (
	"02-601_Project_Fall2025/louvain"
	"bufio"
	"math"
	"os"
	"strconv"
	"strings"
)

func parseFloatList(text string) []float64 {
	lines := strings.Split(text, "\n")
	vals := []float64{}

	for _, line := range lines {
		line = strings.TrimSpace(line)
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		parts := strings.Split(line, ",")
		for _, p := range parts {
			p = strings.TrimSpace(p)
			if p == "" {
				continue
			}

			v, err := strconv.ParseFloat(p, 64)
			if err != nil {
				panic(err)
			}
			vals = append(vals, v)
		}
	}

	return vals
}

// Compare two matrices, treating NaNs as equal
func floatMatrixEqual(a, b [][]float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if len(a[i]) != len(b[i]) {
			return false
		}
		for j := range a[i] {
			if math.Abs(a[i][j]-b[i][j]) > tol {
				return false
			}
		}
	}
	return true
}

// Helper: compare nested float maps
func nestedFloatMapsEqual(a, b ExpressionMap) bool {
	if len(a) != len(b) {
		return false
	}
	for g, rowA := range a {
		rowB, ok := b[g]
		if !ok || len(rowA) != len(rowB) {
			return false
		}
		for s, valA := range rowA {
			valB, ok := rowB[s]
			if !ok || valA != valB {
				return false
			}
		}
	}
	return true
}

// Helper: compare string slices with order
func stringSlicesEqual(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// Helper: compare string slices ignoring order
func stringSlicesEqualUnordered(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	freq := make(map[string]int)
	for _, s := range a {
		freq[s]++
	}
	for _, s := range b {
		freq[s]--
		if freq[s] < 0 {
			return false
		}
	}
	return true
}

// Helper: read files from a directory
func ReadDirectory(dir string) []os.DirEntry {
	files, err := os.ReadDir(dir)
	if err != nil {
		panic(err)
	}
	return files
}

// floatSlicesEqualUnordered checks if two float64 slices contain the same elements,
// ignoring order. Uses a tolerance for floating-point comparison.
func floatSlicesEqualUnordered(a, b []float64) bool {
	if len(a) != len(b) {
		return false
	}

	// Number of decimals to round to (controls tolerance)
	decimals := 9

	countMap := make(map[float64]int)

	// Count elements in a
	for _, val := range a {
		key := roundFloat(val, decimals)
		countMap[key]++
	}

	// Subtract counts using elements from b
	for _, val := range b {
		key := roundFloat(val, decimals)
		countMap[key]--
		if countMap[key] < 0 {
			return false
		}
	}

	return true
}

// roundFloat rounds a float64 to the specified number of decimal places
func roundFloat(f float64, decimals int) float64 {
	factor := math.Pow(10, float64(decimals))
	return math.Round(f*factor) / factor
}

func floatSlicesApproxEqual(a, b []float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if math.Abs(a[i]-b[i]) > tol {
			return false
		}
	}
	return true
}

// floatSlicesEqualOrdered checks if two slices are equal element by element,
// using rounding to avoid floating-point precision issues.
func floatSlicesEqualOrdered(a, b []float64) bool {
	if len(a) != len(b) {
		return false
	}

	decimals := 9 // number of decimals to round to for comparison
	for i := range a {
		if roundFloat(a[i], decimals) != roundFloat(b[i], decimals) {
			return false
		}
	}
	return true
}

func intSlicesEqualOrdered(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func extractEdgesFromGraph(g *louvain.Graph) map[[2]int]float64 {
	edges := make(map[[2]int]float64)

	for u := 0; u < g.NumNodes; u++ {
		for idx, v := range g.Adjacency[u] {
			if u < v { // avoid duplicates
				edges[[2]int{u, v}] = g.Weights[u][idx]
			}
		}
	}
	return edges
}

func parseEdges(r *strings.Reader) [][3]float64 {
	var edges [][3]float64
	scanner := bufio.NewScanner(r)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		fields := strings.Fields(line)
		if len(fields) != 3 {
			continue // skip invalid lines or empty output
		}

		u, err1 := strconv.Atoi(fields[0])
		v, err2 := strconv.Atoi(fields[1])
		w, err3 := strconv.ParseFloat(fields[2], 64)
		if err1 != nil || err2 != nil || err3 != nil {
			panic("invalid edge line: " + line)
		}

		edges = append(edges, [3]float64{float64(u), float64(v), w})
	}

	return edges
}

func compareCommunities(a, b map[int]int) bool {
	if len(a) != len(b) {
		return false
	}

	// Build mapping from a IDs to b IDs
	labelMap := make(map[int]int)
	for node, cidA := range a {
		cidB := b[node]
		if mapped, ok := labelMap[cidA]; ok {
			if mapped != cidB {
				return false
			}
		} else {
			labelMap[cidA] = cidB
		}
	}
	return true
}

// parseGraphNetwork parses input string into GraphNetwork
// Expected input format: one node per line, e.g.:
// 0 G1 1:0.5,2:-0.3
func parseGraphNetwork(input string) GraphNetwork {
	lines := strings.Split(strings.TrimSpace(input), "\n")
	graph := make(GraphNetwork, len(lines))

	for _, line := range lines {
		parts := strings.Fields(line)
		id, _ := strconv.Atoi(parts[0])
		name := parts[1]

		node := &Node{
			ID:       id,
			GeneName: name,
			Edges:    []*Edge{},
		}

		if len(parts) > 2 {
			edgesStr := strings.Split(parts[2], ",")
			for _, e := range edgesStr {
				edgeParts := strings.Split(e, ":")
				toID, _ := strconv.Atoi(edgeParts[0])
				weight, _ := strconv.ParseFloat(edgeParts[1], 64)
				edge := &Edge{To: &Node{ID: toID}, Weight: weight}
				node.Edges = append(node.Edges, edge)
			}
		}
		graph[id] = node
	}
	return graph
}

// parseEdgeStatsOutput parses output string into pos and neg integers
// Expected output format: two integers separated by space, e.g.:
// 3 1
func parseEdgeStatsOutput(output string) (int, int) {
	parts := strings.Fields(strings.TrimSpace(output))
	pos, _ := strconv.Atoi(parts[0])
	neg, _ := strconv.Atoi(parts[1])
	return pos, neg
}

// Helper: parse cluster map from file
func parseClusterMap(data string) map[int]int {
	lines := strings.Split(strings.TrimSpace(data), "\n")
	clusterMap := make(map[int]int)
	for _, line := range lines {
		parts := strings.Fields(line) // "nodeID communityID"
		if len(parts) != 2 {
			continue
		}
		nodeID, _ := strconv.Atoi(parts[0])
		communityID, _ := strconv.Atoi(parts[1])
		clusterMap[nodeID] = communityID
	}
	return clusterMap
}

// Helper: parse expected module sizes from file
func parseModuleSizes(data string) map[int]int {
	lines := strings.Split(strings.TrimSpace(data), "\n")
	sizes := make(map[int]int)
	for _, line := range lines {
		parts := strings.Fields(line) // "communityID size"
		if len(parts) != 2 {
			continue
		}
		communityID, _ := strconv.Atoi(parts[0])
		size, _ := strconv.Atoi(parts[1])
		sizes[communityID] = size
	}
	return sizes
}

func parseStringList(s string) []string {
	parts := strings.Split(strings.TrimSpace(s), ",")
	result := make([]string, 0, len(parts))
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p != "" {
			result = append(result, p)
		}
	}
	return result
}
