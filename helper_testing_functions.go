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

// ParseFloatList takes a string as input and returns a slice of float64 parsed from that string.
// It also ignores the empty lines and comment lines starting with #
func ParseFloatList(text string) []float64 {
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

// FloatMatrixEqual checks whether two 2-D slices of float64 are equal within a given tolerance.
// Returns a boolean variable corresponding to whether or not the matrices are equal.
func FloatMatrixEqual(a, b [][]float64, tol float64) bool {
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

// NestedFloatMapsEqual compares two nested maps of floats for exact equality.
// Returns a boolean variable corresponding to whether or not the maps are equal.
func NestedFloatMapsEqual(a, b ExpressionMap) bool {
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

// StringSlicesEqual checks whether two slices of strings are exactly equal.
// Returns a boolean variable corresponding to whether or not the slices are equal.
func StringSlicesEqual(a, b []string) bool {
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

// StringsSlicesEqualUnordered checks whether 2 slices of strings contain the same elements with the same multiplicities, and ignores the order.
// Returns a boolean variable corresponding to whether or not the slices are equal.
func StringSlicesEqualUnordered(a, b []string) bool {
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

// ReadDirectory reads all entries in a directory and returns them as a slice of os.DirEntry.
func ReadDirectory(dir string) []os.DirEntry {
	files, err := os.ReadDir(dir)
	if err != nil {
		panic(err)
	}
	return files
}

// FloatSlicesEqualUnordered checks if two float64 slices contain the same elements,
// Ignores order. Uses a tolerance for floating-point comparison.
func FloatSlicesEqualUnordered(a, b []float64) bool {
	if len(a) != len(b) {
		return false
	}

	// Number of decimals to round to (controls tolerance)
	decimals := 9

	countMap := make(map[float64]int)

	// Count elements in a
	for _, val := range a {
		key := RoundFloat(val, decimals)
		countMap[key]++
	}

	// Subtract counts using elements from b
	for _, val := range b {
		key := RoundFloat(val, decimals)
		countMap[key]--
		if countMap[key] < 0 {
			return false
		}
	}

	return true
}

// RoundFloat rounds a float64 to the specified number of decimal places
func RoundFloat(f float64, decimals int) float64 {
	factor := math.Pow(10, float64(decimals))
	return math.Round(f*factor) / factor
}

// FloatSlicesApproxEqual checks whether two slices of float64 are approx equal element by element, using a tolerance.
// Returns a boolean variable corresponding to whether or not the slices are equal.
func FloatSlicesApproxEqual(a, b []float64, tol float64) bool {
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

// FloatSlicesEqualOrdered checks if two slices are equal element by element,
// using rounding to avoid floating-point precision issues.
func FloatSlicesEqualOrdered(a, b []float64) bool {
	if len(a) != len(b) {
		return false
	}

	decimals := 9 // number of decimals to round to for comparison
	for i := range a {
		if RoundFloat(a[i], decimals) != RoundFloat(b[i], decimals) {
			return false
		}
	}
	return true
}

// IntSlicesEqualOrdered checks whether 2 []int slices are exactly equal, element-by-element, in the same order.
// Returns a boolean variable corresponding to whether or not the slices are equal.
func IntSlicesEqualOrdered(a, b []int) bool {
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

// ExtractEdgesFromGraph extracts all unique edges from a louvain.Graph and returns them in a map where the key is a 2-element array representing an undirected edge and the value is the weight of that edge
func ExtractEdgesFromGraph(g *louvain.Graph) map[[2]int]float64 {
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

// ParseEdges reads a list of edges from a text source and returns them as a slice of [3]float64.
func ParseEdges(r *strings.Reader) [][3]float64 {
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

// CompareCommunities checks whether 2 community assignments (a and b) of nodes are structurally equivalent, ignoring the actual numeric labels.
// Returns a boolean variable corresponding to whether or not the maps are equal.
func CompareCommunities(a, b map[int]int) bool {
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

// ParseGraphNetwork parses input string into GraphNetwork
// Expected input format: one node per line, e.g.: 0 G1 1:0.5,2:-0.3
func ParseGraphNetwork(input string) GraphNetwork {
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

// ParseEdgeStatsOutput parses output string into pos and neg integers
// Expected output format: two integers separated by space, e.g.:
// 3 1
func ParseEdgeStatsOutput(output string) (int, int) {
	parts := strings.Fields(strings.TrimSpace(output))
	pos, _ := strconv.Atoi(parts[0])
	neg, _ := strconv.Atoi(parts[1])
	return pos, neg
}

// ParseClusterMap reads a textual representation of a community assingment and returns it as a map.
func ParseClusterMap(data string) map[int]int {
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

// ParseModuleSizes reads a textual representation of module sizes and returns it as a map
func ParseModuleSizes(data string) map[int]int {
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

// ParseStringList takes a string contain comma-separated values and returns a slice of trimmed, non-empty strings.
func ParseStringList(s string) []string {
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

// Compare two maps[int][]string for equality (ignoring order)
func CompareModuleMaps(a, b map[int][]string) bool {
	if len(a) != len(b) {
		return false
	}
	for k, va := range a {
		vb, ok := b[k]
		if !ok || len(va) != len(vb) {
			return false
		}
		m := make(map[string]int)
		for _, s := range va {
			m[s]++
		}
		for _, s := range vb {
			if _, ok := m[s]; !ok {
				return false
			}
			m[s]--
			if m[s] == 0 {
				delete(m, s)
			}
		}
		if len(m) != 0 {
			return false
		}
	}
	return true
}

func AdjacencySlicesEqual(a, b [][]int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if len(a[i]) != len(b[i]) {
			return false
		}
		counts := make(map[int]int)
		for _, val := range a[i] {
			counts[val]++
		}
		for _, val := range b[i] {
			if counts[val] == 0 {
				return false
			}
			counts[val]--
		}
		for _, v := range counts {
			if v != 0 {
				return false
			}
		}
	}
	return true
}
