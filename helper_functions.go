//This file has all of the helper functions for the testing code.

package main

import (
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
