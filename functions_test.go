// Gene Co-Expression Network Analysis Pipeline
// Programming for Scientists (Group 2)
// Date: December 12, 2025
package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

type TransformMapToSliceTest struct {
	geneName                string
	geneExpressionBreastMap ExpressionMap
	result                  []float64
}

// TestTransformMapToSlice tests the TransformMapToSlice function.
func TestTransformMapToSlice(t *testing.T) {
	tests := ReadTransformMapToSliceTests("Tests/TransformMapToSlice")
	for _, test := range tests {
		got := TransformMapToSlice(test.geneName, test.geneExpressionBreastMap)
		if !FloatSlicesEqualUnordered(got, test.result) {
			t.Errorf("TransformMapToSlice(%s) failed.\nGot: %v\nWant: %v",
				test.geneName, got, test.result)
		}
	}
}

// ReadTransformMapToSliceTests reads test input and output files for TransformMapToSlice function,
// parsing gene expression data and expected results from the test directory structure.
// It takes a directory path as input and returns a slice of TransformMapToSliceTest structs.
func ReadTransformMapToSliceTests(directory string) []TransformMapToSliceTest {
	// Read all input and output files from test directories
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")
	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	// Initialize slice to store all test cases
	tests := make([]TransformMapToSliceTest, len(inputFiles))

	// Parse each input file to extract gene expression data
	for i, inputFile := range inputFiles {
		f, err := os.Open(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)

		// Initialize test structure for current file
		current := TransformMapToSliceTest{
			geneExpressionBreastMap: make(map[string]map[string]float64),
		}
		var currentMapKey string

		// Parse each line to build gene expression map
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			// Skip empty lines and comments
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}

			switch {
			// Parse gene name and initialize its expression map
			case strings.HasPrefix(line, "geneName:"):
				current.geneName = strings.TrimSpace(strings.TrimPrefix(line, "geneName:"))
				current.geneExpressionBreastMap[current.geneName] = make(map[string]float64)
				currentMapKey = current.geneName

			default:
				if currentMapKey != "" {
					pair := strings.SplitN(line, ":", 2)
					if len(pair) != 2 {
						continue
					}
					key := strings.TrimSpace(pair[0])
					valStr := strings.TrimSpace(pair[1])
					if valStr == "" {
						continue // skip lines like `result:`
					}
					val, err := strconv.ParseFloat(valStr, 64)
					if err != nil {
						panic(err)
					}
					current.geneExpressionBreastMap[currentMapKey][key] = val
				}
			}
		}
		f.Close()
		tests[i] = current
	}

	// Parse output files to get expected results
	for i, outputFile := range outputFiles {
		f, err := os.Open(directory + "/Output/" + outputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)
		// Collect expected float values from output file
		outputVals := []float64{}
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			// Skip empty lines and comments
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			for _, p := range strings.Split(line, ",") {
				valStr := strings.TrimSpace(p)
				if valStr == "" {
					continue
				}
				val, err := strconv.ParseFloat(valStr, 64)
				if err != nil {
					panic(err)
				}
				outputVals = append(outputVals, val)
			}
		}
		f.Close()
		tests[i].result = outputVals
	}

	return tests
}

type GetSampleNamesTest struct {
	geneExpressionBreastMap ExpressionMap
	result                  []string
}

// TestGetSampleNames tests the GetSampleNames function.
func TestGetSampleNames(t *testing.T) {
	tests := ReadGetSampleNamesTests("Tests/GetSampleNames")
	for _, test := range tests {
		got := GetSampleNames(test.geneExpressionBreastMap)
		if !StringSlicesEqualUnordered(got, test.result) {
			t.Errorf("GetSampleNames() failed.\nGot: %v\nWant: %v",
				got, test.result)
		}
	}
}

// ReadGetSampleNamesTests loads test cases for GetSampleNames function by parsing input files containing
// gene expression maps and output files containing expected sample name lists.
// It takes a directory path as input and returns a slice of GetSampleNamesTest structs.
func ReadGetSampleNamesTests(directory string) []GetSampleNamesTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")
	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]GetSampleNamesTest, len(inputFiles))

	// Read input files
	for i, inputFile := range inputFiles {
		f, err := os.Open(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)

		current := GetSampleNamesTest{
			geneExpressionBreastMap: make(map[string]map[string]float64),
		}
		var currentMapKey string

		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}

			switch {
			case strings.HasPrefix(line, "geneName:"):
				currentMapKey = strings.TrimSpace(strings.TrimPrefix(line, "geneName:"))
				current.geneExpressionBreastMap[currentMapKey] = make(map[string]float64)

			default:
				if currentMapKey != "" {
					pair := strings.SplitN(line, ":", 2)
					if len(pair) != 2 {
						continue
					}
					key := strings.TrimSpace(pair[0])
					valStr := strings.TrimSpace(pair[1])
					if valStr == "" {
						continue
					}
					val, err := strconv.ParseFloat(valStr, 64)
					if err != nil {
						panic(err)
					}
					current.geneExpressionBreastMap[currentMapKey][key] = val
				}
			}
		}
		f.Close()
		tests[i] = current
	}

	// Read output files
	for i, outputFile := range outputFiles {
		f, err := os.Open(directory + "/Output/" + outputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)
		outputVals := []string{}
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			for _, s := range strings.Split(line, ",") {
				s = strings.TrimSpace(s)
				if s != "" {
					outputVals = append(outputVals, s)
				}
			}
		}
		f.Close()
		tests[i].result = outputVals
	}

	return tests
}

type GetGeneNamesTest struct {
	geneExpressionBreastMap ExpressionMap
	result                  []string
}

// TestGetGeneNames tests the GetGeneNames function.
func TestGetGeneNames(t *testing.T) {
	tests := ReadGetGeneNamesTests("Tests/GetGeneNames")
	for _, test := range tests {
		got := GetGeneNames(test.geneExpressionBreastMap)
		if !StringSlicesEqual(got, test.result) {
			t.Errorf("GetGeneNames() failed.\nGot: %v\nWant: %v",
				got, test.result)
		}
	}
}

// ReadGetGeneNamesTests reads test data for GetGeneNames function from input files with gene expression data
// and output files with expected sorted gene name lists.
// It takes a directory path as input and returns a slice of GetGeneNamesTest structs.
func ReadGetGeneNamesTests(directory string) []GetGeneNamesTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")
	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]GetGeneNamesTest, len(inputFiles))

	// Read input files
	for i, inputFile := range inputFiles {
		f, err := os.Open(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)

		current := GetGeneNamesTest{
			geneExpressionBreastMap: make(map[string]map[string]float64),
		}
		var currentMapKey string

		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}

			switch {
			case strings.HasPrefix(line, "geneName:"):
				currentMapKey = strings.TrimSpace(strings.TrimPrefix(line, "geneName:"))
				current.geneExpressionBreastMap[currentMapKey] = make(map[string]float64)

			default:
				if currentMapKey != "" {
					pair := strings.SplitN(line, ":", 2)
					if len(pair) != 2 {
						continue
					}
					key := strings.TrimSpace(pair[0])
					valStr := strings.TrimSpace(pair[1])
					if valStr == "" {
						continue
					}
					val, err := strconv.ParseFloat(valStr, 64)
					if err != nil {
						panic(err)
					}
					current.geneExpressionBreastMap[currentMapKey][key] = val
				}
			}
		}
		f.Close()
		tests[i] = current
	}

	// Read output files
	for i, outputFile := range outputFiles {
		f, err := os.Open(directory + "/Output/" + outputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)
		outputVals := []string{}
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			for _, s := range strings.Split(line, ",") {
				s = strings.TrimSpace(s)
				if s != "" {
					outputVals = append(outputVals, s)
				}
			}
		}
		f.Close()
		tests[i].result = outputVals
	}

	return tests
}

type MeanBasedFilterTest struct {
	geneExpressionBreastMap ExpressionMap
	samples                 []string
	meanThresh              float64
	result                  ExpressionMap
}

// TestMeanBasedFilter tests the MeanBasedFilter function.
func TestMeanBasedFilter(t *testing.T) {
	tests := ReadMeanBasedFilterTests("Tests/MeanBasedFilter")
	for _, test := range tests {
		got := MeanBasedFilter(test.geneExpressionBreastMap, test.samples, test.meanThresh)
		if !NestedFloatMapsEqual(got, test.result) {
			t.Errorf("MeanBasedFilter() failed.\nGot: %v\nWant: %v",
				got, test.result)
		}
	}
}

// ReadMeanBasedFilterTests loads test cases for MeanBasedFilter function by parsing input files containing
// gene expression data and thresholds, and output files with filtered results.
// It takes a directory path as input and returns a slice of MeanBasedFilterTest structs.
func ReadMeanBasedFilterTests(directory string) []MeanBasedFilterTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")
	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]MeanBasedFilterTest, len(inputFiles))

	// Parse input files containing gene expression data and filter parameters
	for i, inputFile := range inputFiles {
		f, err := os.Open(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)

		// Initialize test structure for MeanBasedFilter
		current := MeanBasedFilterTest{
			geneExpressionBreastMap: make(ExpressionMap),
		}
		var currentMapKey string

		// Parse structured input with genes, samples, threshold, and expression data
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			// Skip empty lines and comments
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}

			switch {
			// Parse gene name and initialize its expression map
			case strings.HasPrefix(line, "geneName:"):
				currentMapKey = strings.TrimSpace(strings.TrimPrefix(line, "geneName:"))
				current.geneExpressionBreastMap[currentMapKey] = make(map[string]float64)

			// Parse sample names list
			case strings.HasPrefix(line, "samples:"):
				sampleStr := strings.TrimSpace(strings.TrimPrefix(line, "samples:"))
				current.samples = []string{}
				for _, s := range strings.Split(sampleStr, ",") {
					current.samples = append(current.samples, strings.TrimSpace(s))
				}

			// Parse mean threshold for filtering
			case strings.HasPrefix(line, "meanThresh:"):
				threshStr := strings.TrimSpace(strings.TrimPrefix(line, "meanThresh:"))
				thresh, err := strconv.ParseFloat(threshStr, 64)
				if err != nil {
					panic(err)
				}
				current.meanThresh = thresh

			default:
				if currentMapKey != "" {
					pair := strings.SplitN(line, ":", 2)
					if len(pair) != 2 {
						continue
					}
					key := strings.TrimSpace(pair[0])
					valStr := strings.TrimSpace(pair[1])
					if valStr == "" {
						continue
					}
					val, err := strconv.ParseFloat(valStr, 64)
					if err != nil {
						panic(err)
					}
					current.geneExpressionBreastMap[currentMapKey][key] = val
				}
			}
		}
		f.Close()
		tests[i] = current
	}

	// Read output files
	for i, outputFile := range outputFiles {
		f, err := os.Open(directory + "/Output/" + outputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)
		resultMap := make(ExpressionMap)

		var currentGene string
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}

			if strings.HasPrefix(line, "geneName:") {
				currentGene = strings.TrimSpace(strings.TrimPrefix(line, "geneName:"))
				resultMap[currentGene] = make(map[string]float64)
			} else if currentGene != "" {
				pair := strings.SplitN(line, ":", 2)
				if len(pair) != 2 {
					continue
				}
				key := strings.TrimSpace(pair[0])
				valStr := strings.TrimSpace(pair[1])
				if valStr == "" {
					continue
				}
				val, err := strconv.ParseFloat(valStr, 64)
				if err != nil {
					panic(err)
				}
				resultMap[currentGene][key] = val
			}
		}
		f.Close()
		tests[i].result = resultMap
	}

	return tests
}

type SortSampleNamesTest struct {
	samples []string
	result  []string
}

// TestSortSampleNames tests the sortSampleNames function.
func TestSortSampleNames(t *testing.T) {
	tests := ReadSortSampleNamesTests("Tests/SortSampleNames")
	for _, test := range tests {
		// Make a copy so we don't mutate the original test.input
		cpy := make([]string, len(test.samples))
		copy(cpy, test.samples)

		sortSampleNames(cpy)
		if !StringSlicesEqual(cpy, test.result) {
			t.Errorf("sortSampleNames(%v) failed.\nGot: %v\nWant: %v",
				test.samples, cpy, test.result)
		}
	}
}

// ReadSortSampleNamesTests reads test data for the sortSampleNames function by parsing input files with
// sample name lists and output files with expected sorted order.
// It takes a directory path as input and returns a slice of SortSampleNamesTest structs.
func ReadSortSampleNamesTests(directory string) []SortSampleNamesTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")
	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]SortSampleNamesTest, len(inputFiles))

	for i, inputFile := range inputFiles {
		f, err := os.Open(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)
		var input []string
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			input = append(input, line)
		}
		f.Close()
		tests[i].samples = input
	}

	for i, outputFile := range outputFiles {
		f, err := os.Open(directory + "/Output/" + outputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)
		var output []string
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			output = append(output, line)
		}
		f.Close()
		tests[i].result = output
	}

	return tests
}

// Struct for a single test case
type ComputePearsonTest struct {
	sortedGeneNames        []string
	sortedSampleNames      []string
	filteredGeneExpression ExpressionMap
	result                 CorrelationMatrix
}

// TestComputePearsonCorrelation tests the ComputePearsonCorrelation function.
func TestComputePearsonCorrelation(t *testing.T) {
	tests := ReadComputePearsonTests("Tests/ComputePearsonCorrelation")
	for _, test := range tests {
		got := ComputePearsonCorrelation(
			test.sortedGeneNames,
			test.sortedSampleNames,
			test.filteredGeneExpression,
		)
		if !FloatMatrixEqual(got, test.result, 1e-6) {
			t.Errorf("ComputePearsonCorrelation failed.\nGot: %v\nWant: %v",
				got, test.result)
		}
	}
}

// ReadComputePearsonTests loads test cases for ComputePearsonCorrelation function by parsing complex input
// files containing gene expression data and output files with correlation matrices.
// It takes a directory path as input and returns a slice of ComputePearsonTest structs.
func ReadComputePearsonTests(directory string) []ComputePearsonTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")
	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]ComputePearsonTest, len(inputFiles))

	// Parse complex input files containing gene expression data and metadata
	for i, inputFile := range inputFiles {
		f, err := os.Open(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}
		defer f.Close()

		scanner := bufio.NewScanner(f)
		// Initialize test structure with gene expression map
		current := ComputePearsonTest{
			filteredGeneExpression: make(ExpressionMap),
		}
		var currentGene string
		// Parse structured input format with sections for genes, samples, and expression data
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" {
				continue
			}
			switch {
			// Extract sorted gene names from header
			case strings.HasPrefix(line, "geneNames:"):
				current.sortedGeneNames = strings.Fields(strings.TrimPrefix(line, "geneNames:"))
			// Extract sorted sample names from header
			case strings.HasPrefix(line, "sampleNames:"):
				current.sortedSampleNames = strings.Fields(strings.TrimPrefix(line, "sampleNames:"))
			// Start of gene expression block
			case strings.HasSuffix(line, "{"):
				parts := strings.Fields(line)
				currentGene = parts[0]
				current.filteredGeneExpression[currentGene] = make(map[string]float64)
			// End of gene expression block
			case line == "}":
				currentGene = ""
			default:
				if currentGene != "" {
					pair := strings.SplitN(line, ":", 2)
					if len(pair) != 2 {
						continue
					}
					key := strings.TrimSpace(pair[0])
					val, err := strconv.ParseFloat(strings.TrimSpace(pair[1]), 64)
					if err != nil {
						panic(err)
					}
					current.filteredGeneExpression[currentGene][key] = val
				}
			}
		}
		tests[i] = current
	}

	// read output files
	for i, outputFile := range outputFiles {
		f, err := os.Open(directory + "/Output/" + outputFile.Name())
		if err != nil {
			panic(err)
		}
		defer f.Close()

		scanner := bufio.NewScanner(f)
		var outputVals [][]float64
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			parts := strings.Split(line, ",")
			row := make([]float64, 0, len(parts))
			for _, p := range parts {
				val, err := strconv.ParseFloat(strings.TrimSpace(p), 64)
				if err != nil {
					panic(err)
				}
				row = append(row, val)
			}
			outputVals = append(outputVals, row)
		}
		tests[i].result = outputVals
	}

	return tests
}

type TransformMatrixToSliceTest struct {
	matrix CorrelationMatrix
	result []float64
}

// TestTransformMatrixToSlice tests the TransformMatrixToSlice function.
func TestTransformMatrixToSlice(t *testing.T) {
	tests := ReadTransformMatrixToSliceTests("Tests/TransformMatrixToSlice")
	for _, test := range tests {
		got := TransformMatrixToSlice(test.matrix)
		if !FloatSlicesEqualUnordered(got, test.result) {
			t.Errorf("TransformMatrixToSlice failed.\nGot: %v\nWant: %v",
				got, test.result)
		}
	}
}

// ReadTransformMatrixToSliceTests reads test data for TransformMatrixToSlice function by parsing input files
// containing correlation matrices and output files with flattened value arrays.
// It takes a directory path as input and returns a slice of TransformMatrixToSliceTest structs.
func ReadTransformMatrixToSliceTests(directory string) []TransformMatrixToSliceTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")
	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]TransformMatrixToSliceTest, len(inputFiles))

	//read input files containing correlation matrices
	for i, inputFile := range inputFiles {
		f, err := os.Open(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)

		// Initialize test structure for matrix transformation
		current := TransformMatrixToSliceTest{}
		var matrix CorrelationMatrix

		// Parse matrix data line by line
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			// Skip empty lines and comments
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}

			// Expecting matrix rows like: "1.0, 0.2, 0.5"
			rowParts := strings.Split(line, ",")
			row := make([]float64, 0, len(rowParts))

			// Parse each value in the matrix row
			for _, p := range rowParts {
				valStr := strings.TrimSpace(p)
				if valStr == "" {
					continue
				}
				v, err := strconv.ParseFloat(valStr, 64)
				if err != nil {
					panic(err)
				}
				row = append(row, v)
			}

			// Add valid rows to the matrix
			if len(row) > 0 {
				matrix = append(matrix, row)
			}
		}

		f.Close()
		current.matrix = matrix
		tests[i] = current
	}

	//read output files containing expected flattened results
	for i, outputFile := range outputFiles {
		f, err := os.Open(directory + "/Output/" + outputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)

		// Collect expected output values
		outputVals := []float64{}
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			// Skip empty lines and comments
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}

			// Parse comma-separated float values
			for _, p := range strings.Split(line, ",") {
				valStr := strings.TrimSpace(p)
				if valStr == "" {
					continue
				}
				v, err := strconv.ParseFloat(valStr, 64)
				if err != nil {
					panic(err)
				}
				outputVals = append(outputVals, v)
			}
		}

		f.Close()
		tests[i].result = outputVals
	}

	return tests
}

type SortCorrValsTest struct {
	quantileSlice []float64
	result        []float64
}

// TestSortCorrValues tests the SortCorrVals function.
func TestSortCorrValues(t *testing.T) {

	tests := ReadSortCorrValsTests("Tests/SortCorrVals")
	for _, test := range tests {
		got := SortCorrVals(append([]float64{}, test.quantileSlice...))

		if !FloatSlicesEqualUnordered(got, test.result) {
			t.Errorf("SortCorrValues failed .\nquantileSlice: %v\nGot: %v\nWant: %v", test.quantileSlice, got, test.result)
		}
	}

}

// ReadSortCorrValsTests loads test cases for SortCorrVals function by reading input files with unsorted
// correlation values and output files with sorted results.
// It takes a directory path as input and returns a slice of SortCorrValsTest structs.
func ReadSortCorrValsTests(directory string) []SortCorrValsTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")

	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]SortCorrValsTest, len(inputFiles))

	//read input files containing unsorted correlation values
	for i, inputFile := range inputFiles {
		path := directory + "/Input/" + inputFile.Name()
		data, err := os.ReadFile(path)
		if err != nil {
			panic(err)
		}

		// Parse float values from file content
		vals := ParseFloatList(string(data))
		tests[i].quantileSlice = vals
	}

	//read output files containing expected sorted results
	for i, outputFile := range outputFiles {
		path := directory + "/Output/" + outputFile.Name()
		data, err := os.ReadFile(path)
		if err != nil {
			panic(err)
		}

		// Parse expected sorted values
		vals := ParseFloatList(string(data))
		tests[i].result = vals
	}

	return tests
}

// ComputeQuantileTest represents a single test for the ComputeQuntile function.
// It contains the input slice of floats and the expected results (quantiles).
type ComputeQuantileTest struct {
	quantileSlice []float64
	result        []float64
}

// TestComputeQuantile runs tests for the ComputeQuantile function.
func TestComputeQuantile(t *testing.T) {
	//read through all test cases from directory
	tests := ReadComputeQuantileTests("Tests/ComputeQuantile")

	//loop through each test case
	for _, test := range tests {
		//compute quantiles for the input slice
		got := ComputeQuantile(test.quantileSlice)

		// Round each value to 3 decimals for comparison
		for i := range got {
			got[i] = RoundFloat(got[i], 3)
		}

		//Compare computed results with expected results. If they are not the same, an error will print out.
		if !FloatSlicesEqualOrdered(got, test.result) {
			t.Errorf("ComputeQuantile failed.\nInput: %v\nGot: %v\nWant: %v",
				test.quantileSlice, got, test.result)
		}
	}
}

// ReadComputeQuantileTests reads input and output test files from a directory and results a slice of ComputeQuantileTest structs.
// ReadComputeQuantileTests reads test data for ComputeQuantile function by parsing input files containing
// correlation data and output files with computed quantile values.
// It takes a directory path as input and returns a slice of ComputeQuantileTest structs.
func ReadComputeQuantileTests(directory string) []ComputeQuantileTest {
	//list all the input and output files from the directory
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")

	//make sure the number of input and output files match
	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	//initialzie a slice to store all of the test cases
	tests := make([]ComputeQuantileTest, len(inputFiles))

	//Read input files containing correlation data for quantile computation
	for i, inputFile := range inputFiles {
		path := directory + "/Input/" + inputFile.Name()
		data, err := os.ReadFile(path)
		if err != nil {
			panic(err)
		}
		//parse through the input string into a slice of floats
		tests[i].quantileSlice = ParseFloatList(string(data))
	}

	//Read output files containing expected quantile results
	for i, outputFile := range outputFiles {
		path := directory + "/Output/" + outputFile.Name()
		data, err := os.ReadFile(path)
		if err != nil {
			panic(err)
		}
		//parse the expected quantiles into slice of floats
		tests[i].result = ParseFloatList(string(data))
	}

	return tests
}

type BuildGraphTest struct {
	corrMatrix CorrelationMatrix
	geneNames  []string
	threshold  float64
	// expected edges in format: "0:1,3" means node 0 connects to nodes 1 and 3
	expectedEdges map[int][]int
}

// TestBuildGraph runs tests for the BuildGraph function.
func TestBuildGraph(t *testing.T) {
	tests := ReadBuildGraphTests("Tests/BuildGraph")

	for _, test := range tests {
		graph := BuildGraph(test.corrMatrix, test.geneNames, test.threshold)

		// Compare each node’s edges with expected edges
		for nodeID, expected := range test.expectedEdges {
			var got []int
			for _, edge := range graph[nodeID].Edges {
				got = append(got, edge.To.ID)
			}

			if !IntSlicesEqualOrdered(got, expected) {
				t.Errorf("BuildGraph failed.\nNode: %d\nGot: %v\nWant: %v",
					nodeID, got, expected)
			}
		}
	}
}

// ReadBuildGraphTests reads input and output tests from a directory.
// ReadBuildGraphTests loads test cases for BuildGraph function by parsing input files with correlation matrices
// and thresholds, and output files with expected graph structures.
// It takes a directory path as input and returns a slice of BuildGraphTest structs.
func ReadBuildGraphTests(directory string) []BuildGraphTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")

	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]BuildGraphTest, len(inputFiles))

	for i, inputFile := range inputFiles {
		path := directory + "/Input/" + inputFile.Name()
		data, err := os.ReadFile(path)
		if err != nil {
			panic(err)
		}

		tests[i] = parseBuildGraphInput(string(data))
	}

	for i, outputFile := range outputFiles {
		path := directory + "/Output/" + outputFile.Name()
		data, err := os.ReadFile(path)
		if err != nil {
			panic(err)
		}

		tests[i].expectedEdges = parseExpectedEdges(string(data))
	}

	return tests
}

// parseBuildGraphInput parses a single input string containing correlation matrix and threshold data
// into a BuildGraphTest structure for graph construction testing.
// It takes an input data string and returns a BuildGraphTest struct.
func parseBuildGraphInput(data string) BuildGraphTest {
	lines := strings.Split(data, "\n")
	test := BuildGraphTest{}
	stage := "" // Track which section we're currently parsing

	// Parse each line and organize data by section headers
	for _, line := range lines {
		line = strings.TrimSpace(line)
		if line == "" {
			continue
		}

		// Check for section headers that define data structure
		switch line {
		case "threshold:":
			stage = "threshold"
		case "genes:":
			stage = "genes"
		case "matrix:":
			stage = "matrix"
		default:
			// Process data based on current section
			switch stage {
			case "threshold":
				// Parse threshold value for edge creation
				val, _ := strconv.ParseFloat(line, 64)
				test.threshold = val

			case "genes":
				// Parse comma-separated gene names
				test.geneNames = strings.Split(line, ",")

			case "matrix":
				// Parse each matrix row as correlation values
				row := ParseFloatList(line)
				test.corrMatrix = append(test.corrMatrix, row)
			}
		}
	}

	return test
}

// parseExpectedEdges parses expected edge data from string format into a map structure
// representing graph connectivity for test validation.
// It takes a data string and returns a map representing expected edges.
func parseExpectedEdges(data string) map[int][]int {
	result := make(map[int][]int)
	lines := strings.Split(strings.TrimSpace(data), "\n")

	// Parse each line representing node connections
	for _, line := range lines {
		// Split line into node ID and its target connections
		parts := strings.Split(line, ":")
		nodeID, _ := strconv.Atoi(parts[0])

		// Handle nodes with no connections (empty target list)
		if strings.TrimSpace(parts[1]) == "" {
			result[nodeID] = []int{}
			continue
		}

		// Parse comma-separated target node IDs
		targetsStr := strings.Split(parts[1], ",")
		var targets []int
		for _, s := range targetsStr {
			// Convert each target ID string to integer
			id, _ := strconv.Atoi(strings.TrimSpace(s))
			targets = append(targets, id)
		}

		// Store the list of target nodes for this source node
		result[nodeID] = targets
	}

	return result
}

type BuildLouvainGraphTest struct {
	Threshold     float64
	Matrix        [][]float64
	ExpectedEdges [][3]float64 // each edge: {u, v, weight}
}

// TestBuildLouvainGraph tests the BuildLouvainGraph function.
func TestBuildLouvainGraph(t *testing.T) {
	// Read all test cases from input/output directories
	tests := ReadBuildLouvainGraphTests(
		"Tests/BuildLouvainGraph/Input",
		"Tests/BuildLouvainGraph/Output",
	)

	for _, test := range tests {
		// Build the graph from the correlation matrix and threshold
		graph := BuildLouvainGraph(test.Matrix, test.Threshold)

		// Extract edges in a map: key = [2]int{u,v}, value = weight
		gotEdges := ExtractEdgesFromGraph(graph) // map[[2]int]float64

		// Check each expected edge
		for _, e := range test.ExpectedEdges {
			u := int(e[0])
			v := int(e[1])
			wantWeight := e[2]

			// Try both directions for undirected graph
			key1 := [2]int{u, v}
			key2 := [2]int{v, u}

			gotWeight, ok := gotEdges[key1]
			if !ok {
				gotWeight, ok = gotEdges[key2]
			}

			if !ok {
				t.Errorf("missing expected edge %d-%d", u, v)
				continue
			}

			if math.Abs(gotWeight-wantWeight) > 1e-9 {
				t.Errorf("edge weight mismatch for %d-%d: got %.6f, want %.6f",
					u, v, gotWeight, wantWeight)
			}
		}

		// Check for extra edges that were not expected
		expectedEdgeSet := make(map[[2]int]struct{})
		for _, e := range test.ExpectedEdges {
			u := int(e[0])
			v := int(e[1])
			expectedEdgeSet[[2]int{u, v}] = struct{}{}
			expectedEdgeSet[[2]int{v, u}] = struct{}{}
		}

		for key := range gotEdges {
			if _, ok := expectedEdgeSet[key]; !ok {
				t.Errorf("unexpected edge in graph: %v (weight %.6f)", key, gotEdges[key])
			}
		}
	}
}

// ReadBuildLouvainGraphTests loads test cases for BuildLouvainGraph function by reading correlation matrices
// from input files and expected edge structures from output files.
// It takes input and output directory paths and returns a slice of BuildLouvainGraphTest structs.
func ReadBuildLouvainGraphTests(inputDir, outputDir string) []BuildLouvainGraphTest {
	inputFiles := ReadDirectory(inputDir)
	outputFiles := ReadDirectory(outputDir)

	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	var tests []BuildLouvainGraphTest

	// Parse input files containing correlation matrices and thresholds
	for i := range inputFiles {
		test := BuildLouvainGraphTest{}

		// Read and parse input file structure
		inPath := filepath.Join(inputDir, inputFiles[i].Name())
		inData, err := os.ReadFile(inPath)
		if err != nil {
			panic(err)
		}

		scanner := bufio.NewScanner(strings.NewReader(string(inData)))
		section := ""

		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				if strings.HasPrefix(line, "#") {
					section = strings.ToLower(strings.TrimSpace(strings.TrimPrefix(line, "#")))
				}
				continue
			}

			switch section {
			case "threshold":
				val, err := strconv.ParseFloat(line, 64)
				if err != nil {
					panic(err)
				}
				test.Threshold = val
			case "correlation matrix":
				// First line contains matrix dimension
				n, err := strconv.Atoi(line)
				if err != nil {
					panic(err)
				}
				// Initialize matrix and read n rows of n values each
				test.Matrix = make([][]float64, n)
				for r := 0; r < n; r++ {
					if !scanner.Scan() {
						panic("unexpected EOF when reading matrix rows")
					}
					// Parse each row as space-separated float values
					rowStr := strings.Fields(scanner.Text())
					if len(rowStr) != n {
						panic("matrix row does not match expected dimension")
					}
					// Convert string values to float64
					row := make([]float64, n)
					for c, s := range rowStr {
						row[c], err = strconv.ParseFloat(s, 64)
						if err != nil {
							panic(err)
						}
					}
					test.Matrix[r] = row
				}
			}
		}

		//parse output files
		outPath := filepath.Join(outputDir, outputFiles[i].Name())
		outData, err := os.ReadFile(outPath)
		if err != nil {
			panic(err)
		}

		test.ExpectedEdges = ParseEdges(strings.NewReader(string(outData)))

		tests = append(tests, test)
	}

	return tests
}

type RunLouvainOnMatrixTest struct {
	Matrix       CorrelationMatrix // correlation matrix input
	Threshold    float64           // edge inclusion threshold
	ExpectedMod  float64           // expected modularity
	ExpectedComm map[int]int       // expected nodeID -> communityID mapping
}

// TestRunLouvainOnMatrix tests the RunLouvainOnMatrix function.
func TestRunLouvainOnMatrix(t *testing.T) {
	// Load test cases with correlation matrices and expected community assignments
	tests := ReadRunLouvainOnMatrixTests(
		"Tests/RunLouvainOnMatrix/Input",
		"Tests/RunLouvainOnMatrix/Output",
	)

	for _, test := range tests {
		// First build the Louvain graph from correlation matrix
		graph := BuildLouvainGraph(test.Matrix, test.Threshold)
		// Run clustering algorithm to get community assignments and modularity
		gotComm, gotMod := RunLouvainOnMatrix(graph)

		// Ensure modularity is in a reasonable range [-1, 1]
		if gotMod < -1 || gotMod > 1 {
			t.Errorf("modularity out of expected range: %.6f", gotMod)
		}

		// Compare communities: allow relabeling since community IDs can vary
		if !CompareCommunities(gotComm, test.ExpectedComm) {
			t.Errorf("community assignment mismatch.\nGot: %v\nWant: %v", gotComm, test.ExpectedComm)
		}
	}
}

// ReadRunLouvainOnMatrixTests reads test data for RunLouvainOnMatrix function by parsing input files with
// graph data and output files with expected community assignments and modularity values.
// It takes input and output directory paths and returns a slice of RunLouvainOnMatrixTest structs.
func ReadRunLouvainOnMatrixTests(inputDir, outputDir string) []RunLouvainOnMatrixTest {
	inputFiles := ReadDirectory(inputDir)
	outputFiles := ReadDirectory(outputDir)

	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	var tests []RunLouvainOnMatrixTest

	for i := range inputFiles {
		test := RunLouvainOnMatrixTest{}

		//parse input files
		inPath := filepath.Join(inputDir, inputFiles[i].Name())
		inData, err := os.ReadFile(inPath)
		if err != nil {
			panic(err)
		}

		scanner := bufio.NewScanner(strings.NewReader(string(inData)))
		section := ""
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				if strings.HasPrefix(line, "#") {
					section = strings.ToLower(strings.TrimSpace(strings.TrimPrefix(line, "#")))
				}
				continue
			}

			switch section {
			case "threshold":
				val, err := strconv.ParseFloat(line, 64)
				if err != nil {
					panic(err)
				}
				test.Threshold = val
			case "correlation matrix":
				n, err := strconv.Atoi(line)
				if err != nil {
					panic(err)
				}
				test.Matrix = make([][]float64, n)
				for r := 0; r < n; r++ {
					if !scanner.Scan() {
						panic("unexpected EOF reading matrix")
					}
					rowStr := strings.Fields(scanner.Text())
					if len(rowStr) != n {
						panic("matrix row length mismatch")
					}
					row := make([]float64, n)
					for c, s := range rowStr {
						row[c], err = strconv.ParseFloat(s, 64)
						if err != nil {
							panic(err)
						}
					}
					test.Matrix[r] = row
				}
			}
		}

		//parse output files
		outPath := filepath.Join(outputDir, outputFiles[i].Name())
		outData, err := os.ReadFile(outPath)
		if err != nil {
			panic(err)
		}

		test.ExpectedComm = make(map[int]int)
		scanner = bufio.NewScanner(strings.NewReader(string(outData)))
		section = ""
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				if strings.HasPrefix(line, "#") {
					section = strings.ToLower(strings.TrimSpace(strings.TrimPrefix(line, "#")))
				}
				continue
			}

			switch section {
			case "modularity":
				val, err := strconv.ParseFloat(line, 64)
				if err != nil {
					panic(err)
				}
				test.ExpectedMod = val
			case "communities":
				fields := strings.Fields(line)
				if len(fields) != 2 {
					continue
				}
				node, err1 := strconv.Atoi(fields[0])
				cid, err2 := strconv.Atoi(fields[1])
				if err1 != nil || err2 != nil {
					panic("invalid community line: " + line)
				}
				test.ExpectedComm[node] = cid
			}
		}

		tests = append(tests, test)
	}

	return tests
}

// Struct for CountCommunities test case
type CountCommunitiesTest struct {
	CommunityMap map[int]int
	Expected     int
}

// TestCountCommunities tests the CountCommunities function.
func TestCountCommunities(t *testing.T) {
	tests := ReadCountCommunitiesTests("Tests/CountCommunities")

	for _, test := range tests {
		got := CountCommunities(test.CommunityMap)
		if got != test.Expected {
			t.Errorf("CountCommunities failed.\nInput: %v\nGot: %d\nExpected: %d",
				test.CommunityMap, got, test.Expected)
		}
	}
}

// ReadCountCommunitiesTests loads test cases for CountCommunities function by parsing files containing
// community assignments and expected community count values.
// It takes a directory path as input and returns a slice of CountCommunitiesTest structs.
func ReadCountCommunitiesTests(dir string) []CountCommunitiesTest {
	files, err := os.ReadDir(dir)
	if err != nil {
		panic(err)
	}

	var tests []CountCommunitiesTest

	for _, file := range files {
		if file.IsDir() {
			continue
		}
		path := filepath.Join(dir, file.Name())
		data, err := os.ReadFile(path)
		if err != nil {
			panic(err)
		}

		test := CountCommunitiesTest{
			CommunityMap: make(map[int]int),
		}

		scanner := bufio.NewScanner(strings.NewReader(string(data)))
		section := ""
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				if strings.HasPrefix(line, "#") {
					section = strings.ToLower(strings.TrimSpace(strings.TrimPrefix(line, "#")))
				}
				continue
			}

			switch section {
			case "communities":
				// Expected format: nodeID communityID
				fields := strings.Fields(line)
				if len(fields) != 2 {
					continue
				}
				node, err1 := strconv.Atoi(fields[0])
				cid, err2 := strconv.Atoi(fields[1])
				if err1 != nil || err2 != nil {
					panic("invalid community line: " + line)
				}
				test.CommunityMap[node] = cid
			case "expected":
				val, err := strconv.Atoi(line)
				if err != nil {
					panic(err)
				}
				test.Expected = val
			}
		}

		tests = append(tests, test)
	}

	return tests
}

type CalculateNumEdgesTest struct {
	Graph    GraphNetwork
	Expected int
}

// TestCalculateNumEdges tests the CalculateNumEdges function.
func TestCalculateNumEdges(t *testing.T) {
	tests := ReadCalculateNumEdgesTests(
		"Tests/CalculateNumEdges/Input",
		"Tests/CalculateNumEdges/Output",
	)

	for _, test := range tests {
		got := CalculateNumEdges(test.Graph)
		if got != test.Expected {
			t.Errorf("CalculateNumEdges failed.\nGot: %d\nExpected: %d\nGraph: %+v",
				got, test.Expected, test.Graph)
		}
	}
}

// ReadCalculateNumEdgesTests reads test data for CalculateNumEdges function by constructing graph networks
// from input files and reading expected edge counts from output files.
// It takes input and output directory paths and returns a slice of CalculateNumEdgesTest structs.
func ReadCalculateNumEdgesTests(inputDir, outputDir string) []CalculateNumEdgesTest {
	inputFiles, err := os.ReadDir(inputDir)
	if err != nil {
		panic(err)
	}

	outputFiles, err := os.ReadDir(outputDir)
	if err != nil {
		panic(err)
	}

	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]CalculateNumEdgesTest, len(inputFiles))

	// Parse input files to construct graph networks
	for i, inputFile := range inputFiles {
		// Open input file containing graph structure
		path := filepath.Join(inputDir, inputFile.Name())
		file, err := os.Open(path)
		if err != nil {
			panic(err)
		}

		scanner := bufio.NewScanner(file)
		nodeLines := []string{}
		maxNodeID := -1

		// First pass: read all lines and find maximum node ID
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			// Skip empty lines and comments
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			nodeLines = append(nodeLines, line)
			// Extract node ID to determine graph size
			parts := strings.Fields(line)
			nodeID, _ := strconv.Atoi(parts[0])
			if nodeID > maxNodeID {
				maxNodeID = nodeID
			}
		}
		file.Close()

		// Initialize graph with all nodes (create empty nodes first)
		numNodes := maxNodeID + 1
		graph := make(GraphNetwork, numNodes)
		for j := 0; j < numNodes; j++ {
			graph[j] = &Node{
				ID:       j,
				GeneName: "",
				Edges:    []*Edge{},
			}
		}

		// Second pass: parse each line and populate node data and edges
		for _, line := range nodeLines {
			parts := strings.Fields(line)
			nodeID, _ := strconv.Atoi(parts[0])
			// Set gene name for this node
			graph[nodeID].GeneName = parts[1]

			// Parse edges if they exist (format: "toID:weight,toID:weight")
			if len(parts) >= 3 {
				edges := strings.Split(parts[2], ",")
				for _, e := range edges {
					edgeParts := strings.Split(e, ":")
					if len(edgeParts) != 2 {
						continue
					}
					// Create edge with target node and weight
					toID, _ := strconv.Atoi(edgeParts[0])
					weight, _ := strconv.ParseFloat(edgeParts[1], 64)
					edge := &Edge{
						To:     graph[toID], // <-- pointer to actual node
						Weight: weight,
					}
					graph[nodeID].Edges = append(graph[nodeID].Edges, edge)
				}
			}
		}

		tests[i].Graph = graph
	}

	// Read output files
	for i, outputFile := range outputFiles {
		path := filepath.Join(outputDir, outputFile.Name())
		data, err := os.ReadFile(path)
		if err != nil {
			panic(err)
		}
		val, err := strconv.Atoi(strings.TrimSpace(string(data)))
		if err != nil {
			panic(err)
		}
		tests[i].Expected = val
	}

	return tests
}

type ComputeAverageDegreeTest struct {
	NumEdges int
	NumNodes int
	Result   float64
}

// TestComputeAverageDegree runs tests for ComputeAverageDegree
// TestComputeAverageDegree tests the ComputeAverageDegree function.
func TestComputeAverageDegree(t *testing.T) {
	tests := ReadComputeAverageDegreeTests("Tests/ComputeAverageDegree")

	for _, test := range tests {
		got := ComputeAverageDegree(test.NumEdges, test.NumNodes)

		// Round to 6 decimals for comparison
		got = RoundFloat(got, 6)
		want := RoundFloat(test.Result, 6)

		if math.Abs(got-want) > 1e-6 {
			t.Errorf("ComputeAverageDegree failed.\nNumEdges: %d\nNumNodes: %d\nGot: %.6f\nWant: %.6f",
				test.NumEdges, test.NumNodes, got, test.Result)
		}
	}
}

// ReadComputeAverageDegreeTests loads test cases for ComputeAverageDegree function by parsing input files
// with node and edge counts and output files with expected average degree values.
// It takes a directory path as input and returns a slice of ComputeAverageDegreeTest structs.
func ReadComputeAverageDegreeTests(directory string) []ComputeAverageDegreeTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")

	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]ComputeAverageDegreeTest, len(inputFiles))

	for i, inputFile := range inputFiles {
		data, err := os.ReadFile(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}

		lines := strings.Fields(string(data))
		if len(lines) < 2 {
			panic("input file must contain two integers: numEdges numNodes")
		}

		numEdges, err := strconv.Atoi(lines[0])
		if err != nil {
			panic(err)
		}
		numNodes, err := strconv.Atoi(lines[1])
		if err != nil {
			panic(err)
		}

		tests[i].NumEdges = numEdges
		tests[i].NumNodes = numNodes
	}

	for i, outputFile := range outputFiles {
		data, err := os.ReadFile(directory + "/Output/" + outputFile.Name())
		if err != nil {
			panic(err)
		}

		val, err := strconv.ParseFloat(strings.TrimSpace(string(data)), 64)
		if err != nil {
			panic(err)
		}

		tests[i].Result = val
	}

	return tests
}

type ComputeEdgeDensityTest struct {
	NumEdges int
	NumNodes int
	Result   float64
}

// TestComputeEdgeDensity runs tests for ComputeEdgeDensity
// TestComputeEdgeDensity tests the ComputeEdgeDensity function.
func TestComputeEdgeDensity(t *testing.T) {
	tests := ReadComputeEdgeDensityTests("Tests/ComputeEdgeDensity")

	for _, test := range tests {
		got := ComputeEdgeDensity(test.NumEdges, test.NumNodes)

		// Round to 6 decimals for comparison
		got = RoundFloat(got, 6)
		want := RoundFloat(test.Result, 6)

		if math.Abs(got-want) > 1e-6 {
			t.Errorf("ComputeEdgeDensity failed.\nNumEdges: %d\nNumNodes: %d\nGot: %.6f\nWant: %.6f",
				test.NumEdges, test.NumNodes, got, test.Result)
		}
	}
}

// ReadComputeEdgeDensityTests reads test data for ComputeEdgeDensity function by parsing input files
// containing graph parameters and output files with expected density calculations.
// It takes a directory path as input and returns a slice of ComputeEdgeDensityTest structs.
func ReadComputeEdgeDensityTests(directory string) []ComputeEdgeDensityTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")

	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]ComputeEdgeDensityTest, len(inputFiles))

	for i, inputFile := range inputFiles {
		data, err := os.ReadFile(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}

		lines := strings.Fields(string(data))
		if len(lines) < 2 {
			panic("input file must contain two integers: numEdges numNodes")
		}

		numEdges, err := strconv.Atoi(lines[0])
		if err != nil {
			panic(err)
		}
		numNodes, err := strconv.Atoi(lines[1])
		if err != nil {
			panic(err)
		}

		tests[i].NumEdges = numEdges
		tests[i].NumNodes = numNodes
	}

	for i, outputFile := range outputFiles {
		data, err := os.ReadFile(directory + "/Output/" + outputFile.Name())
		if err != nil {
			panic(err)
		}

		val, err := strconv.ParseFloat(strings.TrimSpace(string(data)), 64)
		if err != nil {
			panic(err)
		}

		tests[i].Result = val
	}

	return tests
}

type EdgeStatsTest struct {
	Graph       GraphNetwork
	ExpectedPos int
	ExpectedNeg int
}

// TestEdgeStats tests the EdgeStats function.
func TestEdgeStats(t *testing.T) {
	tests := ReadEdgeStatsTests("Tests/EdgeStats")

	for _, test := range tests {
		gotPos, gotNeg := EdgeStats(test.Graph)

		if gotPos != test.ExpectedPos || gotNeg != test.ExpectedNeg {
			t.Errorf("EdgeStats failed.\nGraph: %+v\nGot pos: %d, neg: %d\nWant pos: %d, neg: %d",
				test.Graph, gotPos, gotNeg, test.ExpectedPos, test.ExpectedNeg)
		}
	}
}

// ReadEdgeStatsTests loads test cases for EdgeStats function by parsing input files with graph network
// data and output files containing positive/negative edge counts.
// It takes a directory path as input and returns a slice of EdgeStatsTest structs.
func ReadEdgeStatsTests(directory string) []EdgeStatsTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")

	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]EdgeStatsTest, len(inputFiles))

	for i, inputFile := range inputFiles {
		path := directory + "/Input/" + inputFile.Name()
		data, err := os.ReadFile(path)
		if err != nil {
			panic(err)
		}
		tests[i].Graph = ParseGraphNetwork(string(data))
	}

	for i, outputFile := range outputFiles {
		path := directory + "/Output/" + outputFile.Name()
		data, err := os.ReadFile(path)
		if err != nil {
			panic(err)
		}
		tests[i].ExpectedPos, tests[i].ExpectedNeg = ParseEdgeStatsOutput(string(data))
	}

	return tests
}

type ComputeModuleSizesTest struct {
	ClusterMap    map[int]int
	ExpectedSizes map[int]int
}

// TestComputeModuleSizes runs all test cases
// TestComputeModuleSizes tests the ComputeModuleSizes function.
func TestComputeModuleSizes(t *testing.T) {
	tests := ReadComputeModuleSizesTests("Tests/ComputeModuleSizes/Input", "Tests/ComputeModuleSizes/Output")

	for _, test := range tests {
		got := ComputeModuleSizes(test.ClusterMap)

		// Compare each expected module size
		for module, expectedSize := range test.ExpectedSizes {
			gotSize, ok := got[module]
			if !ok {
				t.Errorf("missing module %d in result", module)
				continue
			}
			if gotSize != expectedSize {
				t.Errorf("module size mismatch for module %d: got %d, want %d", module, gotSize, expectedSize)
			}
		}

		// Check for extra modules not expected
		for module := range got {
			if _, ok := test.ExpectedSizes[module]; !ok {
				t.Errorf("unexpected module %d in result", module)
			}
		}
	}
}

// ReadComputeModuleSizesTests reads input/output test cases from directories
// ReadComputeModuleSizesTests reads test data for ComputeModuleSizes function by parsing community assignment
// files and expected module size mappings.
// It takes input and output directory paths and returns a slice of ComputeModuleSizesTest structs.
func ReadComputeModuleSizesTests(inputDir, outputDir string) []ComputeModuleSizesTest {
	inputFiles := ReadDirectory(inputDir)
	outputFiles := ReadDirectory(outputDir)

	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]ComputeModuleSizesTest, len(inputFiles))

	for i, inputFile := range inputFiles {
		inputData, err := os.ReadFile(inputDir + "/" + inputFile.Name())
		if err != nil {
			panic(err)
		}
		tests[i].ClusterMap = ParseClusterMap(string(inputData))
	}

	for i, outputFile := range outputFiles {
		outputData, err := os.ReadFile(outputDir + "/" + outputFile.Name())
		if err != nil {
			panic(err)
		}
		tests[i].ExpectedSizes = ParseModuleSizes(string(outputData))
	}

	return tests
}

type TotalNumGenesTest struct {
	BreastGenes  []string
	OvarianGenes []string
	Expected     int
}

// TestTotalNumGenes tests the TotalNumGenes function.
func TestTotalNumGenes(t *testing.T) {
	tests := ReadTotalNumGenesTests("Tests/TotalNumGenes")

	for _, test := range tests {
		got := TotalNumGenes(test.BreastGenes, test.OvarianGenes)
		if got != test.Expected {
			t.Errorf("TotalNumGenes failed.\nBreastGenes: %v\nOvarianGenes: %v\nGot: %d\nWant: %d",
				test.BreastGenes, test.OvarianGenes, got, test.Expected)
		}
	}
}

// ReadTotalNumGenesTests loads test cases for TotalNumGenes function by reading input files with gene name
// lists and output files with expected total unique gene counts.
// It takes a directory path as input and returns a slice of TotalNumGenesTest structs.
func ReadTotalNumGenesTests(directory string) []TotalNumGenesTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")

	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]TotalNumGenesTest, len(inputFiles))

	for i, inputFile := range inputFiles {
		path := directory + "/Input/" + inputFile.Name()
		data, err := os.ReadFile(path)
		if err != nil {
			panic(err)
		}

		// Expecting first line = breast genes, second line = ovarian genes, comma-separated
		lines := strings.Split(strings.TrimSpace(string(data)), "\n")
		if len(lines) != 2 {
			panic("expected 2 lines in input file")
		}

		tests[i].BreastGenes = ParseStringList(lines[0])
		tests[i].OvarianGenes = ParseStringList(lines[1])
	}

	for i, outputFile := range outputFiles {
		path := directory + "/Output/" + outputFile.Name()
		data, err := os.ReadFile(path)
		if err != nil {
			panic(err)
		}
		// single integer
		val, err := strconv.Atoi(strings.TrimSpace(string(data)))
		if err != nil {
			panic(err)
		}
		tests[i].Expected = val
	}

	return tests
}

type LocalClusteringCoeffTest struct {
	GraphData     []string // lines of input representing nodes and edges
	ExpectedCoeff []float64
}

// TestLocalClusteringCoeff tests the LocalClusteringCoeff function.
func TestLocalClusteringCoeff(t *testing.T) {
	tests := ReadLocalClusteringCoeffTests("Tests/LocalClusteringCoeff/Input")

	for i, graph := range tests {
		got := LocalClusteringCoeff(graph)

		t.Logf("Graph %d:", i)
		for j, coeff := range got {
			if math.IsNaN(coeff) {
				t.Logf("Node %d: NaN (less than 2 neighbors)", j)
			} else {
				t.Logf("Node %d: %.4f", j, coeff)
			}
		}
	}
}

// ReadLocalClusteringCoeffTests reads graph network data from input files for testing LocalClusteringCoeff function,
// constructing GraphNetwork structures from text representations.
// It takes an input directory path and returns a slice of GraphNetwork structures.
func ReadLocalClusteringCoeffTests(inputDir string) []GraphNetwork {
	files, _ := os.ReadDir(inputDir)
	var tests []GraphNetwork

	for _, file := range files {
		if file.IsDir() {
			continue
		}
		path := filepath.Join(inputDir, file.Name())
		f, err := os.Open(path)
		if err != nil {
			log.Fatal(err)
		}
		scanner := bufio.NewScanner(f)
		nodes := make(map[int]*Node)
		var nodeOrder []int

		// First pass: create all nodes
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" {
				continue
			}
			parts := strings.Fields(line)
			id, _ := strconv.Atoi(parts[0])
			label := parts[1]
			node := &Node{ID: id, GeneName: label, Edges: []*Edge{}}
			nodes[id] = node
			nodeOrder = append(nodeOrder, id)
		}

		f.Seek(0, 0) // rewind to parse edges
		scanner = bufio.NewScanner(f)
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" {
				continue
			}
			parts := strings.Fields(line)
			id, _ := strconv.Atoi(parts[0])
			if len(parts) < 3 {
				continue
			}
			edgesStr := strings.Split(parts[2], ",")
			for _, e := range edgesStr {
				if e == "" {
					continue
				}
				pair := strings.Split(e, ":")
				neighborID, _ := strconv.Atoi(pair[0])
				weight, _ := strconv.ParseFloat(pair[1], 64)
				nodes[id].Edges = append(nodes[id].Edges, &Edge{
					To:     nodes[neighborID],
					Weight: weight,
				})
			}
		}

		// preserve order
		var graph GraphNetwork
		for _, id := range nodeOrder {
			graph = append(graph, nodes[id])
		}
		tests = append(tests, graph)
		f.Close()
	}
	return tests
}

// BuildGraphFromInputLines constructs a GraphNetwork structure from text lines containing node and edge definitions,
// parsing node IDs, gene names, and edge connections with weights.
// It takes a slice of string lines and returns a GraphNetwork structure.
func BuildGraphFromInputLines(lines []string) GraphNetwork {
	// Initialize graph with space for each node
	graph := make(GraphNetwork, len(lines))

	// Process each line to create nodes and edges
	for _, line := range lines {
		parts := strings.Fields(line)
		id, _ := strconv.Atoi(parts[0]) // Node ID
		name := parts[1]                // Gene name

		// Create node with ID and gene name
		graph[id] = &Node{
			ID:       id,
			GeneName: name,
			Edges:    []*Edge{},
		}

		// Parse edges if they exist (format: "toID:weight,toID:weight")
		if len(parts) > 2 {
			edges := strings.Split(parts[2], ",")
			for _, e := range edges {
				edgeParts := strings.Split(e, ":")
				toID, _ := strconv.Atoi(edgeParts[0])
				weight, _ := strconv.ParseFloat(edgeParts[1], 64)
				// Add edge to current node (note: assumes target node already exists)
				graph[id].Edges = append(graph[id].Edges, &Edge{
					To:     graph[toID],
					Weight: weight,
				})
			}
		}
	}
	return graph
}

type GlobalClusteringCoeffTest struct {
	Input        []float64
	ExpectedMean float64
	ExpectedStd  float64
}

// TestGlobalClusteringCoeff tests the GlobalClusteringCoeff function.
func TestGlobalClusteringCoeff(t *testing.T) {
	inputTests := ReadGlobalClusteringCoeffInputs("Tests/GlobalClusteringCoeff/Input")
	expectedTests := ReadGlobalClusteringCoeffOutputs("Tests/GlobalClusteringCoeff/Output")

	if len(inputTests) != len(expectedTests) {
		t.Fatalf("Mismatch in number of input and output files: %d inputs vs %d outputs",
			len(inputTests), len(expectedTests))
	}

	for i, coeffs := range inputTests {
		gotMean, gotStd := GlobalClusteringCoeff(coeffs)

		// Each output file should contain: mean on line 1, std on line 2
		expectedMean := expectedTests[i][0]
		expectedStd := expectedTests[i][1]

		if math.Abs(gotMean-expectedMean) > 1e-6 {
			t.Errorf("Graph %d: mean mismatch. Got %.6f, want %.6f",
				i, gotMean, expectedMean)
		}

		if math.Abs(gotStd-expectedStd) > 1e-6 {
			t.Errorf("Graph %d: std mismatch. Got %.6f, want %.6f",
				i, gotStd, expectedStd)
		}
	}
}

// ReadGlobalClusteringCoeffInputs reads input files containing local clustering coefficient values for testing
// GlobalClusteringCoeff function, handling NaN values appropriately.
// It takes an input directory path and returns a slice of float64 slices.
func ReadGlobalClusteringCoeffInputs(inputDir string) [][]float64 {
	files, _ := os.ReadDir(inputDir)
	var tests [][]float64

	for _, file := range files {
		if file.IsDir() {
			continue
		}
		path := filepath.Join(inputDir, file.Name())
		f, err := os.Open(path)
		if err != nil {
			log.Fatal(err)
		}
		defer f.Close()

		var coeffs []float64
		scanner := bufio.NewScanner(f)
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" {
				continue
			}
			if strings.ToLower(line) == "nan" {
				coeffs = append(coeffs, math.NaN())
				continue
			}
			val, err := strconv.ParseFloat(line, 64)
			if err != nil {
				log.Fatal(err)
			}
			coeffs = append(coeffs, val)
		}
		tests = append(tests, coeffs)
	}
	return tests
}

// ReadGlobalClusteringCoeffOutputs reads output files containing expected global clustering coefficient values
// for test validation, parsing floating-point results from text files.
// It takes an output directory path and returns a slice of float64 slices.
func ReadGlobalClusteringCoeffOutputs(outputDir string) [][]float64 {
	files, _ := os.ReadDir(outputDir)
	var outputs [][]float64

	for _, file := range files {
		if file.IsDir() {
			continue
		}
		path := filepath.Join(outputDir, file.Name())
		f, err := os.Open(path)
		if err != nil {
			log.Fatal(err)
		}
		defer f.Close()

		var output []float64
		scanner := bufio.NewScanner(f)
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" {
				continue
			}
			val, err := strconv.ParseFloat(line, 64)
			if err != nil {
				log.Fatal(err)
			}
			output = append(output, val)
		}
		outputs = append(outputs, output)
	}
	return outputs
}

type LargestModuleTest struct {
	moduleSizes map[int]int
	resultComm  int
	resultSize  int
}

// TestLargestModule tests the LargestModule function.
func TestLargestModule(t *testing.T) {
	tests := ReadLargestModuleTests("Tests/LargestModules")

	for _, test := range tests {
		gotComm, gotSize := LargestModule(test.moduleSizes)

		if gotComm != test.resultComm || gotSize != test.resultSize {
			t.Errorf("LargestModule(%v) failed.\nGot: (%d, %d)\nWant: (%d, %d)",
				test.moduleSizes, gotComm, gotSize, test.resultComm, test.resultSize)
		}
	}
}

// ReadLargestModuleTests constructs a slice of LargestModuleTest structures from input/output files.
// Parameters are a directory, with a path to the folder containing "Input" and "Output" subfolders with test files.
// It returns a slice of LargestModuleTest containing
//   - moduleSizes: a map[int]int representing node IDs and their module sizes.
//   - resultComm: expected community ID for the largest module.
//   - resultSize: expected size of the largest module.
func ReadLargestModuleTests(directory string) []LargestModuleTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")

	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]LargestModuleTest, len(inputFiles))

	// ---- Read Input Files ----
	for i, inputFile := range inputFiles {
		f, err := os.Open(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)

		current := LargestModuleTest{
			moduleSizes: make(map[int]int),
		}

		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}

			pair := strings.SplitN(line, ":", 2)
			if len(pair) != 2 {
				continue
			}

			keyStr := strings.TrimSpace(pair[0])
			valStr := strings.TrimSpace(pair[1])

			key, err := strconv.Atoi(keyStr)
			if err != nil {
				panic(err)
			}
			val, err := strconv.Atoi(valStr)
			if err != nil {
				panic(err)
			}

			current.moduleSizes[key] = val
		}
		f.Close()
		tests[i] = current
	}

	// ---- Read Output Files ----
	for i, outputFile := range outputFiles {
		f, err := os.Open(directory + "/Output/" + outputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)

		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}

			if strings.HasPrefix(line, "community:") {
				valStr := strings.TrimSpace(strings.TrimPrefix(line, "community:"))
				v, err := strconv.Atoi(valStr)
				if err != nil {
					panic(err)
				}
				tests[i].resultComm = v
			}

			if strings.HasPrefix(line, "size:") {
				valStr := strings.TrimSpace(strings.TrimPrefix(line, "size:"))
				v, err := strconv.Atoi(valStr)
				if err != nil {
					panic(err)
				}
				tests[i].resultSize = v
			}
		}

		f.Close()
	}

	return tests
}

type InvertMapWithGeneNamesTest struct {
	clusterMap map[int]int
	geneNames  []string
	result     map[int][]string
}

// TestInvertMapWithGeneNames tests the InvertMapWithGeneNames function.
func TestInvertMapWithGeneNames(t *testing.T) {
	tests := ReadInvertMapWithGeneNamesTests("Tests/InvertMapWithGeneNames")
	for _, test := range tests {
		got := InvertMapWithGeneNames(test.clusterMap, test.geneNames)
		if !CompareModuleMaps(got, test.result) {
			t.Errorf("InvertMapWithGeneNames(%v, %v) failed.\nGot: %v\nWant: %v",
				test.clusterMap, test.geneNames, got, test.result)
		}
	}
}

// ReadInvertMapWithGeneNamesTests constructs a slice of InvertMapWithGeneNamesTest structures from input/output files.
// Parameters are a directory, with a path to the folder containing "Input" and "Output" subfolders with test files.
// It returns a slice of InvertMapWithGeneNamesTest containing
//   - clusterMap: input map of node IDs to cluster IDs
//   - geneNames: input list of gene names
//   - result: expected inverted map of cluster IDs to gene names
func ReadInvertMapWithGeneNamesTests(directory string) []InvertMapWithGeneNamesTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")
	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]InvertMapWithGeneNamesTest, len(inputFiles))

	for i, inputFile := range inputFiles {
		f, err := os.Open(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)

		current := InvertMapWithGeneNamesTest{
			clusterMap: make(map[int]int),
			geneNames:  []string{},
		}

		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}

			if strings.HasPrefix(line, "geneNames:") {
				namesStr := strings.TrimSpace(strings.TrimPrefix(line, "geneNames:"))
				current.geneNames = strings.Split(namesStr, ",")
				for i := range current.geneNames {
					current.geneNames[i] = strings.TrimSpace(current.geneNames[i])
				}
			} else if strings.HasPrefix(line, "clusterMap:") {
				// subsequent lines are nodeID:communityID
				for scanner.Scan() {
					mapLine := strings.TrimSpace(scanner.Text())
					if mapLine == "" || strings.HasPrefix(mapLine, "#") {
						continue
					}
					pair := strings.SplitN(mapLine, ":", 2)
					if len(pair) != 2 {
						break
					}
					nodeID, err1 := strconv.Atoi(strings.TrimSpace(pair[0]))
					commID, err2 := strconv.Atoi(strings.TrimSpace(pair[1]))
					if err1 != nil || err2 != nil {
						panic("Invalid clusterMap entry")
					}
					current.clusterMap[nodeID] = commID
				}
				break
			}
		}
		f.Close()
		tests[i] = current
	}

	// read output files
	for i, outputFile := range outputFiles {
		f, err := os.Open(directory + "/Output/" + outputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)
		outMap := make(map[int][]string)
		var currentKey int

		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}

			if strings.HasPrefix(line, "communityID:") {
				keyStr := strings.TrimSpace(strings.TrimPrefix(line, "communityID:"))
				key, err := strconv.Atoi(keyStr)
				if err != nil {
					panic(err)
				}
				currentKey = key
				outMap[currentKey] = []string{}
			} else {
				genes := strings.Split(line, ",")
				for _, g := range genes {
					g = strings.TrimSpace(g)
					if g != "" {
						outMap[currentKey] = append(outMap[currentKey], g)
					}
				}
			}
		}
		f.Close()
		tests[i].result = outMap
	}

	return tests
}

type MapIntValuesToFloatSliceTest struct {
	input  map[int]int
	result []float64
}

// TestMapIntValuesToFloatSlice tests the MapIntValuesToFloatSlice function.
func TestMapIntValuesToFloatSlice(t *testing.T) {
	tests := ReadMapIntValuesToFloatSliceTests("Tests/MapIntValuesToFloatSlice")
	for _, test := range tests {
		got := MapIntValuesToFloatSlice(test.input)
		if !FloatSlicesEqualUnordered(got, test.result) {
			t.Errorf("MapIntValuesToFloatSlice(%v) failed.\nGot: %v\nWant: %v",
				test.input, got, test.result)
		}
	}
}

// ReadMapIntValuesToFloatSliceTests constructs a slice of MapIntValuesToFloatSliceTest structures from input/output files.
// Parameters are a directory, with a path to the folder containing "Input" and "Output" subfolders with test files.
// It returns a slice of MapIntValuesToFloatSliceTest containing
//   - input: map of integers representing node IDs and cluster IDs
//   - result: expected slice of float64 after mapping the values
func ReadMapIntValuesToFloatSliceTests(directory string) []MapIntValuesToFloatSliceTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")
	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]MapIntValuesToFloatSliceTest, len(inputFiles))

	for i, inputFile := range inputFiles {
		f, err := os.Open(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)
		current := MapIntValuesToFloatSliceTest{
			input: make(map[int]int),
		}

		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			// expected format: key:value
			parts := strings.SplitN(line, ":", 2)
			if len(parts) != 2 {
				continue
			}
			key, err1 := strconv.Atoi(strings.TrimSpace(parts[0]))
			val, err2 := strconv.Atoi(strings.TrimSpace(parts[1]))
			if err1 != nil || err2 != nil {
				panic("invalid integer in input file")
			}
			current.input[key] = val
		}
		f.Close()
		tests[i] = current
	}

	// read output files
	for i, outputFile := range outputFiles {
		f, err := os.Open(directory + "/Output/" + outputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)
		outputVals := []float64{}
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			for _, p := range strings.Split(line, ",") {
				valStr := strings.TrimSpace(p)
				if valStr == "" {
					continue
				}
				val, err := strconv.ParseFloat(valStr, 64)
				if err != nil {
					panic(err)
				}
				outputVals = append(outputVals, val)
			}
		}
		f.Close()
		tests[i].result = outputVals
	}

	return tests
}

type FilterNaNsTest struct {
	input  []float64
	result []float64
}

// TestFilterNaNs tests the FilterNaNs function.
func TestFilterNaNs(t *testing.T) {
	tests := ReadFilterNaNsTests("Tests/FilterNaNs")
	for _, test := range tests {
		got := FilterNaNs(test.input)
		if !FloatSlicesEqualUnordered(got, test.result) {
			t.Errorf("FilterNaNs(%v) failed.\nGot: %v\nWant: %v",
				test.input, got, test.result)
		}
	}
}

// ReadFilterNaNsTests constructs a slice of FilterNaNsTest structures from input/output files.
// Parameters are a directory, with a path to the folder containing "Input" and "Output" subfolders with test files.
// It returns a slice of FilterNaNsTest containing
//   - input: a slice of float64 values corresponding to clustering coefficient values
//   - result: a slice of float64 values with removed NaN values from the clustering coefficient values
func ReadFilterNaNsTests(directory string) []FilterNaNsTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")

	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]FilterNaNsTest, len(inputFiles))

	// read input files
	for i, inputFile := range inputFiles {
		f, err := os.Open(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}
		defer f.Close()

		scanner := bufio.NewScanner(f)
		inputVals := []float64{}

		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			for _, p := range strings.Split(line, ",") {
				valStr := strings.TrimSpace(p)
				if valStr == "" {
					continue
				}
				if valStr == "NaN" || valStr == "nan" {
					inputVals = append(inputVals, math.NaN())
					continue
				}
				val, err := strconv.ParseFloat(valStr, 64)
				if err != nil {
					panic(err)
				}
				inputVals = append(inputVals, val)
			}
		}
		tests[i].input = inputVals
		f.Close()
	}

	// read output files
	for i, outputFile := range outputFiles {
		f, err := os.Open(directory + "/Output/" + outputFile.Name())
		if err != nil {
			panic(err)
		}
		defer f.Close()

		scanner := bufio.NewScanner(f)
		outputVals := []float64{}

		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			for _, p := range strings.Split(line, ",") {
				valStr := strings.TrimSpace(p)
				if valStr == "" {
					continue
				}
				val, err := strconv.ParseFloat(valStr, 64)
				if err != nil {
					panic(err)
				}
				outputVals = append(outputVals, val)
			}
		}
		tests[i].result = outputVals
		f.Close()
	}

	return tests
}

type ComputeDegreesTest struct {
	graph  GraphNetwork
	result []float64
}

// TestComputeDegrees tests the ComputeDegrees function.
func TestComputeDegrees(t *testing.T) {
	tests := ReadComputeDegreesTests("Tests/ComputeDegrees")
	for _, test := range tests {
		got := ComputeDegrees(test.graph)
		if !FloatSlicesEqualUnordered(got, test.result) {
			t.Errorf("ComputeDegrees() failed.\nGot: %v\nWant: %v", got, test.result)
		}
	}
}

// ReadComputeDegreesTests constructs a slice of ComputeDegreesTest structures from input/output files.
// Parameters are a directory, with a path to the folder containing "Input" and "Output" subfolders with test files.
// It returns a slice of ComputeDegreesTest containing
//   - graph: an input GraphNetwork that represents a slice of Node pointers
//   - result: a slice of float64 values corresponding to the expected degrees of each node
func ReadComputeDegreesTests(directory string) []ComputeDegreesTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")
	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]ComputeDegreesTest, len(inputFiles))

	for i, inputFile := range inputFiles {
		f, err := os.Open(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)

		nodes := GraphNetwork{}
		var nodeMap = make(map[int]*Node)

		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			// line format: NodeID,GeneName,Neighbor1ID,Neighbor2ID,...
			parts := strings.Split(line, ",")
			nodeID, _ := strconv.Atoi(parts[0])
			node := &Node{ID: nodeID, GeneName: parts[1]}
			nodes = append(nodes, node)
			nodeMap[nodeID] = node
		}
		f.Close()

		// second pass to add edges
		f, _ = os.Open(directory + "/Input/" + inputFile.Name())
		scanner = bufio.NewScanner(f)
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			parts := strings.Split(line, ",")
			nodeID, _ := strconv.Atoi(parts[0])
			node := nodeMap[nodeID]
			for _, neighborStr := range parts[2:] {
				if neighborStr == "" {
					continue
				}
				neighborID, _ := strconv.Atoi(neighborStr)
				neighborNode, exists := nodeMap[neighborID]
				if exists {
					node.Edges = append(node.Edges, &Edge{To: neighborNode})
				}
			}
		}
		f.Close()

		tests[i].graph = nodes
	}

	// read output files
	for i, outputFile := range outputFiles {
		f, err := os.Open(directory + "/Output/" + outputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)
		outputVals := []float64{}
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			for _, p := range strings.Split(line, ",") {
				valStr := strings.TrimSpace(p)
				if valStr == "" {
					continue
				}
				val, err := strconv.ParseFloat(valStr, 64)
				if err != nil {
					panic(err)
				}
				outputVals = append(outputVals, val)
			}
		}
		f.Close()
		tests[i].result = outputVals
	}

	return tests
}

type ECDFTest struct {
	Sample []float64
	X      float64
	Result float64
}

// TestCalculateECDF tests the calculateECDF function.
func TestCalculateECDF(t *testing.T) {
	tests, err := ReadECDFTests("calculateECDF")
	if err != nil {
		t.Fatalf("Failed to read ECDF tests: %v", err)
	}

	for i, test := range tests {
		got := calculateECDF(test.Sample, test.X)
		if math.Abs(got-test.Result) > 1e-6 {
			t.Errorf("Test %d failed: got %f, want %f", i+1, got, test.Result)
		}
	}
}

// ReadECDFTests constructs a slice of ECDFTest structures from input/output files.
// Parameters are a directory, with a path to the folder containing "Input" and "Output" subfolders with test files.
// It returns a slice of ECDFTest containing
//   - Sample: slice of int input values regarding to the sample
//   - X: float64 value which is the point at which the ECDF is evaluated
//   - result: expected ECDF value of X
func ReadECDFTests(folder string) ([]ECDFTest, error) {
	inputFiles, err := filepath.Glob(filepath.Join(folder, "input", "*.txt"))
	if err != nil {
		return nil, err
	}

	var tests []ECDFTest
	for _, inputFile := range inputFiles {
		base := filepath.Base(inputFile)
		outFile := filepath.Join(folder, "output", base)

		// Read input
		data, err := os.ReadFile(inputFile)
		if err != nil {
			return nil, err
		}
		lines := strings.Split(strings.TrimSpace(string(data)), "\n")
		if len(lines) < 2 {
			return nil, fmt.Errorf("input file %s must have 2 lines", inputFile)
		}

		// Parse sample
		strVals := strings.Split(lines[0], ",")
		sample := make([]float64, len(strVals))
		for i, s := range strVals {
			sample[i], err = strconv.ParseFloat(strings.TrimSpace(s), 64)
			if err != nil {
				return nil, fmt.Errorf("invalid float in %s: %v", inputFile, err)
			}
		}

		// Parse x
		x, err := strconv.ParseFloat(strings.TrimSpace(lines[1]), 64)
		if err != nil {
			return nil, fmt.Errorf("invalid x value in %s: %v", inputFile, err)
		}

		// Read expected output
		outData, err := os.ReadFile(outFile)
		if err != nil {
			return nil, err
		}
		expVal, err := strconv.ParseFloat(strings.TrimSpace(string(outData)), 64)
		if err != nil {
			return nil, fmt.Errorf("invalid output in %s: %v", outFile, err)
		}

		tests = append(tests, ECDFTest{
			Sample: sample,
			X:      x,
			Result: expVal,
		})
	}
	return tests, nil
}

type KSTwoSampleTest struct {
	Sample1 []float64
	Sample2 []float64
	Result  float64
}

// TestKSTwoSampleStatistic tests the KSTwoSampleStatistic function.
func TestKSTwoSampleStatistic(t *testing.T) {
	tests, err := ReadKSTwoSampleStatisticTest("Tests/KSTwoSampleStatistic")
	if err != nil {
		t.Fatalf("Failed to read KS tests: %v", err)
	}

	for i, test := range tests {
		got := KSTwoSampleStatistic(test.Sample1, test.Sample2)
		if math.Abs(got-test.Result) > 1e-6 {
			t.Errorf("Test %d failed: got %f, want %f", i+1, got, test.Result)
		}
	}
}

// ReadKSTwoSampleStatisticTest constructs a slice of KSTwoSampleTest structures from input/output files.
// Parameters are a directory, with a path to the folder containing "Input" and "Output" subfolders with test files.
// It returns a slice of KSTwoSampleTest containing
//   - Sample1: slice of int input values regarding to dataset1
//   - Sample2: slice of int input values regarding to dataset2
//   - Result: expected KS Test statistic D in float64 format
func ReadKSTwoSampleStatisticTest(folder string) ([]KSTwoSampleTest, error) {
	// Find all input files
	inputFiles, err := filepath.Glob(filepath.Join(folder, "input", "*.txt"))
	if err != nil {
		return nil, err
	}

	var tests []KSTwoSampleTest
	for i, inputFile := range inputFiles {
		// Construct the corresponding output file name (output1.txt, output2.txt, ...)
		outFile := filepath.Join(folder, "output", fmt.Sprintf("output%d.txt", i+1))

		// --- Read input file ---
		data, err := os.ReadFile(inputFile)
		if err != nil {
			return nil, err
		}
		lines := strings.Split(strings.TrimSpace(string(data)), "\n")
		if len(lines) < 2 {
			return nil, fmt.Errorf("input file %s must have 2 lines", inputFile)
		}

		// Parse a line into []float64
		parseSample := func(line string) ([]float64, error) {
			parts := strings.Split(line, ",")
			sample := make([]float64, len(parts))
			for i, s := range parts {
				val, err := strconv.ParseFloat(strings.TrimSpace(s), 64)
				if err != nil {
					return nil, fmt.Errorf("invalid float in %s: %v", inputFile, err)
				}
				sample[i] = val
			}
			return sample, nil
		}

		sample1, err := parseSample(lines[0])
		if err != nil {
			return nil, err
		}
		sample2, err := parseSample(lines[1])
		if err != nil {
			return nil, err
		}

		// --- Read expected output file ---
		outData, err := os.ReadFile(outFile)
		if err != nil {
			return nil, err
		}
		expected, err := strconv.ParseFloat(strings.TrimSpace(string(outData)), 64)
		if err != nil {
			return nil, fmt.Errorf("invalid float in output file %s: %v", outFile, err)
		}

		// Append test case
		tests = append(tests, KSTwoSampleTest{
			Sample1: sample1,
			Sample2: sample2,
			Result:  expected,
		})
	}

	return tests, nil
}

type KSTestCase struct {
	Dataset1 []float64
	Dataset2 []float64
	DStat    float64
	PValue   float64
}

// TestKSTwoSampleStatistic tests the KSTwoSampleStatistic function.
func TestKSTest(t *testing.T) {
	tests, err := ReadKSTestCases("Tests/KSTest")
	if err != nil {
		t.Fatalf("Failed to read KS tests: %v", err)
	}

	for i, test := range tests {
		gotD, gotP := KSTest(test.Dataset1, test.Dataset2)
		if math.Abs(gotD-test.DStat) > 1e-6 || math.Abs(gotP-test.PValue) > 1e-6 {
			t.Errorf("Test %d failed:\nGot: D=%.6f, P=%.6f\nWant: D=%.6f, P=%.6f",
				i+1, gotD, gotP, test.DStat, test.PValue)
		}
	}
}

// ReadKSTestCases constructs a slice of KSTestCase structures from input/output files.
// Parameters are a directory, with a path to the folder containing "Input" and "Output" subfolders with test files.
// It returns a slice of KSTestCase containing
//   - Dataset1: slice of float64 values for the first empirical sample
//   - Dataset2: slice of float64 values for the second empirical sample
//   - DStat: expected Kolmogorov-Smirnov D Statistic
//   - PValue: expected p=value corresponding to the KS Test
func ReadKSTestCases(folder string) ([]KSTestCase, error) {
	inputFiles, err := filepath.Glob(filepath.Join(folder, "input", "*.txt"))
	if err != nil {
		return nil, err
	}

	var tests []KSTestCase
	for i, inputFile := range inputFiles {
		// corresponding output file with same base name
		outFile := filepath.Join(folder, "output", fmt.Sprintf("output%d.txt", i+1))

		// Read input
		data, err := os.ReadFile(inputFile)
		if err != nil {
			return nil, err
		}
		lines := strings.Split(strings.TrimSpace(string(data)), "\n")
		if len(lines) < 2 {
			return nil, fmt.Errorf("input file %s must have at least 2 lines", inputFile)
		}

		parseSample := func(line string) ([]float64, error) {
			parts := strings.Split(line, ",")
			sample := make([]float64, len(parts))
			for i, s := range parts {
				val, err := strconv.ParseFloat(strings.TrimSpace(s), 64)
				if err != nil {
					return nil, fmt.Errorf("invalid float in %s: %v", inputFile, err)
				}
				sample[i] = val
			}
			return sample, nil
		}

		dataset1, err := parseSample(lines[0])
		if err != nil {
			return nil, err
		}
		dataset2, err := parseSample(lines[1])
		if err != nil {
			return nil, err
		}

		// Read expected output (two lines: first D statistic, second p-value)
		outData, err := os.ReadFile(outFile)
		if err != nil {
			return nil, err
		}
		outLines := strings.Split(strings.TrimSpace(string(outData)), "\n")
		if len(outLines) < 2 {
			return nil, fmt.Errorf("output file %s must have 2 lines", outFile)
		}
		dStat, err := strconv.ParseFloat(strings.TrimSpace(outLines[0]), 64)
		if err != nil {
			return nil, fmt.Errorf("invalid D statistic in %s: %v", outFile, err)
		}
		pVal, err := strconv.ParseFloat(strings.TrimSpace(outLines[1]), 64)
		if err != nil {
			return nil, fmt.Errorf("invalid p-value in %s: %v", outFile, err)
		}

		tests = append(tests, KSTestCase{
			Dataset1: dataset1,
			Dataset2: dataset2,
			DStat:    dStat,
			PValue:   pVal,
		})
	}

	return tests, nil
}
