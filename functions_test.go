package main

import (
	"bufio"
	"os"
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
		if !floatSlicesEqualUnordered(got, test.result) {
			t.Errorf("TransformMapToSlice(%s) failed.\nGot: %v\nWant: %v",
				test.geneName, got, test.result)
		}
	}
}

// ReadTransformMapToSliceTests reads Input/Output files for TransformMapToSlice tests
func ReadTransformMapToSliceTests(directory string) []TransformMapToSliceTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")
	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]TransformMapToSliceTest, len(inputFiles))

	for i, inputFile := range inputFiles {
		f, err := os.Open(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)

		current := TransformMapToSliceTest{
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

type GetSampleNamesTest struct {
	geneExpressionBreastMap ExpressionMap
	result                  []string
}

// TestGetSampleNames tests the GetSampleNames function.
func TestGetSampleNames(t *testing.T) {
	tests := ReadGetSampleNamesTests("Tests/GetSampleNames")
	for _, test := range tests {
		got := GetSampleNames(test.geneExpressionBreastMap)
		if !stringSlicesEqualUnordered(got, test.result) {
			t.Errorf("GetSampleNames() failed.\nGot: %v\nWant: %v",
				got, test.result)
		}
	}
}

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

func TestGetGeneNames(t *testing.T) {
	tests := ReadGetGeneNamesTests("Tests/GetGeneNames")
	for _, test := range tests {
		got := GetGeneNames(test.geneExpressionBreastMap)
		if !stringSlicesEqual(got, test.result) {
			t.Errorf("GetGeneNames() failed.\nGot: %v\nWant: %v",
				got, test.result)
		}
	}
}

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

func TestMeanBasedFilter(t *testing.T) {
	tests := ReadMeanBasedFilterTests("Tests/MeanBasedFilter")
	for _, test := range tests {
		got := MeanBasedFilter(test.geneExpressionBreastMap, test.samples, test.meanThresh)
		if !nestedFloatMapsEqual(got, test.result) {
			t.Errorf("MeanBasedFilter() failed.\nGot: %v\nWant: %v",
				got, test.result)
		}
	}
}

func ReadMeanBasedFilterTests(directory string) []MeanBasedFilterTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")
	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]MeanBasedFilterTest, len(inputFiles))

	for i, inputFile := range inputFiles {
		f, err := os.Open(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)

		current := MeanBasedFilterTest{
			geneExpressionBreastMap: make(ExpressionMap),
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

			case strings.HasPrefix(line, "samples:"):
				sampleStr := strings.TrimSpace(strings.TrimPrefix(line, "samples:"))
				current.samples = []string{}
				for _, s := range strings.Split(sampleStr, ",") {
					current.samples = append(current.samples, strings.TrimSpace(s))
				}

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

func TestSortSampleNames(t *testing.T) {
	tests := ReadSortSampleNamesTests("Tests/SortSampleNames")
	for _, test := range tests {
		// Make a copy so we don't mutate the original test.input
		cpy := make([]string, len(test.samples))
		copy(cpy, test.samples)

		sortSampleNames(cpy)
		if !stringSlicesEqual(cpy, test.result) {
			t.Errorf("sortSampleNames(%v) failed.\nGot: %v\nWant: %v",
				test.samples, cpy, test.result)
		}
	}
}

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

func TestComputePearsonCorrelation(t *testing.T) {
	tests := ReadComputePearsonTests("Tests/ComputePearsonCorrelation")
	for _, test := range tests {
		got := ComputePearsonCorrelation(
			test.sortedGeneNames,
			test.sortedSampleNames,
			test.filteredGeneExpression,
		)
		if !floatMatrixEqual(got, test.result, 1e-6) {
			t.Errorf("ComputePearsonCorrelation failed.\nGot: %v\nWant: %v",
				got, test.result)
		}
	}
}

// ReadComputePearsonTests parses the input/output test files
func ReadComputePearsonTests(directory string) []ComputePearsonTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")
	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]ComputePearsonTest, len(inputFiles))

	for i, inputFile := range inputFiles {
		f, err := os.Open(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}
		defer f.Close()

		scanner := bufio.NewScanner(f)
		current := ComputePearsonTest{
			filteredGeneExpression: make(ExpressionMap),
		}
		var currentGene string
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" {
				continue
			}
			switch {
			case strings.HasPrefix(line, "geneNames:"):
				current.sortedGeneNames = strings.Fields(strings.TrimPrefix(line, "geneNames:"))
			case strings.HasPrefix(line, "sampleNames:"):
				current.sortedSampleNames = strings.Fields(strings.TrimPrefix(line, "sampleNames:"))
			case strings.HasSuffix(line, "{"):
				parts := strings.Fields(line)
				currentGene = parts[0]
				current.filteredGeneExpression[currentGene] = make(map[string]float64)
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
		if !floatSlicesEqualUnordered(got, test.result) {
			t.Errorf("TransformMatrixToSlice failed.\nGot: %v\nWant: %v",
				got, test.result)
		}
	}
}

// ReadTransformMatrixToSliceTests reads Input/Output files for matrix tests
func ReadTransformMatrixToSliceTests(directory string) []TransformMatrixToSliceTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")
	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]TransformMatrixToSliceTest, len(inputFiles))

	//read input files
	for i, inputFile := range inputFiles {
		f, err := os.Open(directory + "/Input/" + inputFile.Name())
		if err != nil {
			panic(err)
		}
		scanner := bufio.NewScanner(f)

		current := TransformMatrixToSliceTest{}
		var matrix CorrelationMatrix

		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}

			// Expecting matrix rows like:
			// 1.0, 0.2, 0.5
			rowParts := strings.Split(line, ",")
			row := make([]float64, 0, len(rowParts))

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

			if len(row) > 0 {
				matrix = append(matrix, row)
			}
		}

		f.Close()
		current.matrix = matrix
		tests[i] = current
	}

	//read output files
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

func TestSortCorrValues(t *testing.T) {

	tests := ReadSortCorrValsTests("Tests/SortCorrVals")
	for _, test := range tests {
		got := SortCorrVals(append([]float64{}, test.quantileSlice...))

		if !floatSlicesEqualUnordered(got, test.result) {
			t.Errorf("SortCorrValues failed .\nquantileSlice: %v\nGot: %v\nWant: %v", test.quantileSlice, got, test.result)
		}
	}

}

func ReadSortCorrValsTests(directory string) []SortCorrValsTest {
	inputFiles := ReadDirectory(directory + "/Input")
	outputFiles := ReadDirectory(directory + "/Output")

	if len(inputFiles) != len(outputFiles) {
		panic("number of input and output files do not match")
	}

	tests := make([]SortCorrValsTest, len(inputFiles))

	//read input files
	for i, inputFile := range inputFiles {
		path := directory + "/Input/" + inputFile.Name()
		data, err := os.ReadFile(path)
		if err != nil {
			panic(err)
		}

		vals := parseFloatList(string(data))
		tests[i].quantileSlice = vals
	}

	//read output files
	for i, outputFile := range outputFiles {
		path := directory + "/Output/" + outputFile.Name()
		data, err := os.ReadFile(path)
		if err != nil {
			panic(err)
		}

		vals := parseFloatList(string(data))
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
			got[i] = roundFloat(got[i], 3)
		}

		//Compare computed results with expected results. If they are not the same, an error will print out.
		if !floatSlicesEqualOrdered(got, test.result) {
			t.Errorf("ComputeQuantile failed.\nInput: %v\nGot: %v\nWant: %v",
				test.quantileSlice, got, test.result)
		}
	}
}

// ReadComputeQuantileTests reads input and output test files from a directory and results a slice of ComputeQuantileTest structs.
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

	//Read input files
	for i, inputFile := range inputFiles {
		path := directory + "/Input/" + inputFile.Name()
		data, err := os.ReadFile(path)
		if err != nil {
			panic(err)
		}
		//parse through the input string into a slice of floats
		tests[i].quantileSlice = parseFloatList(string(data))
	}

	//Read output files
	for i, outputFile := range outputFiles {
		path := directory + "/Output/" + outputFile.Name()
		data, err := os.ReadFile(path)
		if err != nil {
			panic(err)
		}
		//parse the expected quantiles into slice of floats
		tests[i].result = parseFloatList(string(data))
	}

	return tests
}
