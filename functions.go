package main

// Function TransformMapToSlice takes a map's outerkey (geneName) and the map its in as input
// and returns a slice corresponding to the float64 tpm values for each of the samples for that gene
func TransformMapToSlice(geneName string, geneExpressionBreastMap map[string]map[string]float64) []float64 {
	//initialize an empty slice to for a geneName, the float64 values are the tpm values for each of the samples
	geneNameslice := make([]float64, 0)

	//lookup that input geneName key in the map. The nested structure returns an inner map corresponding to
	//the inner keys being the sample names and values are the tpm values
	innerMap := geneExpressionBreastMap[geneName]

	//range through the inner map, and create the slice of tpm values
	for _, tpmVal := range innerMap {
		geneNameslice = append(geneNameslice, tpmVal)
	}

	//return the slice (pseudo vector) for one gene
	return geneNameslice

}

// Store all of the sample names from the Data
// Function GetSampleNames takes a map[string]map[string]float64 as input and grabs the inner string key corresponding to the
// sample names and returns it as a slice of sample names
func GetSampleNames(geneExpressionBreastMap map[string]map[string]float64) []string {

	samples := make([]string, 0)

	for _, innerMap := range geneExpressionBreastMap {
		for sample := range innerMap {
			// check if it's already in the list
			found := false
			for _, s := range samples {
				if s == sample {
					found = true
					break
				}
			}
			if !found {
				samples = append(samples, sample)
			}
		}
	}
	return samples //return slice of string names
}

// Use Mean-Based Filtering in order to normalize data while maintaining fixed number of samples across each gene
// Function MeanBasedFilter takes map[string]map[string]float64 as input and returns a filtered map[string]map[string]float64 as output
// Step 1: Find/confirm number of samples for each gene (21 samples)
// Step 2: For each gene, Calculate the mean of the tpm values across its samples and then,
// Step 3: If the mean of that gene's tpm values is <=10 then remove the gene from further analyses
func MeanBasedFilter(geneExpressionBreastMap map[string]map[string]float64, samples []string, meanThresh float64) map[string]map[string]float64 {

	//initialize a nested map for return corresponding to the input map filtered for only genes with a mean tpm greater than 10.0
	out := make(map[string]map[string]float64, len(geneExpressionBreastMap))

	//find the number of samples
	n := len(samples)

	//range through the input nested map
	for gene, row := range geneExpressionBreastMap {
		// Sum the tpm values over all samples for a given gene
		var sum float64 //declare the sum within the loop so that it resets to 0 with each loop for every gene
		for _, s := range samples {
			sum += row[s] //sum the tpm values for all the samples for a given gene
		}
		//calculate the mean by dividing the tpm value sum by number of samples
		mean := sum / float64(n)

		// Keep gene only if mean TPM > threshold(meanThresh)
		if mean > meanThresh {
			cpy := make(map[string]float64, len(row))
			for s, v := range row {
				cpy[s] = v
			}
			out[gene] = cpy
		}
	}
	return out
}
