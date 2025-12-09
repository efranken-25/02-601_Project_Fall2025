package main

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

/*

Step 2: filter low-expression genes

1. Read in the data from a file. Parse the file
2. Only keep data columns and map string to float64: gene_name, tpm_unstranded
3. Filter gene_type and keep only those labeled as "protein coding" in column gene_type
//make map of maps preserving sample column inner key as BRCA#
*/

// ParseGeneExpressionFile reads a tab-separated gene expression file and returns
// a map from gene_name -> tpm_unstranded for rows where gene_type == "protein_coding"
// *Copilot Generated*
func ParseGeneExpressionFile(path string) (map[string]float64, error) {
	geneTpm := make(map[string]float64)

	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)

	headerFound := false
	idxGeneName, idxGeneType, idxTpm := -1, -1, -1

	for scanner.Scan() {
		line := scanner.Text()
		if strings.TrimSpace(line) == "" {
			continue
		}
		if strings.HasPrefix(line, "#") {
			continue
		}

		fields := strings.Split(line, "\t")

		if !headerFound {
			for i, h := range fields {
				switch strings.TrimSpace(h) {
				case "gene_name":
					idxGeneName = i
				case "gene_type":
					idxGeneType = i
				case "tpm_unstranded":
					idxTpm = i
				}
			}
			if idxGeneName == -1 || idxGeneType == -1 || idxTpm == -1 {
				return nil, fmt.Errorf("required columns (gene_name, gene_type, tpm_unstranded) not found in header")
			}
			headerFound = true
			continue
		}

		if idxGeneType >= len(fields) || idxGeneName >= len(fields) || idxTpm >= len(fields) {
			continue
		}

		if strings.TrimSpace(fields[idxGeneType]) != "protein_coding" {
			continue
		}

		name := strings.TrimSpace(fields[idxGeneName])
		if name == "" {
			continue
		}

		tpmStr := strings.TrimSpace(fields[idxTpm])
		if tpmStr == "" {
			continue
		}

		tpm, err := strconv.ParseFloat(tpmStr, 64)
		if err != nil {
			continue
		}

		geneTpm[name] = tpm
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return geneTpm, nil
}

// ReadGeneExpressionDirToGeneMap reads all files in dir that match "*_gene_counts.tsv"
// and returns a map[gene_name]map[sampleID]tpm_unstranded. It uses ParseGeneExpressionFile
// (which already filters to protein_coding).
// *Copilot Generated*
func ReadGeneExpressionDirToGeneMap(dir string) (map[string]map[string]float64, error) {
	geneMap := make(map[string]map[string]float64)

	entries, err := os.ReadDir(dir)
	if err != nil {
		return nil, err
	}

	for _, e := range entries {
		if e.IsDir() {
			continue
		}
		name := e.Name()
		if strings.HasPrefix(name, ".") {
			continue
		}
		if !strings.HasSuffix(name, "_gene_counts.tsv") {
			continue
		}

		full := filepath.Join(dir, name)
		sample := strings.TrimSuffix(name, "_gene_counts.tsv")
		// fallback: if trimming didn't change name (unexpected), strip extension
		if sample == name {
			sample = strings.TrimSuffix(name, filepath.Ext(name))
		}

		m, err := ParseGeneExpressionFile(full)
		if err != nil {
			return nil, fmt.Errorf("parsing %s: %w", full, err)
		}

		for gene, tpm := range m {
			if _, ok := geneMap[gene]; !ok {
				geneMap[gene] = make(map[string]float64)
			}
			geneMap[gene][sample] = tpm
		}
	}

	return geneMap, nil
}

// *ChatGPT generated*
// WriteEdgesCSV writes the edges of a graph network to a CSV file.
// The CSV contains columns: from, to, weight.
// To avoid duplicate edges in undirected graphs, only writes edges where source ID < target ID.
func WriteEdgesCSV(filename string, breastGraph GraphNetwork) error {
	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	w := csv.NewWriter(f)
	defer w.Flush()

	// header
	if err := w.Write([]string{"from", "to", "weight"}); err != nil {
		return err
	}

	// Avoid duplicate edges: since graph is undirected and you stored edges both ways,
	// only write when u.ID < v.ID.
	for _, u := range breastGraph {
		for _, e := range u.Edges {
			v := e.To
			if u.ID < v.ID {
				record := []string{
					u.GeneName, // from (gene name)
					v.GeneName, // to   (gene name)
					strconv.FormatFloat(e.Weight, 'f', 6, 64),
				}
				if err := w.Write(record); err != nil {
					return err
				}
			}
		}
	}

	return nil
}

// *ChatGPT generated*
// WriteCommunitiesCSV writes node community assignments to a CSV file.
// The CSV contains columns: id, label, community.
// Each row represents a gene with its assigned community ID.
func WriteCommunitiesCSV(filename string, breastGeneNames []string, communities map[int]int) error {
	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	w := csv.NewWriter(f)
	defer w.Flush()

	// header
	if err := w.Write([]string{"id", "label", "community"}); err != nil {
		return err
	}

	for nodeID, gene := range breastGeneNames {
		commID, ok := communities[nodeID]
		if !ok {
			// should not usually happen, but skip if not present
			continue
		}

		record := []string{
			gene,                 // id (used by visNetwork internally)
			gene,                 // label (what you see)
			strconv.Itoa(commID), // community (cluster ID)
		}

		if err := w.Write(record); err != nil {
			return err
		}
	}

	return nil
}

// *ChatGPT generated*
// WriteCommunityStatsCSV writes per-community statistics to a CSV file.
// Columns: community_id, num_nodes, num_edges, density
func WriteCommunityStatsCSV(
	filename string,
	moduleSizes map[int]int,
	moduleEdges map[int]int,
	moduleDensities map[int]float64,
) error {

	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	w := csv.NewWriter(f)
	defer w.Flush()

	// header
	if err := w.Write([]string{"community_id", "num_nodes", "num_edges", "density"}); err != nil {
		return err
	}

	// sort community IDs for deterministic output
	keys := make([]int, 0, len(moduleSizes))
	for k := range moduleSizes {
		keys = append(keys, k)
	}
	sort.Ints(keys)

	for _, commID := range keys {
		nNodes := moduleSizes[commID]
		nEdges := moduleEdges[commID]
		dens := moduleDensities[commID]

		record := []string{
			strconv.Itoa(commID),
			strconv.Itoa(nNodes),
			strconv.Itoa(nEdges),
			strconv.FormatFloat(dens, 'f', 6, 64),
		}
		if err := w.Write(record); err != nil {
			return err
		}
	}

	return nil
}

// *ChatGPT*
// WriteNodeStatsCSV writes per-node statistics to a CSV file.
// Columns: node_id, gene_name, community_id, degree, clustering
func WriteNodeStatsCSV(
	filename string,
	geneNames []string,
	clusterMap map[int]int,
	degrees []float64,
	clusteringCoeffs []float64,
) error {

	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	w := csv.NewWriter(f)
	defer w.Flush()

	// header
	if err := w.Write([]string{"node_id", "gene_name", "community_id", "degree", "clustering"}); err != nil {
		return err
	}

	n := len(geneNames)
	for i := 0; i < n; i++ {
		commID, ok := clusterMap[i]
		if !ok {
			// if for some reason this node isn't in the clusterMap, skip it
			continue
		}

		deg := 0.0
		if i < len(degrees) {
			deg = degrees[i]
		}

		c := math.NaN()
		if i < len(clusteringCoeffs) {
			c = clusteringCoeffs[i]
		}

		record := []string{
			strconv.Itoa(i),                      // node_id
			geneNames[i],                         // gene_name
			strconv.Itoa(commID),                 // community_id
			strconv.FormatFloat(deg, 'f', 6, 64), // degree
			strconv.FormatFloat(c, 'f', 6, 64),   // clustering
		}

		if err := w.Write(record); err != nil {
			return err
		}
	}

	return nil
}
