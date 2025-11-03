
// ParseGeneExpressionFile reads a tab-separated gene expression file and returns
// a map from gene_name -> tpm_unstranded for rows where gene_type == "protein_coding"

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
