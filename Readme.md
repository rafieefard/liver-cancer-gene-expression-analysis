# Gene Expression Analysis Tool for Liver Cancer

A Python tool for analyzing gene expression data from liver cancer samples to identify differentially expressed genes, particularly oncogenes and tumor suppressor genes, using the fold change approach.

## Description

This tool analyzes gene expression data from liver cancer (HCC) and normal samples to identify genes that show significant differences in expression. It uses the fold change approach (HCC/normal ratio) which is the standard in the field according to Ashekyan et al. (2023) because:

### Overview
A comprehensive Python tool designed for analyzing gene expression patterns in liver cancer, specifically focusing on identifying differentially expressed genes between Hepatocellular Carcinoma (HCC) and normal tissue samples. The tool implements object-oriented programming principles to process, analyze, and visualize gene expression data.

#### Fold Change Approach
The tool uses the fold change (FC) approach to calculate differential expression:
```python
FC = HCC_mean / normal_mean
```
This approach was chosen based on recent research (Ashekyan et al., 2023) because:
1. Scale Independence
   - Works equally well for both highly and lowly expressed genes
   - 2-fold change means the same thing whether expression values are large or small
   - Example: FC=2 has same biological significance for values (100 vs 50) and (10 vs 5)

2. Biological Interpretation
   - FC > 1: Indicates overexpression in cancer (potential oncogenes)
   - FC < 1: Indicates underexpression in cancer (potential tumor suppressors)
   - Directly interpretable as "X times more/less expressed"

3. Industry Standard
   - Common in gene expression studies
   - Well-understood by biologists
   - Consistent with literature conventions

#### Gene Classification
The tool identifies:

1. Potential Oncogenes:
   - FC > 1 (overexpressed in cancer)
   - Example: FC = 2.0 means gene is expressed 2x higher in cancer
   - These genes might promote cancer growth when overexpressed

2. Potential Tumor Suppressor Genes (TSG):
   - FC < 1 (underexpressed in cancer)
   - Example: FC = 0.5 means gene is expressed 2x lower in cancer
   - These genes might normally prevent cancer when properly expressed

### Data Processing Pipeline

1. Data Loading
   - Reads CSV file containing gene expression values
   - Organizes data by both genes and samples
   - Handles both HCC and normal tissue samples

2. Statistical Analysis
   - Calculates basic statistics (mean, median, standard deviation)
   - Computes differential expression using fold change
   - Identifies most differentially expressed genes
   - Example calculation:
     ```python
     # For a gene that is 2x overexpressed in cancer:
     HCC_mean = 100
     normal_mean = 50
     fold_change = 100/50 = 2.0  # Potential oncogene
     
     # For a gene that is 2x underexpressed in cancer:
     HCC_mean = 30
     normal_mean = 60
     fold_change = 30/60 = 0.5  # Potential tumor suppressor
     ```

3. Advanced Analysis Features
   - Ranks genes by magnitude of differential expression
   - Handles both over and underexpression equally
   - For comparison purposes, underexpressed genes (FC < 1) are converted to negative reciprocals
   - Example:
     ```python
     FC = 2.0 → +2.0 (2-fold overexpression)
     FC = 0.5 → -2.0 (2-fold underexpression)
     ```

4. Output Generation
   - Generates formatted reports
   - Provides statistical summaries
   - Lists top differentially expressed genes
   - Shows gene-sample combinations above threshold values

   
### Data Requirements
Input data should be in CSV format with:
- Samples in rows
- Genes in columns
- Type column specifying "HCC" or "normal"
- Expression values as numeric data
- Headers for all columns

## Getting Started

### Dependencies

* Python 3.7+
* No external libraries required
* Windows/Linux/Mac OS compatible

### Installing

```bash
git clone https://github.com/rafieefard
cd gene-expression-tool
```

### Executing program

Run from command line with these arguments:
```bash
python main.py [data_file.csv] [gene_names] [sample_names] [n_threshold] [value_threshold] [width] [destination] [result_file.txt]
```

Example:
```bash
python main.py "./Data/Liver_GSE14520_U133A.csv" "1007_s_at,1053_at" "GSM362958.CEL.gz" 10 14 105 Display results.txt
```

Arguments:
1. Path to liver cancer data CSV
2. Comma-separated gene names
3. Comma-separated sample names
4. Number of top differentially expressed genes to find
5. Expression value threshold for filtering
6. Width of the output display
7. Output destination ("Display" or path to output file)
8. (Optional) Result file name if "Display" was chosen

## Features

### 1. Statistical Analysis

#### Basic Statistics
- **Mean Expression**
  ```python
  mean = sum(expressions) / len(expressions)
  ```
  - Calculated separately for HCC and normal samples
  - Used in fold change calculations
  - Rounded to 3 decimal places for readability

- **Median Expression**
  - Handles both odd and even number of samples
  - Provides robust measure less affected by outliers
  - Implementation:
    ```python
    sorted_expr = sorted(expressions)
    if len(sorted_expr) % 2 == 1:
        median = sorted_expr[len(sorted_expr) // 2]
    else:
        median = (sorted_expr[len(sorted_expr) // 2 - 1] + 
                 sorted_expr[len(sorted_expr) // 2]) / 2
    ```

- **Standard Deviation**
  - Measures expression variability
  - Calculated using the formula:
    ```python
    squared_diff_sum = sum((expr - mean_value) ** 2 for expr in expressions)
    std = sqrt(squared_diff_sum / len(expressions))
    ```

#### Differential Expression Analysis

- **Top Differential Genes**
  - Identifies most significantly changed genes
  - Handles both over/underexpression using comparable values:
    ```python
    comparable_value = value if value >= 1 else -1/value
    ```
  - Returns original fold change values for interpretation

### 2. Data Handling and Organization

#### Gene-centric Analysis
- Organizes expression data by genes
- Separates HCC and normal expressions
- Example structure:
  ```python
  gene_data = {
      "gene1": {
          "HCC": [expr1, expr2, ...],
          "normal": [expr1, expr2, ...]
      }
  }
  ```

#### Sample-centric Analysis
- Organizes data by samples
- Maintains sample type information
- Structure:
  ```python
  sample_data = {
      "sample1": (sample_type, [expressions]),
      "sample2": (sample_type, [expressions])
  }
  ```

### 3. Report Generation

#### Statistical Reports
- Formatted tables with:
  - Gene/Sample names
  - Mean values
  - Median values
  - Standard deviations
  - Differential expression values
  - Maximum expression values

#### Gene Classification Reports
```
****************************************
Gene "GENE1" with differential 2.500 is Overexpressed and probably is an oncogene
Gene "GENE2" with differential 0.400 is Underexpressed and probably is a tumor suppressor gene
****************************************
```

#### Filtered Expression Reports
- Shows gene-sample combinations above threshold
- Organized by genes
- Sorted by expression values
- Example:
  ```
  ========= Combination for Gene "GENE1" =========
  Gene name    Sample name    Expression value
  ------------------------------------------------
  GENE1        SAMPLE1           15.234
  GENE1        SAMPLE2           14.567
  ```

### 4. Error Handling and Validation

#### Input Validation
- File format checking
- Data type validation
- Range checking for numerical inputs
- Example:
  ```python
  if not file_path.endswith('.csv'):
      raise ValueError("Input must be CSV file")
  ```

#### Error Recovery
- Handles file permission issues
- Provides meaningful error messages
- Implements fallback options
- Example:
  ```python
  try:
      with open(result_file, "w") as f:
          f.write("")
  except PermissionError:
      result_file = f"{result_file[:-4]}_{timestamp}.txt"
  ```

### 5. Customization Options

#### Output Formatting
- Adjustable display width
- Configurable decimal places
- Choice of output destination:
  - Display to console
  - Write to file
  - Both simultaneously

#### Analysis Parameters
- Adjustable thresholds for:
  - Number of top genes
  - Expression value filtering
  - Statistical significance
- Customizable gene/sample lists


## File Structure

### Project Structure
```
gene-expression-tool/
│
├── main.py                 # Main program orchestrator
├── gene_expression_data.py # Data handling class
├── statistical_analysis.py # Statistical computations
├── analysis_report.py      # Report generation
├── constants.py           # Configuration constants
└── data/                  # Data directory
    └── GSE14520_U133A.csv # Sample data file
```

### Dependencies
- Python 3.7 or higher
- Core Python libraries only:
  ```python
  import sys
  from math import sqrt
  from collections import defaultdict
  from time import time
  ```
- No external dependencies required for basic functionality

### Class Architecture

#### 1. GeneExpressionData Class
```python
class GeneExpressionData:
    def __init__(self, liver_file_path):
        self.liver_file_path = liver_file_path
        self.data = self.read_cancer_data(liver_file_path)
```
Key Methods:
- `read_cancer_data()`: Parses CSV file
- `get_gene_exp()`: Returns gene expressions
- `get_sample_exp()`: Returns sample expressions

Data Structure:
```python
self.data = {
    'by_gene': {
        'gene1': [expressions],
        'type': [sample_types]
    },
    'by_sample': {
        'sample1': (type, [expressions])
    }
}
```

#### 2. StatisticalAnalysis Class
```python
class StatisticalAnalysis:
    def __init__(self, gene_expr_data):
        self.gene_expr_data = gene_expr_data
```
Key Methods:
- Statistical calculations:
  ```python
  def get_gene_mean(self, genes_list, sample_types):
      exprs = {}
      for gene in genes_list:
          exprs[gene] = []
          for sample_type in sample_types:
              exprs[gene].extend(
                  self.gene_expr_data.get_gene_exp(gene)[sample_type]
              )
      return {
          gene: round(sum(exprs[gene])/len(exprs[gene]), 3)
          for gene in genes_list
      }
  ```
- Differential analysis:
  ```python
  def get_gene_differential(self, genes_list):
      mean_hcc = self.get_gene_mean(genes_list, ['HCC'])
      mean_normal = self.get_gene_mean(genes_list, ['normal'])
      return {
          gene: round(mean_hcc[gene]/mean_normal[gene], 3)
          for gene in genes_list
      }
  ```

#### 3. AnalysisReport Class
```python
class AnalysisReport:
    def __init__(self, destination, result_file, width=170):
        self.destination = destination
        self.result_file = result_file
        self.width = width
```
Output Methods:
- Text formatting:
  ```python
  def wrap_text(self, text):
      words = text.split(" ")
      lines = []
      line = ""
      for word in words:
          if len(line + word + " ") < self.width - 4:
              line += word + " "
          else:
              lines.append(line)
              line = word + " "
      return lines
  ```
- File handling:
  ```python
  def write_content(self, content):
      if self.destination == "Display":
          print(content)
      if self.result_file:
          with open(self.result_file, "a") as f:
              f.write(content + "\n")
  ```

### Data Flow
1. Input Processing:
   ```python
   gene_expr_data = GeneExpressionData(file_path)
   stat = StatisticalAnalysis(gene_expr_data)
   report = AnalysisReport(destination, result_file)
   ```

2. Analysis Pipeline:
   ```python
   gene_stats = stat.create_aggregate_stats(
       stat.get_gene_mean(genes),
       stat.get_gene_median(genes),
       stat.get_gene_standard_deviation(genes),
       stat.get_gene_differential(genes)
   )
   ```

3. Output Generation:
   ```python
   report.write_stat("Gene", header, gene_stats, footer)
   ```

### Error Handling
1. File Operations:
   ```python
   try:
       with open(file_path, "r") as f:
           # File operations
   except FileNotFoundError:
       print(f"File not found: {file_path}")
       sys.exit(1)
   except PermissionError:
       # Handle busy files
       new_file = f"{file_path}_{int(time())}"
   ```

2. Input Validation:
   ```python
   if not isinstance(threshold, (int, float)):
       raise TypeError("Threshold must be numeric")
   if threshold < 0:
       raise ValueError("Threshold must be positive")
   ```

### Performance Considerations
1. Memory Efficiency:
   - Uses generators where appropriate
   - Minimizes data copying
   ```python
   expressions = (
       float(expr) 
       for expr in row_data[2:]
   )
   ```

2. Computation Optimization:
   - Caches frequently used values
   - Pre-calculates statistics when possible
   ```python
   gene_means = self.get_gene_mean(genes_list)  # Calculate once
   for gene in genes_list:
       # Use pre-calculated means
   ```

## Help

### Basic Usage

#### Command Line Syntax
```bash
python main.py [data_file] [genes] [samples] [n] [threshold] [width] [destination] [result_file]
```

#### Parameters Explained
1. **data_file**: Path to CSV file
   ```bash
   "./Data/Liver_GSE14520_U133A.csv"
   ```

2. **genes**: Comma-separated gene names (no spaces)
   ```bash
   "1007_s_at,1053_at,117_at"
   ```

3. **samples**: Comma-separated sample names
   ```bash
   "GSM362958.CEL.gz,GSM362959.CEL.gz"
   ```

4. **n**: Number of top differential genes to find
   ```bash
   10  # Will find 10 most differentially expressed genes
   ```

5. **threshold**: Expression value cutoff
   ```bash
   14.0  # Will filter expressions above 14.0
   ```

6. **width**: Output display width
   ```bash
   155  # Characters per line for formatting
   ```

7. **destination**: Output mode
   ```bash
   "Display"  # Show in console
   "./output.txt"  # Write to file
   ```

8. **result_file**: (Optional) Additional output file when using Display
   ```bash
   "results.txt"
   ```

### Example Usage Scenarios

#### 1. Basic Analysis
```bash
python main.py "./data/liver_data.csv" "1007_s_at" "GSM362958.CEL.gz" 5 14 155 Display
```
Output Example:
```
******************************************
Stats of genes
******************************************
Gene name    Mean    Median    Std Dev    Differential
1007_s_at    7.123   7.000     0.456      2.345
--------------------------------------------------
End of genes stats
==========================================
```

#### 2. Multiple Gene Analysis
```bash
python main.py "./data/liver_data.csv" "1007_s_at,1053_at,117_at" "GSM362958.CEL.gz" 10 14 155 Display
```
Example Output:
```
******************************************
Gene "1007_s_at" with differential 2.345 is Overexpressed and probably is an oncogene
Gene "1053_at" with differential 0.432 is Underexpressed and probably is a tumor suppressor gene
Gene "117_at" with differential 1.789 is Overexpressed and probably is an oncogene
******************************************
```

#### 3. File Output with Detailed Statistics
```bash
python main.py "./data/liver_data.csv" "1007_s_at" "GSM362958.CEL.gz" 5 14 155 results.txt
```

### Common Use Cases

#### 1. Finding Potential Oncogenes
```bash
# Looking for genes with high differential expression (> 1)
python main.py data.csv "all" "all" 20 14 155 Display
```
This will show:
- Top 20 differentially expressed genes
- Both over and underexpressed genes
- Their potential role (oncogene/TSG)

#### 2. Sample Analysis
```bash
# Analyzing specific samples
python main.py data.csv "all" "GSM362958.CEL.gz,GSM362959.CEL.gz" 10 14 155 stats.txt
```
Produces:
- Sample statistics
- Expression profiles
- Comparisons between samples

#### 3. High Expression Filter
```bash
# Finding highly expressed genes
python main.py data.csv "all" "all" 10 20 155 high_expr.txt
```
Shows:
- Genes expressed above threshold 20
- Which samples show high expression
- Sorted by expression value

### Error Handling Examples

#### 1. Invalid File Path
```bash
python main.py "nonexistent.csv" "1007_s_at" "GSM362958.CEL.gz" 5 14 155 Display
```
Output:
```
Error: File not found at nonexistent.csv
```

#### 2. Invalid Parameters
```bash
# Invalid threshold
python main.py data.csv "1007_s_at" "GSM362958.CEL.gz" -5 14 155 Display
```
Output:
```
Error: n_threshold must be a positive integer.
Program will continue with default value of 10.
```

### Best Practices

1. **Data Preparation**
   - Ensure CSV format is correct
   - Check for missing values
   - Verify gene/sample names

2. **Analysis Strategy**
   - Start with known genes
   - Use reasonable thresholds
   - Validate results

3. **Output Management**
   - Use descriptive file names
   - Set appropriate width for readability
   - Keep backups of important results

## Authors

Javad Rafieefard (https://github.com/rafieefard)
javadrafieefard@gmail.com

## Version

* 1.0
    * Initial Release

## Acknowledgments

* Based on methodology from Ashekyan et al. (2023) - "Transcriptomic Maps of Colorectal Liver Metastasis"
* Fold change approach validated by current research in cancer genomics

