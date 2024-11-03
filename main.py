"""
A Gene Expression Analysis Tool for Liver Cancer Data

This script helps you analyze gene expression data from liver cancer samples. 
It processes gene expression data from a CSV file and generate some cool stats 

Input arguments:
    1. Path to liver cancer data CSV
    2. Comma-separated gene names
    3. Comma-separated sample names
    4. Number of top differentially expressed genes to find
    5. Expression value threshold for filtering
    6. Width of the output display
    7. Output destination ("Display" or path to output file)
    8. (Optional) Result file name if "Display" was chosen

Example:
    python main.py \
        "./Data/Liver_GSE14520_U133A.csv" "1007_s_at, 1053_at" \
            "GSM362958.CEL.gz, GSM362959.CEL.gz" 10 14 155 ./result.txt

The program does the following tasks:
- Loads the liver cancer gene expression data
- Calculates various statistics (mean, median, standard deviation)
- Finds genes that are expressed differently in cancer vs normal samples
- Shows the n top differentially expressed genes
- Outputs everything in a nice, readable format

The program works with the following classes:
- GeneExpressionData: Handles the raw data
- StatisticalAnalysis: Calculates various statistics
- AnalysisReport: Reports the results
"""

import sys
from unittest.mock import DEFAULT
from statistical_analysis import StatisticalAnalysis
from gene_expression_data import GeneExpressionData
from analysis_report import AnalysisReport
from constants import *

def analize_args(args):
    """Handles the command line arguments."""
    # Check the number of arguments
    if len(args) != 8 and len(args) != 9:
        print("Usage: python main.py liver_file_path gene_names sample_names "
              "n_threshold value_threshold width destination [result_file]")
        sys.exit(1)
    liver_file_path = args[1] # Get the path to the liver cancer data
    if not liver_file_path.endswith('.csv'): # Check if it's a CSV file
        print("Input file must be a CSV file.")
        sys.exit(1)
    input_genes = args[2] # Get the gene names
    if not input_genes: # Check if there are any gene names
        print("No gene names provided.")
        sys.exit(1)
    # Put the gene names in a list and clean them
    gene_names = [gene.strip() for gene in input_genes.split(",")]
    input_samples = args[3] # Get the sample names
    if not input_samples:
        print("No sample names provided.")
        sys.exit(1)
    # Put the sample names in a list and clean them
    sample_names = [sample.strip() for sample in input_samples.split(",")]
    # get the threshold for the most differentially expressed genes and check its integrity
    # if it's not a valid number, set it to DEFAULT_N_THRESHOLD = 10
    try:
        n_threshold = int(args[4])
    except ValueError:
        print("n_threshold must be an integer.")
        print("Program will continue with default value of 10.")
        n_threshold = DEFAULT_N_THRESHOLD
    if n_threshold < 0:
        print("n_threshold must be a positive integer.")
        print("Program will continue with default value of 10.")
        n_threshold = DEFAULT_N_THRESHOLD
    # get the threshold for the most gene values and check its integrity
    # if it's not a valid number, set it to DEFAULT_VALUE_THRESHOLD = 14
    try:
        value_threshold = float(args[5])
    except ValueError:
        print("value_threshold must be a float.")
        print("Program will continue with default value of 14.")
        value_threshold = DEFAULT_VALUE_THRESHOLD
    if value_threshold < 0:
        print("value_threshold must be a positive float.")
        print("Program will continue with default value of 14.")
        value_threshold = DEFAULT_VALUE_THRESHOLD
    # get the width of the display and check its integrity
    # if it's not a valid number, set it to DEFAULT_WIDTH = 105
    try:
        width = int(args[6])
    except ValueError:
        print("width must be an integer number.")
        print("Program will continue with default value of 105.")
        width = DEFAULT_WIDTH
    if width < 5:
        print("width must be at least 95.")
        print("Program will continue with default value of 95.")
        width = MIN_WIDTH
    destination = args[7] # get the output destination
    if destination.capitalize() == "Display":
        # Check if it's a valid txt file
        try:
            result_file = args[8]
            if not result_file.endswith('.txt'):
                print("Result file must be a .txt file.")
                print("Program will continue without a result file.")
                result_file = ""
        except IndexError:
            result_file = ""
    else:
        result_file = destination
        destination = "File"
        if not destination.endswith('.txt'):
            print("Output file must be a .txt file.")
            print("Program will continue without a result "
                  "file and showing the report in the display.")
            destination = "Display"
    return (liver_file_path, gene_names, sample_names, n_threshold,
            value_threshold, width, destination, result_file)

def do_gene_analysis(stat, gene_names):
    """
    Perform statistical analysis on genes.
    """
    # calculate the mean of list of genes
    gene_means = stat.get_gene_mean(gene_names, ['HCC', 'normal'])
    gene_median = stat.get_gene_median(gene_names)
    gene_standard_deviation = stat.get_gene_standard_deviation(gene_names)
    gene_differential = stat.get_gene_differential(gene_names)
    max_gene_exp = stat.get_gene_max_exp(gene_names)

    return stat.create_aggregate_stats(
        gene_means, gene_median, gene_standard_deviation,
        gene_differential, max_gene_exp
    )

def do_sample_analysis(stat, sample_names):
    """
    Perform statistical analysis on samples.

    Args:
        stat (StatisticalAnalysis): Statistical analysis object
        sample_names (list): List of sample names to analyze

    Returns:
        dict: Aggregated sample statistics
    """
    max_sample_exp = stat.get_sample_max_exp(sample_names)
    sample_mean = stat.get_sample_mean(sample_names)
    sample_median = stat.get_sample_median(sample_names)
    sample_standard_deviation = stat.get_sample_standard_deviation(sample_names)

    return stat.create_aggregate_stats(
        sample_mean, sample_median,
        sample_standard_deviation, max_sample_exp
    )

def main():
    """
    Main function that coordinates the analysis.

    Takes command line arguments, sets up the analysis objects.
    """

    # Handle command line arguments
    (liver_file_path, gene_names, sample_names, n_threshold,
     value_threshold, width, destination, result_file) = analize_args(sys.argv)

    # Create objects
    report = AnalysisReport(destination, result_file, width)
    gene_expr_data = GeneExpressionData(liver_file_path)
    stat = StatisticalAnalysis(gene_expr_data)

    # Calculate statistical data
    # Get Mean, Median, SD, Diff, Max for each user gene
    gene_stats = do_gene_analysis(stat, gene_names)
    # Get Mean, Median, SD, Max for each user sample
    sample_stats = do_sample_analysis(stat, sample_names)
    # Get n highest differential genes
    n_high_diff_genes = stat.get_n_highest_diff(n_threshold)
    # Get combinations of genes and samples with expression above threshold
    filtered_gene_by_value = stat.filter_expr_by_value(value_threshold)
    # Get list of all genes
    genes_list = gene_expr_data.get_gene_list()
    # Get list of all samples
    samples_list = gene_expr_data.get_sample_list()

    # Report user genes statistics
    gene_stats_header = "Stats of genes"
    gene_stats_footer = "End of genes stats"
    report.write_stat("Gene", gene_stats_header, gene_stats, gene_stats_footer)

    # Report user samples statistics
    sample_stats_header = "Stats of samples"
    sample_stats_footer = "End of samples stats"
    report.write_stat("Sample", sample_stats_header, sample_stats, sample_stats_footer)

    # Report n highest differentially expressed genes
    n_highest_genes_header = f" \
    {n_threshold} genes that are expressed most differently in normal versus HCC."
    n_highest_genes_footer = "End of filtered genes"
    report.write_n_highest_genes(
        n_highest_genes_header, n_high_diff_genes, n_highest_genes_footer)

    # Report combinations of genes and samples with expression above threshold
    filtered_gene_header = f" {len(filtered_gene_by_value)} combinations of genes and samples "
    f"which expressed above for each gene {value_threshold} in order of expression value"
    filtered_gene_footer = "End of filtered genes"
    report.write_filtered_gene_by_value(
        filtered_gene_header, filtered_gene_by_value, filtered_gene_footer)

    # Report list of all genes
    genes_list_header = "List of genes"
    genes_list_footer = "\nEnd of list of genes"
    report.write_list(genes_list_header, genes_list, genes_list_footer)

    # Report list of all samples
    samples_list_header = "List of samples"
    samples_list_footer = "\nEnd of list of samples"
    report.write_list(samples_list_header, samples_list, samples_list_footer)


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError:
        print(f"File not found at {sys.argv[1]}")
        sys.exit()
