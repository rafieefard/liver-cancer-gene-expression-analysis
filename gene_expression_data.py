"""
Gene Expression Data Module

This module provides the GeneExpressionData class for loading and managing gene expression data
from liver cancer samples. It reads CSV files containing gene expression values
and provides methods to access data by gene and sample.
"""

from collections import defaultdict

class GeneExpressionData:
    """
    Attributes:
        liver_file_path (str): The path to the liver cancer data file.
        data (dict): A dictionary containing data by gene and sample.
        example: {'by_gene': {[samples], 'type': [types], gene: [expressions]}, \\
        'by_sample': {[genes], sample: [expressions]}}

    Methods:
        get_gene_list(): Returns a list of gene names.
        get_gene_exp(gene: str): Returns a dictionary of gene expressions for a given gene, 
        by HCC and normal samples.
        Example: {"HCC": [expressions], "normal": [expressions]}
        get_sample_exp(sample: str): Returns a list of expressions for a given sample.
        get_sample_list(): Returns a list of sample names.

    Example:
        data = GeneExpressionData("./Data/Liver_GSE14520_U133A.csv")
        gene_list = data.get_gene_list()
        gene_expressions = data.get_gene_exp("1007_s_at")
        sample_expressions = data.get_sample_exp("GSM362958.CEL.gz")
    """
    def __init__(self, liver_file_path):
        self.liver_file_path = liver_file_path
        self.data = self.read_cancer_data(liver_file_path)

    @staticmethod
    def read_cancer_data(liver_file_path: str) -> dict:
        """Reads the cancer data file and returns a dictionary of gene expression data."""
        data_by_samples = {} # Example: {sample: (sample_type, [expressions])}
        data_by_gene = defaultdict(list) # Example: {gene: [expressions]}
        with open(liver_file_path, "r", encoding="UTF-8") as cancer_file:
            print("Reading file...")
            first_line = cancer_file.readline().strip()
            headers = first_line.split(",") # put the headers in a list
            for line in cancer_file:
                row_data = line.split(",")
                sample = row_data[0]
                sample_type = row_data[1]
                # first 2 columns are sample name and type
                exprs_by_sample = [float(expr) for expr in row_data[2:]]
                data_by_samples[sample] = (sample_type, exprs_by_sample)
                for col_no, header in enumerate(headers):
                    # Handle first 2 columns which are sample name and type not float values
                    try:
                        data_by_gene[header].append(float(row_data[col_no]))
                    except ValueError:
                        data_by_gene[header].append(row_data[col_no])
            data_by_gene = dict(data_by_gene) # Convert defaultdict to regular dict
            print("done.")
            data = {'by_gene': data_by_gene, 'by_sample': data_by_samples}
        return data
        # Example:
        # {'by_gene':{
        # 'samples': [samples], 'type': [types], gene_name: [expressions]
        #           },
        # 'by_sample:{
        # 'samples': [gene_names], sample_name: [expressions]
        #           }
        # }

    def get_gene_list(self) -> list:
        """Returns a list of gene names."""
        return list(self.data['by_gene'].keys())[2:]

    def get_gene_exp(self, gene: str) -> dict:
        """Returns a dictionary of gene expressions for a given gene, by HCC and normal samples."""
        gene_expressions = defaultdict(list)
        for row, exp in enumerate(self.data['by_gene'][gene]):
            sample_type = self.data['by_gene']['type'][row]
            gene_expressions[sample_type].append(exp)
        gene_expressions = dict(gene_expressions) # Convert defaultdict to regular dict
        return gene_expressions # Sample:{"HCC": [expressions], "normal": [expressions]}

    def get_sample_exp(self, sample_name: str) -> list:
        """Returns a list of expressions for a given sample."""
        sample_expressions = self.data['by_sample'][sample_name][1]
        return sample_expressions

    def get_sample_list(self) -> list:
        """Returns a list of sample names."""
        return list(self.data['by_sample'].keys())
    