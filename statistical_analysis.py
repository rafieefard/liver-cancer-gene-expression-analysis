"""
Statistical Analysis Module

This module provides the StatisticalAnalysis class for 
performing statistical analysis on gene expression data.
It includes methods for calculating mean, median, standard deviation, 
and differential between HCC and normal samples.
"""

from math import sqrt
from gene_expression_data import GeneExpressionData

class StatisticalAnalysis:
    """A class for performing statistical analysis on gene expression data.
    Attributes:
        gene_expr_data (GeneExpressionData): A GeneExpressionData object.

    Methods:
        get_gene_mean(genes_list: list, types = ['HCC', 'normal'])
        get_gene_median(genes_list: list)
        get_gene_standard_deviation(genes_list: list)
        get_gene_differential(genes_list: list)
        create_aggregate_stats(gene_mean, gene_median, gene_standard_deviation, gene_differential)
        get_n_highest_diff(threshold_n)
        get_gene_max_exp(genes_list: list)

    """
    def __init__(self, gene_expr_data: GeneExpressionData):
        """Initializes the StatisticalAnalysis object."""
        self.gene_expr_data = gene_expr_data

    def get_gene_mean(self, genes_list: list, sample_types: list) -> dict:
        """Returns a dictionary of gene means for a list of genes."""
        exprs = {}
        for gene in genes_list:
            exprs[gene] = []
            # Aggregate expressions for each sample type, HCC and normal
            for sample_type in sample_types:
                exprs[gene].extend(self.gene_expr_data.get_gene_exp(gene)[sample_type])

        genes_mean = {"Gene name": "Mean"} # Put headers in the dictionary for output
        genes_mean.update({
            gene: round(sum(exprs[gene]) / len(exprs[gene]), 3) for gene in genes_list})
        return genes_mean # Example: {'Gene name': 'Mean', '121_at': 7.012, '117_at': 4.013, ...}

    def get_gene_median(self, genes_list: list) -> dict:
        """Returns a dictionary of gene medians for a list of genes."""
        gene_median = {"Gene name": "Median"} # Put headers in the dictionary for output
        for gene in genes_list:
            gen_exprs = self.gene_expr_data.get_gene_exp(gene)
            # Aggregate expressions for each sample type, HCC and normal
            expressions = sum((gen_exprs[sample_type] for sample_type in gen_exprs), [])
            exprs_sorted = sorted(expressions)
            n = len(sorted(expressions))
            if n % 2 == 1:
                gene_median[gene] = round(exprs_sorted[n // 2], 3)
            else:
                gene_median[gene] = round((exprs_sorted[n // 2 - 1] + exprs_sorted[n // 2]) / 2, 3)
        return gene_median
        # Example: {'Gene name': 'Median', '121_at': 7.012, '117_at': 4.013, ...}

    def get_gene_standard_deviation(self, genes_list: list):
        """Returns a dictionary of gene standard deviations for a list of genes."""
        gene_std = {"Gene name": "Standard dev."} # Put headers in the dictionary for output
        gene_means = self.get_gene_mean(genes_list, ['HCC', 'normal'])
        for gene in genes_list:
            gen_exp = self.gene_expr_data.get_gene_exp(gene)
            # Aggregate expressions for each sample type, HCC and normal
            exprs = sum((gen_exp[sample_type] for sample_type in gen_exp), [])
            mean_value = gene_means[gene]
            squared_diff_sum = sum((expr - mean_value) ** 2 for expr in exprs)
            std = sqrt(squared_diff_sum / len(exprs))
            gene_std[gene] = round(std, 3)
        return gene_std
        # Example: {'Gene name': 'Standard dev.', '121_at': 7.012, '117_at': 4.013, ...}

    def get_gene_differential(self, genes_list):
        """Returns a dictionary of differential between HCC and normal for a list of genes.
        The differential is the ratio of the mean expression of 
        HCC samples to the mean expression of normal samples.
        this approach name is Fold Change, which is easy to calculate and is scale independent
        """
        gene_differential = {"Gene name": "HCC Dif. / normal Dif."}
        mean_hcc = self.get_gene_mean(genes_list, sample_types = ['HCC'])
        mean_normal = self.get_gene_mean(genes_list, sample_types = ['normal'])
        gene_differential.update({
            gene: round(mean_hcc[gene] / mean_normal[gene], 3) for gene in genes_list})
        return gene_differential
        # Example: {'Gene name': 'HCC Dif. / normal Dif.', '121_at': 7.012, '117_at': 4.013, ...}

    def get_gene_max_exp(self, genes_list):
        """Returns a dictionary of gene max expressions for a list of genes."""
        samples = self.gene_expr_data.get_sample_list()
        # put headers in dictionary for output
        gene_max_exp = {"Gene name": "Max Value(Sample name)"}
        for gene in genes_list:
            gen_exp_by_type = self.gene_expr_data.get_gene_exp(gene)
            # Aggregate expressions for each sample type
            expressions = sum(gen_exp_by_type.values(), [])
            max_value = max(expressions)
            idx = expressions.index(max_value) # Get index of max value
            sample_name = samples[idx] # Get corresponding sample
            gene_max_exp[gene] = f"{round(max_value, 3)}({sample_name})"
        return gene_max_exp
        # Example: {'Gene name': 'Max Value(Sample name)',
        # '121_at': '7.012(GSM362958.CEL.gz)', '117_at': '4.013(GSM712539.CEL.gz)', ...}

    def get_n_highest_diff(self, threshold_n):
        """Returns a dictionary of n highest differential 
        between HCC and normal for a list of genes."""
        # Get differential of all genes between HCC and normal
        gene_diffs = self.get_gene_differential(self.gene_expr_data.get_gene_list())
        # Since some genes maybe are ancogene with gene_diff > 1 and some others
        # are TSG with gene_diff < 1, we need to invert their differential for comparison
        gene_diffs_inv = {
            gene: (diff if diff >= 1 else 1 / diff)
            for gene, diff in gene_diffs.items()
            if gene != "Gene name"
        }
        n_highest_genes = {} # dict
        # Extract n highest differential genes
        while threshold_n:
            max_diff_gene_name = max(gene_diffs_inv.keys(), key = lambda x: gene_diffs_inv[x])
            n_highest_genes[max_diff_gene_name] = gene_diffs[max_diff_gene_name]
            gene_diffs_inv.pop(max_diff_gene_name)
            threshold_n -= 1
        return n_highest_genes

    def filter_expr_by_value(self, threshold_value):
        """Returns a dictionary of filtered data by value for a list of genes."""
        filtered_data_by_gene = {} # dict
        sample_list = self.gene_expr_data.get_sample_list() # Get list of all sample names
        for gene in self.gene_expr_data.get_gene_list(): # Iterate over all genes
            gene_exprs = self.gene_expr_data.data["by_gene"][gene]
            # Get expressions more than threshold for each gene
            high_expressions = [
                (sample_list[i], round(expr, 3))
                for i, expr in enumerate(gene_exprs)
                if expr > threshold_value
            ] # Example: [('GSM362958.CEL.gz', 7.012), ('GSM712539.CEL.gz', 4.013), ...]

            # Add gene and its expressions if high_expressions is not empty
            if high_expressions:
                filtered_data_by_gene[gene] = high_expressions

        # Sort expressions of each gene by value in descending order
        sorted_filtered_data_by_gene = {
            gene_name: sorted(sam_exp_list, key=lambda x:x[1], reverse=True)
            for gene_name, sam_exp_list in filtered_data_by_gene.items()
        }
        return sorted_filtered_data_by_gene
        # Example: {
        # '1007_s_at': [('GSM362960.CEL.gz', 7.803), ('GSM362959.CEL.gz', 7.586)],
        # '117_at': [('GSM362958.CEL.gz', 6.788)],
        #        }

    def get_sample_mean(self, samples_list):
        """Returns a dictionary of sample means for a list of samples."""
        samples_mean = {"Sample name": "Mean"} # Put headers in dictionary for output
        for sample in samples_list:
            expressions = self.gene_expr_data.data['by_sample'][sample][1]
            samples_mean[sample] = round(sum(expressions) / len(expressions), 3)
        return samples_mean
        # Example: {'Sample name': 'Mean', 'GSM3958.CEL.gz': 7.012, 'GSM719.CEL.gz': 4.013, ...}

    def get_sample_median(self, samples_list: list) -> dict:
        """Returns a dictionary of sample medians for a list of samples."""
        sample_median = {"Sample name": "Median"} # Put headers in dictionary for output
        for sample_name in samples_list:
            expressions = self.gene_expr_data.get_sample_exp(sample_name)
            expr_sorted = sorted(expressions)
            n = len(expr_sorted)
            if n % 2 == 1:
                sample_median[sample_name] = round(expr_sorted[n // 2], 3)
            else:
                sample_median[sample_name] = round(
                    (expr_sorted[n // 2 - 1] + expr_sorted[n // 2]) / 2, 3)
        return sample_median
        # Example: {'Sample name': 'Median', 'GSM3958.CEL.gz': 7.012, 'GSM719.CEL.gz': 4.013, ...}

    def get_sample_standard_deviation(self, samples_list: list):
        """Returns a dictionary of sample standard deviations for a list of samples."""
        sample_std = {"Sample name": "Standard dev."} # Put headers in dictionary for output
        sample_means = self.get_sample_mean(samples_list)
        for sample_name in samples_list:
            exprs = self.gene_expr_data.get_sample_exp(sample_name)
            mean_value = sample_means[sample_name]
            squared_diff_sum = sum((expr - mean_value) ** 2 for expr in exprs)
            std = sqrt(squared_diff_sum / len(exprs))
            sample_std[sample_name] = round(std, 3)
        return sample_std
        # Example: {'Sample name': 'Standard dev.', 'GSM358.CEL.gz': 0.0, 'GSM79.CEL.gz': 0.0, ...}

    def get_sample_max_exp(self, samples_list):
        """Returns a dictionary of sample max expressions for a list of samples."""

        # Put headers in dictionary for output
        sample_max_exp = {"Sample name": "Max value(Gene name)"}
        genes = self.gene_expr_data.get_gene_list() # Get list of all gene names
        for sample in samples_list:
            sample_exprs = self.gene_expr_data.data['by_sample'][sample][1]
            max_value = max(sample_exprs)
            idx = sample_exprs.index(max_value) # Get index of max value
            gene_name = genes[idx] # Get corresponding gene
            sample_max_exp[sample] = f"{round(max_value, 3)}({gene_name})"
        return sample_max_exp
        # {'Sample name': 'Max value(Gene name)',
        # 'GSM362958.CEL.gz': '6.801(1007_s_at)',
        # 'GSM362959.CEL.gz': '7.586(117_at)',
        # 'GSM712538.CEL.gz': '7.291(1007_s_at)'
        # }

    def create_aggregate_stats(self, *stats) -> dict:
        """Returns a dictionary of aggregated statistics for a list of statistics."""
        aggregated_stat = {} # dict
        instance_names = stats[0].keys() # Get genes or samples names
        for instance in instance_names:
            aggregated_stat[instance] = [stat[instance] for stat in stats]
        return aggregated_stat
        # Example: {
        # 'Gene name': ['Mean', 'Median', 'Standard dev.', 'HCC Dif. / normal Dif.', 'Max Value(Sample name)'],
        # '117_at': [4.147, 3.777, 1.011, 1.149, '6.788(GSM362958.CEL.gz)'],
        # '1255_g_at': [3.257, 3.263, 0.127, 1.063, '3.477(GSM362960.CEL.gz)'],
        # }
        # Or
        # {'Sample name': ['Mean', 'Median', 'Standard dev.', 'Max value(Gene name)'],
        # 'GSM362958.CEL.gz': [5.216, 5.431, 1.399, '6.801(1007_s_at)'],
        # 'GSM362959.CEL.gz': [4.986, 4.194, 1.51, '7.586(1007_s_at)'],
        # }
