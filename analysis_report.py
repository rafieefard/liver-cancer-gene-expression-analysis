"""This module contains a class for writing analysis reports.

The class provides methods for writing header, footer, and contents
of analysis reports.
"""
import sys
from time import time
class AnalysisReport:
    """Class for generating formatted analysis reports

    The class provides methods for writing header, footer, and contents
    of analysis reports.

    Attributes:
        destination (str): The destination of the report
        result_file (str): The path of the result file
        width (int): The width of the report

    Methods:
        wrap_text(self, text)
        create_empty_file(self)
        write_content(self, content)
        write_header(self, header)
        write_footer(self, footer)
        write_list(self, header, result, footer)
        write_n_highest_genes(self, header, result, footer)
        write_stat(self, instance, header, results, footer)
        write_filtered_gene_by_value(self, header, result, footer)
    """
    def __init__(self, destination: str, result_file: str, width = 170):
        self.destination = destination
        self.output_file_path = destination
        self.width = width
        self.result_file = result_file
        self.create_empty_file()

    def wrap_text(self, text):
        """Wrap text into lines of a given width."""
        words = text.split(" ")
        lines = []
        line= ""
        for word in words:
            next_word = word + " "
            if len(line + next_word) < self.width - 4: # 2 characters margin from each side
                line += next_word
            else:
                lines.append(line)
                line = next_word
        if sum(len(line) for line in lines) < len(text): # Add last line to list of lines
            lines.append(line)
        lines[-1]= lines[-1][:-1] # Remove space character from end of the last line
        return lines

    def create_empty_file(self):
        """Clean previous file content for new content and ready it to write"""
        if self.result_file: # Open file if result_file is not empty
            # Check if the file is use by another program and rename it
            try:
                with open(self.result_file, "w", encoding="utf-8") as report_file:
                    report_file.write("")
            except PermissionError:
                self.result_file = f"{self.result_file[:-4]}_{str(int(time()))}.txt"
                print("The report file is being used by another program.")
                print(f"The report will be written in {self.result_file}")

    def write_content(self, content):
        """Write content to destination."""
        if self.destination.capitalize() == "Display":
            print(content)
        if self.result_file:
            # Open file to add new content to it
            with open(self.result_file, "a", encoding="utf-8") as report_file:
                report_file.write(content + "\n")

    def write_header(self, header):
        """Write header in a box created with asterisks to destination."""
        wrapped_text = self.wrap_text(header)
        output = "*" * self.width + "\n"
        for line in wrapped_text:
            output += "* " + line.center(self.width - 4) + " *\n"
        output += "*" * self.width
        self.write_content(output) # ready header for use by write_content function

    def write_footer(self, header):
        """Write footer to destination."""
        wrapped_text = self.wrap_text(header)
        output = ""
        for line in wrapped_text:
            output += line.center(self.width - 4) + " \n"
        output += "=" * self.width + "\n\f"
        self.write_content(output) # ready footer for use by write_content function

    def write_stat(self, instance, header, results, footer):
        """Write statistics in several columns to destination."""
        self.write_header(header)
        output = ""
        # Get the number of padding between columns instance can be "Gene" or "Sample"
        padding_no = len(results[f"{instance} name"])
        # results = {'Gene name': ['Mean', 'Median', 'Standard dev.', 'Max value(Gene name)'], ...
        # Put lenght of header titles for each column to list as primary value
        instance_cols_width = \
            [len(list(results.items())[0][0])] + [len(col) for col in list(results.items())[0][1]]
        # Iterate over results items to find maximum content length for each column
        # and update columns width
        for instance, stat_value_list in results.items():
            if len(instance) > instance_cols_width[0]: # Check first column (Instance name)
                instance_cols_width[0] = len(instance)
            for idx, expr in enumerate(stat_value_list): # Check stats value columns
                if len(str(expr)) > instance_cols_width[idx + 1]:
                    instance_cols_width[idx + 1] = len(str(expr))
        total_stats_col_width = sum(instance_cols_width)

        # Calculate distance between columns which can not less than 4
        padding = max((self.width - total_stats_col_width) // padding_no, 4)
        # Total space for each column
        stats_col_width = [width + padding for width in instance_cols_width]

        for idx, instance_name_stats in enumerate(results.items()):
            instance_name = instance_name_stats[0]
            stats = instance_name_stats[1]
            # Put gene or sample name at the center of first column for each row
            instance_name_str = instance_name.center(instance_cols_width[0], " ")
            # Put statistical values at the center of columns for each row
            str_stats = [
                str(value).center(stats_col_width[i + 1], " ") for i, value in enumerate(stats)
            ]
            concat_stats = "".join(str_stats)
            output += f"{instance_name_str}{concat_stats}\n"
            output += f"{(sum(stats_col_width) - padding) * '-'}\n" # Draw line between rows
        output = output[:-1]
        self.write_content(output)
        self.write_footer(footer)

    def write_n_highest_genes(self, header, result, footer):
        """Write n highest genes to destination."""
        self.write_header(header)
        output, output_for_line = "", ""
        max_gene_name_width = max((len(gene) for gene in result)) # Get maximum gene name length

        # Genes with differential more than 1 are Overexpressed and probabely oncogene
        # and Genes with differential less than 1 are Underexpressed and probabely TSG
        for gene,diff in result.items():
            output_for_line += f'Gene "{gene.center(max_gene_name_width, " ")}" with differential '
            output_for_line += f"{str(diff).center(5, ' ')} is "
            output_for_line += f"{'Overexpressed ' if diff > 1 else 'Underexpressed'} "
            output_for_line += "and probabely is a "
            output_for_line += f"{'oncogene' if diff > 1 else 'tumor suppressor gene(TSG)'}"
            output_for_line = output_for_line.center(self.width, " ")
            output += output_for_line + "\n\n"
            output_for_line = ""

        output += "-" * (len(output) // len(result))
        self.write_content(output)
        self.write_footer(footer)

    def write_filtered_gene_by_value(self, header, result, footer):
        """Write filtered gene by value to destination."""
        self.write_header(header)
        output = ""
        # Calculate column widths
        gene_width = max((len(gene) for gene in result)) + 5
        # Find maximum sample name length and add padding
        sample_width = max((len(sample_exp[0]) for sample_exp_list in result.values() \
                            for sample_exp in sample_exp_list)) + 10
        expr_width = 20 # Fixed width for expression values
        total_width = gene_width + sample_width + expr_width
        output = ""
        # Iterate through each gene and its sample expressions
        for i, gene_sample_exp_list in enumerate(result.items()):
            gene_name = gene_sample_exp_list[0]
            sample_exp_list = gene_sample_exp_list[1]
            # Create section header for each gene
            com_no = f" Combination No. {i + 1} with gene \"{gene_name}\" "
            output += f"\n{com_no.center(total_width, '/')}\n"

            # Create table header
            output += f"{total_width * '-'}\n"
            output += f"{'Gene name'.center(gene_width, ' ')}"
            output += f"{'Sample name'.center(sample_width, ' ')}"
            output += f"{'Expression value'.center(expr_width, ' ')}\n"
            output += f"{total_width * '-'}\n"

            # Add each sample and its expression value
            for sample_exp in sample_exp_list:
                output += f"{gene_name.center(gene_width, ' ')}"
                output += f"{sample_exp[0].center(sample_width, ' ')}"
                output += f"{str(sample_exp[1]).center(expr_width, ' ')}"
                output += f"\n{total_width * '-'}\n"
        self.write_content(output)
        self.write_footer(footer)

    def write_list(self, header, result, footer):
        """Write list to destination."""
        self.write_header(header)
        output = ""
        # Lambda function to calculate padding between columns
        padding = lambda n: (self.width - (max((len(str(item)) for item in result))) * n) // (n - 1)

        # Calculate maximum possible number of columns based on width
        n_max = self.width // max((len(str(item)) for item in result))

        # Calculate padding between columns
        col_padding = padding(n_max)

        # Reduce number of columns until padding is at least 5 spaces
        while col_padding < 5:
            n_max -= 1
            col_padding = padding(n_max)
        n = n_max
        col_width = self.width // n
        string_result = [str(item) for item in result]

        # Format output in columns
        for i, item in enumerate(string_result):
            output += item + " " * (col_width - len(str(item)))
            # Add newline after each row is complete
            if (i + 1) % n == 0:
                output += "\n"
        self.write_content(output)
        self.write_footer(footer)
