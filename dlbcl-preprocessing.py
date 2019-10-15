import copy
import math
import random
import progressbar
import statistics
import snr
import zscore
import sort

raw_data_matrix = list()
data_matrix = list()
normalized_data_matrix = list()
selected_data_matrix = list()
gene_attributes = list()
selected_gene_attributes = list()
snr_tuples = list()
attribute_selection_count = 100
DLBCL_count = 0
FL_count = 0
progressbar_total = 9

def take_user_input():

    global attribute_selection_count
    global progressbar_total

    selection_count = input("Enter the number of gene attributes to select (default is " + str(attribute_selection_count) + "): ")

    if not selection_count:
        selection_count = "100"

    attribute_selection_count = int(selection_count)

    progressbar.show(0, progressbar_total, prefix = 'Progress:', suffix = 'Complete', length = 50)

def read_dlbcl_raw_dataset():

    global raw_data_matrix
    global progressbar_total

    with open("DLBCL.txt", 'r') as datafile:
        for line in datafile:
            line = line.rstrip()
            splitted_line_list = line.split("\t")
            raw_data_matrix.append(splitted_line_list)

    progressbar.show(1, progressbar_total, prefix = 'Progress:', suffix = 'Complete', length = 50)

def tidy_raw_dataset():

    global raw_data_matrix
    global data_matrix
    global gene_attributes

    gene_attributes = raw_data_matrix[0][:-1]
    data_matrix = raw_data_matrix[3:]

    progressbar.show(2, 9, prefix = 'Progress:', suffix = 'Complete', length = 50)

def write_as_csv():

    global data_matrix
    global gene_attributes
    global progressbar_total

    filename = "dlbcl-fl.csv"

    writefile = open(filename, 'w+')
    write_file_content = ""

    for attribute in gene_attributes:
        write_file_content += attribute + ","

    write_file_content += "class\n"

    for sample in data_matrix:

        line = ""

        for value in sample[:-1]:
            line += str(value) + ","

        line += sample[-1] + "\n"
        write_file_content += line

    writefile.write(write_file_content)
    writefile.close()

    progressbar.show(3, progressbar_total, prefix = 'Progress:', suffix = 'Complete', length = 50)

def count_class_strength():

    global data_matrix
    global DLBCL_count
    global FL_count
    global progressbar_total

    for sample in data_matrix:
        if sample[-1] == "DLBCL":
            DLBCL_count += 1
        elif sample[-1] == "FL":
            FL_count += 1

    progressbar.show(4, progressbar_total, prefix = 'Progress:', suffix = 'Complete', length = 50)
    
def convert_datapoints_to_number():

    global data_matrix
    global progressbar_total

    for sample in data_matrix:
        for index in range(len(sample) - 1):
            sample[index] = float(sample[index])

    progressbar.show(5, progressbar_total, prefix = 'Progress:', suffix = 'Complete', length = 50)

def normalize_data():

    global data_matrix
    global progressbar_total

    sample_length = len(data_matrix[0])
    sample_count = len(data_matrix)

    for attribute_index in range(sample_length - 1):

        attribute_list = list()
        normalized_attribute_list = list()

        for sample_index in range(sample_count):
            attribute_list.append(data_matrix[sample_index][attribute_index])

        for attribute in attribute_list:
            
            z_score = zscore.calculate_zscore(attribute, attribute_list)
            rounded_zscore = math.ceil(z_score * 10000) / 10000
            normalized_attribute_list.append(rounded_zscore)

        for sample_index in range(sample_count):
            data_matrix[sample_index][attribute_index] = normalized_attribute_list[sample_index]

    progressbar.show(6, progressbar_total, prefix = 'Progress:', suffix = 'Complete', length = 50)

def sort_by_SNR():

    global data_matrix
    global DLBCL_count
    global snr_tuples
    global progressbar_total

    sample_length = len(data_matrix[0])
    sample_count = len(data_matrix)

    for attribute_index in range(sample_length - 1):

        attribute_list = list()

        for sample_index in range(sample_count):
            attribute_list.append(data_matrix[sample_index][attribute_index])

        snr_value = snr.mod_SNR(attribute_list[:DLBCL_count], attribute_list[DLBCL_count:])
        rounded_snr = math.ceil(snr_value * 1000) / 1000
        snr_tuples.append((attribute_index, rounded_snr))

    sort.randomized_quick_sort_for_tuples(snr_tuples, 0, len(snr_tuples) - 1)

    progressbar.show(7, progressbar_total, prefix = 'Progress:', suffix = 'Complete', length = 50)

def prepare_selected_dataset():

    global data_matrix
    global gene_attributes
    global selected_gene_attributes
    global snr_tuples
    global attribute_selection_count
    global progressbar_total

    flag = 0

    for sample in data_matrix:

        selected_attribute_sample = list()

        for index in range(attribute_selection_count):

            snr_tuple = snr_tuples[index]
            selected_index = snr_tuple[0]

            if flag == 0:
                selected_gene_attributes.append(gene_attributes[selected_index])

            selected_attribute_sample.append(sample[selected_index])

        selected_attribute_sample.append(sample[-1])
        selected_data_matrix.append(copy.deepcopy(selected_attribute_sample))
        
        flag = 1

    progressbar.show(8, progressbar_total, prefix = 'Progress:', suffix = 'Complete', length = 50)

def write_to_file():

    global selected_data_matrix
    global selected_gene_attributes
    global attribute_selection_count
    global progressbar_total

    filename = "dlbcl-selected-" + str(attribute_selection_count) + ".csv"

    writefile = open(filename, 'w+')
    write_file_content = ""

    for attribute in selected_gene_attributes:
        write_file_content += attribute + ","

    write_file_content += "class\n"

    for sample in selected_data_matrix:

        line = ""

        for value in sample[:-1]:
            line += str(value) + ","

        line += sample[-1] + "\n"
        write_file_content += line

    writefile.write(write_file_content)
    writefile.close()

    progressbar.show(9, progressbar_total, prefix = 'Progress:', suffix = 'Complete', length = 50)
    print("The pre-processed data has been saved in the file: ", filename, "!", end = "\n")

def main():

    take_user_input()
    read_dlbcl_raw_dataset()
    tidy_raw_dataset()
    write_as_csv()
    count_class_strength()
    convert_datapoints_to_number()
    normalize_data()
    sort_by_SNR()
    prepare_selected_dataset()
    write_to_file()

    return

main()