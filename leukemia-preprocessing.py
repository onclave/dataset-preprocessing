import copy
import math
import random

raw_data_matrix = list()
data_matrix = list()
normalized_data_matrix = list()
selected_data_matrix = list()
gene_attributes = list()
selected_gene_attributes = list()
snr_tuples = list()
attribute_selection_count = 100
ALL_count = 0
AML_count = 0

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')

    if iteration == total: 
        print()

def mean(values):
    return (sum(values) / len(values))

def variance(values):

    mean_value = mean(values)
    sum_mean_diff_sqr = 0

    for value in values:
        sum_mean_diff_sqr += math.pow((value - mean_value), 2)

    return sum_mean_diff_sqr / (len(values) - 1)

def standard_deviation(values):
    return math.sqrt(variance(values))

def mod_SNR(class_1_values, class_2_values):
    return abs((mean(class_1_values) + mean(class_2_values)) / (standard_deviation(class_1_values) + standard_deviation(class_2_values)))

def ZScore(value, value_list):
    return (value - mean(value_list)) / standard_deviation(value_list)

def swap_for_tuples(tuples, first_index, second_index):

    temp = tuples[second_index]
    tuples[second_index] = tuples[first_index]
    tuples[first_index] = temp

def partition_for_tuples(tuples, head, tail):

    pivot = tuples[tail][1]
    pivot_index = head
    for i in range(head, tail):
        if tuples[i][1] >= pivot:
            swap_for_tuples(tuples, pivot_index, i)
            pivot_index += 1

    swap_for_tuples(tuples, pivot_index, tail)

    return pivot_index

def randomized_partition_for_tuples(tuples, head, tail):

    swap_for_tuples(tuples, head, random.randint(head, tail))
    return partition_for_tuples(tuples, head, tail)

def randomized_quick_sort_for_tuples(tuples, head, tail):

    if(head < tail):
        pivot = randomized_partition_for_tuples(tuples, head, tail)
        randomized_quick_sort_for_tuples(tuples, head, pivot - 1)
        randomized_quick_sort_for_tuples(tuples, pivot + 1, tail)

def transfer_list(source, destination):

    for item in source:
        destination.append(item)

def take_user_input():

    global attribute_selection_count

    selection_count = input("Enter the number of gene attributes to select (default is " + str(attribute_selection_count) + "): ")

    if not selection_count:
        selection_count = "100"

    attribute_selection_count = int(selection_count)

    printProgressBar(0, 9, prefix = 'Progress:', suffix = 'Complete', length = 50)

def read_leukemia_raw_dataset():

    global raw_data_matrix

    with open("leukemia.txt", 'r') as datafile:
        for line in datafile:

            splitted_line_list = line.split("\t")
            raw_data_matrix.append(splitted_line_list)

    printProgressBar(1, 9, prefix = 'Progress:', suffix = 'Complete', length = 50)

def tidy_raw_dataset():

    global raw_data_matrix
    global data_matrix
    global gene_attributes

    gene_attributes = raw_data_matrix[0][1:]
    data_matrix = raw_data_matrix[3:-1]

    printProgressBar(2, 9, prefix = 'Progress:', suffix = 'Complete', length = 50)

def aggregate_same_class_samples():

    global data_matrix
    global ALL_count
    global AML_count

    ALL_list = list()
    AML_list = list()

    for sample in data_matrix:
        if sample[0] == "ALL":
            ALL_list.append(copy.deepcopy(sample))
        elif sample[0] == "AML":
            AML_list.append(copy.deepcopy(sample))

    ALL_count = len(ALL_list)
    AML_count = len(AML_list)
    data_matrix = list()

    transfer_list(ALL_list, data_matrix)
    transfer_list(AML_list, data_matrix)

    printProgressBar(3, 9, prefix = 'Progress:', suffix = 'Complete', length = 50)

def remove_newline_regex():

    global data_matrix
    global gene_attributes

    gene_attributes[-1] = gene_attributes[-1][:-1]

    for sample in data_matrix:
        sample[-1] = sample[-1][:-1]

    printProgressBar(4, 9, prefix = 'Progress:', suffix = 'Complete', length = 50)

def convert_datapoints_to_number():

    global data_matrix

    for sample in data_matrix:
        for index in range(1, len(sample)):
            sample[index] = float(sample[index])

    printProgressBar(5, 9, prefix = 'Progress:', suffix = 'Complete', length = 50)

def calculate_ZScore():

    global data_matrix

    sample_length = len(data_matrix[0])
    sample_count = len(data_matrix)

    for attribute_index in range(1, sample_length):

        attribute_list = list()
        normalized_attribute_list = list()

        for sample_index in range(sample_count):
            attribute_list.append(data_matrix[sample_index][attribute_index])

        for attribute in attribute_list:
            
            z_score = ZScore(attribute, attribute_list)
            rounded_zscore = math.ceil(z_score * 10000) / 10000
            normalized_attribute_list.append(rounded_zscore)

        for sample_index in range(sample_count):
            data_matrix[sample_index][attribute_index] = normalized_attribute_list[sample_index]

    printProgressBar(6, 9, prefix = 'Progress:', suffix = 'Complete', length = 50)

def sort_by_SNR():

    global data_matrix
    global ALL_count
    global snr_tuples

    sample_length = len(data_matrix[0])
    sample_count = len(data_matrix)

    for attribute_index in range(1, sample_length):

        attribute_list = list()

        for sample_index in range(sample_count):
            attribute_list.append(data_matrix[sample_index][attribute_index])

        snr_value = mod_SNR(attribute_list[:ALL_count], attribute_list[ALL_count:])
        rounded_snr = math.ceil(snr_value * 1000) / 1000
        snr_tuples.append((attribute_index, rounded_snr))

    randomized_quick_sort_for_tuples(snr_tuples, 0, len(snr_tuples) - 1)

    printProgressBar(7, 9, prefix = 'Progress:', suffix = 'Complete', length = 50)

def prepare_selected_dataset():

    global data_matrix
    global gene_attributes
    global selected_gene_attributes
    global snr_tuples
    global attribute_selection_count

    flag = 0

    for sample in data_matrix:

        selected_attribute_sample = list()

        for index in range(attribute_selection_count):

            snr_tuple = snr_tuples[index]
            selected_index = snr_tuple[0]

            if flag == 0:
                selected_gene_attributes.append(gene_attributes[selected_index - 1])

            selected_attribute_sample.append(sample[selected_index])

        selected_attribute_sample.append(sample[0])
        selected_data_matrix.append(copy.deepcopy(selected_attribute_sample))
        
        flag = 1

        printProgressBar(8, 9, prefix = 'Progress:', suffix = 'Complete', length = 50)

def write_to_file():

    global selected_data_matrix
    global selected_gene_attributes
    global attribute_selection_count


    filename = "leukemia-selected-" + str(attribute_selection_count) + ".csv"

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

    printProgressBar(9, 9, prefix = 'Progress:', suffix = 'Complete', length = 50)
    print("The pre-processed data has been saved in the file: ", filename, "!", end = "\n")

def main():

    take_user_input()
    read_leukemia_raw_dataset()
    tidy_raw_dataset()
    aggregate_same_class_samples()
    remove_newline_regex()
    convert_datapoints_to_number()
    calculate_ZScore()
    sort_by_SNR()
    prepare_selected_dataset()
    write_to_file()

    return

main()