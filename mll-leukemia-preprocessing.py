import copy
import math
import progressbar
import statistics
import snr
import zscore
import sort

raw_data_matrix = list()
transposed_raw_data_matrix = list()
selected_data_matrix = list()
gene_attributes = list()
ALL_samples = list()
AML_samples = list()
MLL_samples = list()
ALL_count = 0
AML_count = 0
MLL_count = 0
attribute_selection_count = 100
selected_gene_attributes = list()
snr_tuples = list()
progressbar_total = 11


def take_user_input():

    global attribute_selection_count
    global progressbar_total

    selection_count = input(
        "Enter the number of gene attributes to select (default is "
        + str(attribute_selection_count)
        + "): "
    )

    if not selection_count:
        selection_count = "100"

    attribute_selection_count = int(selection_count)

    progressbar.show(
        0, progressbar_total, prefix="Progress:", suffix="Complete", length=50
    )


def read_mll_leukemia_raw_dataset():

    global raw_data_matrix
    global progressbar_total

    with open("mll_leukemia.txt", "r") as dataset:
        for line in dataset:
            raw_data_matrix.append(line.rstrip())

    progressbar.show(
        1, progressbar_total, prefix="Progress:", suffix="Complete", length=50
    )


def count_class_strength():

    global raw_data_matrix
    global ALL_count
    global AML_count
    global MLL_count
    global progressbar_total

    splitted_attributes = raw_data_matrix[0].split("\t")[2:]

    for word in splitted_attributes:
        if "ALL" in word:
            ALL_count = ALL_count + 1
        elif "AML" in word:
            AML_count = AML_count + 1
        elif "MLL" in word:
            MLL_count = MLL_count + 1

    progressbar.show(
        2, progressbar_total, prefix="Progress:", suffix="Complete", length=50
    )


def arrange_dataset():

    global raw_data_matrix
    global gene_attributes
    global ALL_samples
    global MLL_samples
    global AML_samples
    global ALL_count
    global AML_count
    global MLL_count
    global progressbar_total

    linecount = 0

    for line in raw_data_matrix[1:]:

        splitted_line_list = line.split("\t")
        gene_attributes.append(splitted_line_list[0])

        ALL_start_index = 2
        ALL_end_index = ALL_count + 2
        MLL_start_index = ALL_end_index
        MLL_end_index = ALL_count + 2 + MLL_count
        AML_start_index = MLL_end_index
        AML_end_index = ALL_count + MLL_count + AML_count + 2

        for index in range(len(splitted_line_list[ALL_start_index:ALL_end_index])):
            if linecount == 0:
                ALL_samples.append(list())
            ALL_samples[index].append(splitted_line_list[index + ALL_start_index])

        for index in range(len(splitted_line_list[MLL_start_index:MLL_end_index])):
            if linecount == 0:
                MLL_samples.append(list())
            MLL_samples[index].append(splitted_line_list[index + MLL_start_index])

        for index in range(len(splitted_line_list[AML_start_index:AML_end_index])):
            if linecount == 0:
                AML_samples.append(list())
            AML_samples[index].append(splitted_line_list[index + AML_start_index])

        linecount = linecount + 1

    progressbar.show(
        3, progressbar_total, prefix="Progress:", suffix="Complete", length=50
    )


def add_class_labels():

    global ALL_samples
    global MLL_samples
    global AML_samples
    global progressbar_total

    for sample in ALL_samples:
        sample.append("ALL")

    for sample in MLL_samples:
        sample.append("MLL")

    for sample in AML_samples:
        sample.append("AML")

    progressbar.show(
        4, progressbar_total, prefix="Progress:", suffix="Complete", length=50
    )


def transpose_dataset():

    global raw_data_matrix
    global gene_attributes
    global transposed_raw_data_matrix
    global ALL_samples
    global MLL_samples
    global AML_samples
    global progressbar_total

    for sample in ALL_samples:
        transposed_raw_data_matrix.append(sample)

    for sample in MLL_samples:
        transposed_raw_data_matrix.append(sample)

    for sample in AML_samples:
        transposed_raw_data_matrix.append(sample)

    progressbar.show(
        5, progressbar_total, prefix="Progress:", suffix="Complete", length=50
    )


def write_as_csv():

    global transposed_raw_data_matrix
    global gene_attributes
    global progressbar_total

    filename = "mll-leukemia.csv"

    writefile = open(filename, "w+")
    write_file_content = ""

    for attribute in gene_attributes:
        write_file_content += attribute + ","

    write_file_content += "class\n"

    for sample in transposed_raw_data_matrix:

        line = ""

        for value in sample[:-1]:
            line += str(value) + ","

        line += sample[-1] + "\n"
        write_file_content += line

    writefile.write(write_file_content)
    writefile.close()

    progressbar.show(
        6, progressbar_total, prefix="Progress:", suffix="Complete", length=50
    )


def convert_datapoints_to_number():

    global transposed_raw_data_matrix
    global progressbar_total

    for sample in transposed_raw_data_matrix:
        for index in range(len(sample) - 1):
            sample[index] = float(sample[index])

    progressbar.show(
        7, progressbar_total, prefix="Progress:", suffix="Complete", length=50
    )


def normalize_data():

    global transposed_raw_data_matrix
    global progressbar_total

    sample_length = len(transposed_raw_data_matrix[0])
    sample_count = len(transposed_raw_data_matrix)

    for attribute_index in range(sample_length - 1):

        attribute_list = list()
        normalized_attribute_list = list()

        for sample_index in range(sample_count):
            attribute_list.append(
                transposed_raw_data_matrix[sample_index][attribute_index]
            )

        for attribute in attribute_list:

            z_score = zscore.calculate_zscore(attribute, attribute_list)
            rounded_zscore = math.ceil(z_score * 10000) / 10000
            normalized_attribute_list.append(rounded_zscore)

        for sample_index in range(sample_count):
            transposed_raw_data_matrix[sample_index][attribute_index] = normalized_attribute_list[sample_index]

    progressbar.show(
        8, progressbar_total, prefix="Progress:", suffix="Complete", length=50
    )


def sort_by_standard_deviation():

    global transposed_raw_data_matrix
    global ALL_count
    global MLL_count
    global AML_count
    global snr_tuples
    global progressbar_total

    sample_length = len(transposed_raw_data_matrix[0])
    sample_count = len(transposed_raw_data_matrix)
    MLL_end_index = ALL_count + MLL_count

    for attribute_index in range(sample_length - 1):

        attribute_list = list()

        for sample_index in range(sample_count):
            attribute_list.append(transposed_raw_data_matrix[sample_index][attribute_index])

        sd_value = statistics.standard_deviation(attribute_list[:ALL_count]) + statistics.standard_deviation(attribute_list[ALL_count: MLL_end_index]) + statistics.standard_deviation(attribute_list[MLL_end_index:])
        rounded_sd = math.ceil(sd_value * 10000) / 10000
        snr_tuples.append((attribute_index, rounded_sd))

    sort.randomized_quick_sort_for_tuples(snr_tuples, 0, len(snr_tuples) - 1)

    progressbar.show(
        9, progressbar_total, prefix="Progress:", suffix="Complete", length=50
    )

def prepare_selected_dataset():

    global transposed_raw_data_matrix
    global selected_data_matrix
    global gene_attributes
    global selected_gene_attributes
    global snr_tuples
    global attribute_selection_count
    global progressbar_total

    flag = 0

    for sample in transposed_raw_data_matrix:

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

    progressbar.show(
        10, progressbar_total, prefix = 'Progress:', suffix = 'Complete', length = 50
    )

def write_to_file():

    global selected_data_matrix
    global selected_gene_attributes
    global attribute_selection_count
    global progressbar_total

    filename = "mll-leukemia-selected-" + str(attribute_selection_count) + ".csv"

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

    progressbar.show(11, progressbar_total, prefix = 'Progress:', suffix = 'Complete', length = 50)
    print("The pre-processed data has been saved in the file: ", filename, "!", end = "\n")


def main():

    take_user_input()
    read_mll_leukemia_raw_dataset()
    count_class_strength()
    arrange_dataset()
    add_class_labels()
    transpose_dataset()
    write_as_csv()
    convert_datapoints_to_number()
    normalize_data()
    sort_by_standard_deviation()
    prepare_selected_dataset()
    write_to_file()

    return


main()
