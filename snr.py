import math
import statistics

def mod_SNR(class_1_values, class_2_values):
    return abs((statistics.mean(class_1_values) + statistics.mean(class_2_values)) / (statistics.standard_deviation(class_1_values) + statistics.standard_deviation(class_2_values)))

def mod_SNR_3(class_1_values, class_2_values, class_3_values):

	return abs((statistics.mean(class_1_values) + statistics.mean(class_2_values) + statistics.mean(class_3_values)) / (statistics.standard_deviation(class_1_values) + statistics.standard_deviation(class_2_values) + statistics.standard_deviation(class_3_values)))