import statistics

def calculate_zscore(value, value_list):
    return (value - statistics.mean(value_list)) / statistics.standard_deviation(value_list)