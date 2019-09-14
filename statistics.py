import math

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