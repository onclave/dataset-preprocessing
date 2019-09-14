import random

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