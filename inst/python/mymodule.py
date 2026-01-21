import numpy as np

def add_two_numbers(a, b):
    return a + b

def sum_with_numpy(values):
    """
    Take a 1D list/array of numbers, convert to NumPy, and return the sum.
    """
    arr = np.array(values)
    return float(arr.sum())
