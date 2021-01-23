'''
Create a list containing 10,000 random integers.  Write a function that
returns the position(s) in this list which contain the third-largest value
'''

"""
@author: Celine Becquet
@creation_date: 01/16/2021
"""

__author__ = "Celine Becquet"
__email__ = "celine.becquet@gmail.com"
__status__ = "dev"
__version__ = 1.0

import os
import doctest
import numpy as np



# The function that identifies the *position(s)* of the 3rd-largest value
def index_of_3rd_largest(numbers_list):
    """ Find index(es) of 3rd largest value in numbers_list
    Args:
        numbers_list (ndarray(dtype=int, ndim=1): List of random integers
    Example: Test index_of_3rd_largest with 10 random integer between 0-4 with set random seed
        >>> np.random.seed(seed=100)
        >>> l = np.random.randint(low=0, high=5, size=10)
        >>> print(l)
        [0 0 3 0 2 4 2 2 2 2]
        >>> index_of_3rd_largest(l)
        Index(es) of 3rd largest value (2):
        [4 6 7 8 9]
    """
    final_list = np.asarray(np.where(numbers_list == np.unique(numbers_list)[-3]))[0]
    print("Index(es) of 3rd largest value (" +str(np.unique(numbers_list)[-3]) +"):")
    print(final_list)
    return None

# Main function: initializes the list of 10000 random integers,
# calls the index_of_3rd_largest function, and prints the result to STDOUT
def main():
    # Specify seed for random generator
    seed = input("Enter a seed number to set random generator (press enter for no specific seed): ")

    if not seed.isdigit() :
        print("No seed specified for random generator\n")
    else :
        seed=int(seed)
        print("I will use seed '"+ str(seed)+"' to generate my list of 10K values\n")
        np.random.seed(seed=int(seed))

    # Specify maximum value within the list of integer
    int_range = input("Enter a number to define the highest number in the list (press enter for default of 10k): ")
    if not int_range.isdigit():
        int_range = 10000
    else:
        int_range = int(int_range)
    print("I will generate 10K values between 0 - " + str(int_range) + "\n")

    # initializes the list of 10000 random integers
    l = np.random.randint(low=0, high=int_range, size=10000)

    index_of_3rd_largest(l)
    pass

if __name__ == "__main__":
    main()
    doctest.testmod(verbose=True)

