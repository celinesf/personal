#!/usr/bin/python

'''
Created on Oct 11, 2013 2:36 pm

@author: Celine Becquet
'''


def problem2(x, input_array):

    j = 0
    while input_array[j] < x :
        j +=1
        if j == len(input_array): break
    if j < len(input_array):
        return j
    else: return -1


if __name__ == '__main__':
    
    print problem2(10,[1,2,2,3,3,5,5,6,7,9])
    
    
    pass