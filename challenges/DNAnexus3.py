#!/usr/bin/python

'''
Created on Oct 11, 2013 2:36 pm

@author: Celine Becquet
'''


import sys
import os.path

def problem3(position, filename):
    if os.path.exists(filename + '.idx'):
        idx = open(filename + '.idx', 'r')
        index = idx.read().split(' ')
        if int(position)+1 >= len(index):
            print 'The file stop at line ' + str(len(index)-1)
        else:
            offset = index[int(position)]
            input_file = open(filename , 'r')
            input_file.seek(int(offset))
            print  input_file.readline()
            input_file.close()
        idx.close()
       
    else:
        print "Writing index to " + filename+ ".idx..." ,
        idx = open(filename + '.idx', 'wb')
        input_file =open(filename , 'r')
        c= input_file.read(1)
        i = 0 # position of character
        idx.write(str(i)+' ')
        while c:
            if '\n' in c:
                idx.write(str(i+1)+' ')
            c= input_file.read(1)
            i+=1
        idx.close()
        input_file.close() 
        print 'done.'
        problem3(position, filename)
              

if __name__ == '__main__':
    position = sys.argv[1]
    filename = sys.argv[2]

    problem3(position, filename)
      
    pass