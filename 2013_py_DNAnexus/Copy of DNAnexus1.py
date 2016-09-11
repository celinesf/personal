#!/usr/bin/python

'''
Created on Oct 11, 2013 2:36 pm

@author: Celine Becquet
'''

def getDNA(bit_pair):
    if bit_pair == '00':
        return 'A'
    if bit_pair == '01':
        return 'C'
    if bit_pair == '10':
        return 'G'
    if bit_pair == '11':
        return 'T' 
    else:
        print'I got a problem with getDNA', bit_pair
        return None



def problem1(L,filename):
    input_file = open(filename,'rb') 
    out_file = open(filename+str(L),'w')
    c= input_file.read(1)
    n_dna = 1 ## counter on number of sequences
    n_c = 0 # counter on characters considered for this sequence/ pieace of DNA
    piece = {"@READ_":"","+READ_":""}
    while c:
        ### get bit in a byte
        byte = ord(c)
        byte = bin(byte)[2:].rjust(8, '0')
        ### convert bits to data
        letter = getDNA(byte[0:2])
        score = chr(33+int(byte[2:8], 2))
        
        ### append letter for this piece of DNA
        if n_c < L:
            piece['@READ_'] += letter
            piece['+READ_'] += score
        n_c +=1
        
        ### output full piece of DNA
        if n_c == L:
            out_file.write('@READ_'+str(n_dna)+'\n'+piece['@READ_']+"\n"+'+READ_'+str(n_dna)+'\n'+piece['+READ_']+"\n")
            n_c = 0 
            piece = {"@READ_":"","+READ_":""}
            n_dna += 1
        

        ### read next character
        c = input_file.read(1)

    input_file.close()
    out_file.close() 



if __name__ == '__main__':
    
    problem1(80,'input')
      
    pass