//
//  main.cpp
//  Personalis_C
//
//  Created by Celine Becquet on 1/16/21.
//  Copyright © 2021 Celine Becquet. All rights reserved.
//
//https://github.com/pizzard/ska_sort
//https://probablydance.com/2016/12/27/i-wrote-a-faster-sorting-algorithm/

#include <iostream>
#include <functional>
#include <vector>
#include <random>
#include <algorithm>
#include <climits>
#include <iterator>
 
// A function to do counting sort of values according to
// the digit represented by exp.
void countSort( std::vector<unsigned short int> &values,  int exp)
{
    int vector_size = values.size();
    //int output[n]; // output array
    std::vector<unsigned short int> output(values.size());
    int i, count[10] = { 0 };

    // Store count of occurrences in count[]
    for (i = 0; i < vector_size; i++)
        count[(values[i] / exp) % 10]++;


    // Change count[i] so that count[i] now contains actual
    //  position of this digit in output[]
    for (i = 1; i < 10; i++)
        count[i] += count[i - 1];

    // Build the output array
    for (int i = vector_size - 1; i >= 0; i--) {
        output[count[(values[i] / exp) % 10] - 1] = values[i];
        count[(values[i] / exp) % 10]--;
    }

    // Copy the output vector to values, so that values now
    // contains sorted numbers according to current digit
    for (size_t i = 0; i < vector_size; i++)
        values[i] = output[i];
}

  
std::vector<unsigned short int> sort_short_ints(const std::vector<unsigned short int> &values)
{
    int vector_size = values.size();
    std::vector<unsigned short int> sorted_values(values.size());
    sorted_values = values;
    //TODO:  Implement a linear (O(n)) sort algorithm for unsigned 16-bit integers

    // Find the maximum number to know number of digits
    unsigned short  m = *max_element(values.begin(),values.end());
 
    // Do counting sort for every digit. Note that instead
    // of passing digit number, exp is passed. exp is 10^i
    // where i is current digit number
    for (int exp = 1; m / exp > 0; exp *= 10)
    {
       //countSort(sorted_values,  exp);
        std::vector<unsigned short int> sorted_tmp(values.size());
        int i, count[10] = { 0 };

        // Store count of occurrences in count[]
        for (i = 0; i < vector_size; i++)
            count[(sorted_values[i] / exp) % 10]++;

        // Change count[i] so that count[i] now contains actual
        //  position of this digit in output[]
        for (i = 1; i < 10; i++)
            count[i] += count[i - 1];

        // Build the output array
        for ( i = vector_size - 1; i >= 0; i--) {
            sorted_tmp[count[(sorted_values[i] / exp) % 10] - 1] = sorted_values[i];
            count[(sorted_values[i] / exp) % 10]--;
        }

        // Copy the output vector to values, so that values now
        // contains sorted numbers according to current digit
        for ( i = 0; i < vector_size; i++)
            sorted_values[i] = sorted_tmp[i];
    }
    return sorted_values;
}

int main(int argc, char *argv[])
{
    //Generate a vector of random unsigned short integers
    int seed = 12345;
    int vector_size = 1000 ; //1000000 ;
    std::vector<unsigned short int> random_values(vector_size);
    std::default_random_engine engine(seed);
    std::uniform_int_distribution<unsigned short int> dist(0, USHRT_MAX);
  
    for (size_t i=0; i<vector_size; ++i)
    {
        random_values[i] = dist(engine);
    }
//    random_values =  { 5, 6, 19, 2, 5, 7, 0, 23, 6, 256, 255, 8, 99, 1024, 65535, 65534 };
//    vector_size = random_values.size();
    
    //Uncomment to show random values
//    for (size_t i=0; i<vector_size; i++)
//    {
//        std::cerr << random_values[i] << ",";
//    }
//    std::cerr << std::endl;
    

    //Do not modify this line:
    std::vector<unsigned short int> sorted_values = sort_short_ints(random_values);

    //Uncomment to show sorted values
//    for (size_t i=0; i<sorted_values.size(); i++)
//       {
//           std::cerr << sorted_values[i] << ",";
//       }
//       std::cerr << std::endl;

    //Test
    if (sorted_values.size() != random_values.size())
    {
        std::cerr << "Error: sorted vector is not the same length as input vector" << std::endl;
        exit(1);
    }

    for (size_t i=1; i<sorted_values.size(); ++i)
    {
        if (sorted_values[i] < sorted_values[i-1])
        {
            std::cerr << "Error: sorted vector is not sorted:" << std::endl;

            std::cerr << "    [" << i-1 << "] = " << sorted_values[i-1] << std::endl;
            std::cerr << "    [" << i   << "] = " << sorted_values[i]   << std::endl;
            exit(1);
        }
    }

    std::cout << "Done sorting vector of "  << vector_size << " unsigned 16-bit integers " << std::endl;
    return 0;
}

