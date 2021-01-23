//
//  main.cpp
//  RadixInt
//
//  Created by Celine Becquet on 1/22/21.
//  Copyright © 2021 Celine Becquet. All rights reserved.
//


#include <iostream>
#include <functional>
#include <vector>
#include <random>
#include <algorithm>
#include <climits>
#include <iterator>


// C++ implementation of Radix Sort
#include <iostream>
using namespace std;
  
// A utility function to get maximum value in arr[]
int getMax(int arr[], int n)
{
    int mx = arr[0];
    for (int i = 1; i < n; i++)
        if (arr[i] > mx)
            mx = arr[i];
    return mx;
}
  
// A function to do counting sort of arr[] according to
// the digit represented by exp.
void countSort(int arr[], int n, int exp)
{
    int output[n]; // output array
    int i, count[10] = { 0 };
  
    // Store count of occurrences in count[]
    for (i = 0; i < n; i++)
        count[(arr[i] / exp) % 10]++;
  
    // Change count[i] so that count[i] now contains actual
    //  position of this digit in output[]
    for (i = 1; i < 10; i++)
        count[i] += count[i - 1];
  
    // Build the output array
    for (i = n - 1; i >= 0; i--) {
        output[count[(arr[i] / exp) % 10] - 1] = arr[i];
        count[(arr[i] / exp) % 10]--;
    }
  
    // Copy the output array to arr[], so that arr[] now
    // contains sorted numbers according to current digit
    for (i = 0; i < n; i++)
        arr[i] = output[i];
}



  
// The main function to that sorts arr[] of size n using
// Radix Sort
void radixsort(int arr[], int n)
//void radixsort(int arr[], int n)
{
    // Find the maximum number to know number of digits
    int m = getMax(arr, n);
  
    // Do counting sort for every digit. Note that instead
    // of passing digit number, exp is passed. exp is 10^i
    // where i is current digit number
    for (int exp = 1; m / exp > 0; exp *= 10)
        countSort(arr, n, exp);
}
  
// A utility function to print an array
void print(int arr[], int n)
{
    for (int i = 0; i < n; i++)
        cout << arr[i] << " ";
    cout  << endl;
}
  
// Driver Code
int main()
{
    int vector_size = 1000000 ;
    int arr[vector_size]; // = { 170, 45, 75, 90, 802, 24, 2, 66 };
    // randome value of it
    int seed = 12345;

    std::default_random_engine engine(seed);
    std::uniform_int_distribution<int> dist(0, INT_MAX);
   
   
    for (size_t i=0; i< vector_size; ++i)
    {
      arr[i] = dist(engine);
    }
    int n = sizeof(arr) / sizeof(arr[0]);
     cout  << n<<endl;
     cout  << sizeof(arr)/sizeof(*arr)<<endl;
    //print(arr, n);
    
      // Function Call
    radixsort(arr, n);
    //print(arr, n);
    


    //Test
    if (sizeof(arr)/sizeof(*arr) != vector_size)
    {
        std::cerr << "Error: sorted vector is not the same length as input vector" << std::endl;
        exit(1);
    }

    for (size_t i=1; i< n; ++i)
    {
        if (arr[i] < arr[i-1])
        {
            std::cerr << "Error: sorted vector is not sorted:" << std::endl;

            std::cerr << "    [" << i-1 << "] = " << arr[i-1] << std::endl;
            std::cerr << "    [" << i   << "] = " << arr[i]   << std::endl;
            exit(1);
        }
    }

    std::cout << "Done sorting vector of "  << vector_size << " unsigned 16-bit integers " << std::endl;
    
    return 0;
}

