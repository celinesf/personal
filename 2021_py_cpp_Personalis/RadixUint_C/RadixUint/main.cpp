//
//  main.cpp
//  RadixUint
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
//
//  radix1.cpp
//  Personalis_C
//
//  Created by Celine Becquet on 1/22/21.
//  Copyright © 2021 Celine Becquet. All rights reserved.
//


// C++ implementation of Radix Sort
#include <iostream>
using namespace std;

// A utility function to print an array
void print(const std::vector<unsigned short int> &arr, size_t n)
{
    for (size_t  i = 0; i < n ; i++)
        cout << arr[i] << " ";
    cout  << endl;
}

void printint(int arr[], size_t n)
{
    for (size_t  i = 0; i < n ; i++)
        cout << arr[i] << " ";
    cout  << endl;
}

// A utility function to get maximum value in arr[]
int getMax(int arr[], size_t n)
{
    int mx = arr[0];
    for (int i = 1; i < n; i++)
        if (arr[i] > mx)
            mx = arr[i];
    return mx;
}
  
// A function to do counting sort of arr[] according to
void countSort( std::vector<unsigned short int> &arr, size_t n, int exp)
{
    //int output[n]; // output array
    std::vector<unsigned short int> output(arr.size());
    int i, count[10] = { 0 };

    // Store count of occurrences in count[]
    for (i = 0; i < n; i++)
        count[(arr[i] / exp) % 10]++;


    // Change count[i] so that count[i] now contains actual
    //  position of this digit in output[]
    for (i = 1; i < 10; i++)
        count[i] += count[i - 1];

    // Build the output array
    for (int i = int(n) - 1; i >= 0; i--) {
        output[count[(arr[i] / exp) % 10] - 1] = arr[i];
        count[(arr[i] / exp) % 10]--;
         //cout << "i " << i <<  " out " << output[count[(arr[i] / exp) % 10] - 1] << " count " << count[(arr[i] / exp) % 10] << " index " << count[(arr[i] / exp) % 10] <<endl;
    }
//    cout << "outut " << endl;
//    print(output,n);

    // Copy the output array to arr[], so that arr[] now
    // contains sorted numbers according to current digit
    for (size_t i = 0; i < n; i++)
        arr[i] = output[i];
//  cout << "arr " << endl;
//    print(arr,n);
}



  
// The main function to that sorts arr[] of size n using
// Radix Sort
void radixsort( std::vector<unsigned short int> &arr, size_t n)
{
    // Find the maximum number to know number of digits
    //int m = getMax(arr, n);
    unsigned short  m = *max_element(arr.begin(),arr.end());
 
    cout << "HERE " << m  << endl;
        // Do counting sort for every digit. Note that instead
    // of passing digit number, exp is passed. exp is 10^i
    // where i is current digit number
    for (int exp = 1; m / exp > 0; exp *= 10)
    {
        countSort(arr, n, exp);
        cout << "HERE2 "<< exp << " exp "<<  m / exp << endl;
    }
}
  

  
// Driver Code
int main()
{
    int vector_size = 1000000 ;

    std::vector<unsigned short int> arr(vector_size);
    // randome value of it
    int seed = 12345;
    std::default_random_engine engine(seed);
    std::uniform_int_distribution<uint> dist(0, USHRT_MAX);

    for (size_t i=0; i< vector_size; ++i)
    {
      arr[i] = dist(engine);
    }
    size_t n = arr.size();
     cout  << n<<endl;
     
    //print(arr, n);
    
    
    // Function Call
    radixsort(arr, n);
    //print(arr, n);
    


    //Test
    if (arr.size() != vector_size)
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

