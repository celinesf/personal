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



struct PartitionInfo
{
    PartitionInfo()
        : count(0)
    {
    }
 
    union
    {
        size_t count;
        size_t offset;
    };
    size_t next_offset;
};

template<typename It, typename F>
inline It custom_std_partition(It begin, It end, F && func)
{
    for (;; ++begin)
    {
        if (begin == end)
            return end;
        if (!func(*begin))
            break;
    }
    It it = begin;
    for(++it; it != end; ++it)
    {
        if (!func(*it))
            continue;

        std::iter_swap(begin, it);
        ++begin;
    }
    return begin;
}

std::vector<unsigned short int> sort_short_ints(const std::vector<unsigned short int> &values)
{
    std::vector<unsigned short int> sorted_values(values.size());
  
    //TODO:  Implement a linear (O(n)) sort algorithm for unsigned 16-bit integers
    auto begin = values.begin();
    auto end = values.end();
   
    
//    unsigned short min=65535;
//    unsigned short max = 0;
//
//    size_t counts[ 65535] = {};
//    for (auto it = begin; it != end; ++it)
//    {
//        ++counts[*it];
//    }
//   //std::cout << ' ' << min << ' ' << max << std::endl;
//    size_t total = 0;
//    for (size_t & count : counts)
//    {
//        //std::cout << ' ' << count;
//        size_t old_count = count;
//        count = total;
//        total += old_count;
//    }
//    //std::cout << ' ' << total << std::endl;
//    for (; begin != end; ++begin)
//    {
//        auto key = *begin;
//        //std::cout << ' ' << *begin << ' ' << total <<  ' ' << counts[key] << std::endl ;
//        sorted_values[counts[key]++] = std::move(*begin);
//    }
//    //std::cout << std::endl;
    
    
    PartitionInfo partitions[65535];
    for (auto it = begin; it != end; ++it)
    {
       ++partitions[*it].count;
    }
    uint16_t remaining_partitions[65535];
    size_t total = 0;
    int num_partitions = 0;
    for (int i = 0; i < 65535; ++i)
    {
       size_t count = partitions[i].count;
       if (count)
       {
           partitions[i].offset = total;
           total += count;
           remaining_partitions[num_partitions] = i;
           ++num_partitions;
       }
       partitions[i].next_offset = total;
    }
    for (uint16_t * last_remaining = remaining_partitions + num_partitions, * end_partition = remaining_partitions + 1; last_remaining > end_partition;)
    {
       last_remaining = custom_std_partition(remaining_partitions, last_remaining, [&](uint16_t partition)
       {
           size_t & begin_offset = partitions[partition].offset;
           size_t & end_offset = partitions[partition].next_offset;
           if (begin_offset == end_offset)
               return false;

           unroll_loop_four_times(begin + begin_offset, end_offset - begin_offset, [partitions = partitions, begin, &extract_key, sort_data](It it)
           {
               uint16_t this_partition = extract_key(*it);
               size_t offset = partitions[this_partition].offset++;
               std::iter_swap(it, begin + offset);
           });
           return begin_offset != end_offset;
       });
    }
    
    return sorted_values;
}

int main(int argc, char *argv[])
{
    //Generate a vector of random unsigned short integers
    int seed = 12345;
    int vector_size = 8 ; //1000000 ;
    std::vector<unsigned short int> random_values(vector_size);
    std::default_random_engine engine(seed);
    std::uniform_int_distribution<unsigned short int> dist(0, USHRT_MAX);
  
    for (size_t i=0; i<vector_size; ++i)
    {
        random_values[i] = dist(engine);
    }
random_values = { 4, 4, 2, 4, 1, 1, 4, 5, 4};

    //Uncomment to show random values
    for (size_t i=0; i<vector_size; i++)
    {
        std::cerr << random_values[i] << ",";
    }
    std::cerr << std::endl;
    

    //Do not modify this line:
    std::vector<unsigned short int> sorted_values = sort_short_ints(random_values);
   
    //Uncomment to show sorted values
    for (size_t i=0; i<sorted_values.size(); i++)
       {
           std::cerr << sorted_values[i] << ",";
       }
       std::cerr << std::endl;

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

