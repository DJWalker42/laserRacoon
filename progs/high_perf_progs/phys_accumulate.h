#ifndef LASER_RACOON_PHYS_ACCUMULATE_H
#define LASER_RACOON_PHYS_ACCUMULATE_H

#include <vector>
#include <future>
#include <thread>
#include <numeric>


/*
 *  Taken from stack overflow answer to:
 *
 *  https://stackoverflow.com/questions/28048539/calculating-the-sum-of-a-large-vector-in-parallel
 *
 *  with minor alterations for C++11
 *
 *  The part_size computation allows us to distribute the sizes to the threads as
 *  evenly as possible when the size of the input vector (container) is not an
 *  integer multiple of the number of logical threads (integer division truncation).
 *
 */

//Minimum number of elements for multithreaded algorithm.
#define MT_MIN_SIZE 10000

namespace phys {


template <typename T, typename U>
U omp_accumulate(const std::vector<T>& v, U init) {
    U sum = init;

    #pragma omp parallel for reduction(+:sum)
    for(std::size_t i = 0; i < v.size(); i++) {
        sum += v[i];
    }

    return sum;
}


template <typename InputIt, typename T>
T accumulate_p(InputIt first, InputIt last, T init) {
    // Determine total size.
    const auto size = std::distance(first, last);
    // Determine how many parts the work shall be split into.
    const auto parts = (size < MT_MIN_SIZE)? 1 : std::thread::hardware_concurrency();

    std::vector<std::future<T>> futures;

    // For each part, calculate size and run accumulate on a separate thread.
    for (std::size_t i = 0; i != parts; ++i) {
        const auto part_size = (size * i + size) / parts - (size * i) / parts;
        futures.emplace_back(std::async(std::launch::async,
            [=] { return std::accumulate(first, std::next(first, part_size), T{}); }));
        std::advance(first, part_size);
    }

    // Wait for all threads to finish execution and accumulate results.
    return std::accumulate(std::begin(futures), std::end(futures), init,
        [] (const T prev, std::future<T>& future) { return prev + future.get(); });
}

} //namespace

#endif //header guard
