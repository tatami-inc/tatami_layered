#ifndef MOCK_LAYERED_SPARSE_DATA_H
#define MOCK_LAYERED_SPARSE_DATA_H

#include <random>
#include <vector>

inline void mock_layered_sparse_data(std::size_t NR, std::size_t NC, std::vector<std::size_t>& rows, std::vector<std::size_t>& cols, std::vector<int>& vals) {
    std::mt19937_64 rng(NR * NC + 1);

    // This idea is to always have some values in all categories;
    // we wouldn't get this if we uniformly sampled.
    auto rand = [&]() -> int {
        int upper = 1;
        switch (rng() % 3) {
            case 0:
                upper = 256;
                break;
            case 1:
                upper = 65536;
                break;
            case 2:
                upper = 10000000; // that's big enough.
                break;
        }
        return rng() % upper;
    };

    std::vector<unsigned char> present(NC);
    for (decltype(NR) r = 0; r < NR; ++r) {
        auto original = cols.size();

        // This averages out to around about two values per row. This gives
        // us something interesting w.r.t. multiple values at different limits
        // (e.g., about 1/9th of the rows will have both values at the lowest limit).
        std::size_t number = (rng() % 3) + 1; 

        for (decltype(number) i = 0; i < number; ++i) {
            std::size_t candidate;
            do {
                candidate = rng() % NC;
            } while (present[candidate]);
            present[candidate] = 1;

            cols.push_back(candidate);
            rows.push_back(r);
            vals.push_back(rand());
        }

        std::sort(cols.begin() + original, cols.end());
        for (size_t i = 0; i < number; ++i) {
            present[cols[original + i]] = 0;
        } 
    }

    return;
}

#endif
