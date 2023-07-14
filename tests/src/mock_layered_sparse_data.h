#ifndef MOCK_LAYERED_SPARSE_DATA_H
#define MOCK_LAYERED_SPARSE_DATA_H

#include <random>
#include <vector>

template<bool ByRow>
void mock_layered_sparse_data(size_t NR, size_t NC, std::vector<size_t>& rows, std::vector<size_t>& cols, std::vector<int>& vals) {
    std::mt19937_64 rng(ByRow ? NR : NC);

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

    auto secondary_dim = (ByRow ? NC : NR);
    std::vector<unsigned char> present(secondary_dim);
    auto choose = [&](auto& primary_store, auto& secondary_store, int primary_index) -> void {
        // This averages out to around about two values per row. This gives
        // us something interesting w.r.t. multiple values at different limits
        // (e.g., about 1/9th of the rows will have both values at the lowest limit).
        size_t number = (rng() % 3) + 1; 

        size_t original = secondary_store.size();
        for (size_t i = 0; i < number; ++i) {
            size_t candidate;
            do {
                candidate = rng() % secondary_dim;
            } while (present[candidate]);
            present[candidate] = 1;

            secondary_store.push_back(candidate);
            primary_store.push_back(primary_index);
            vals.push_back(rand());
        }

        std::sort(secondary_store.begin() + original, secondary_store.end());
        for (size_t i = 0; i < number; ++i) {
            present[secondary_store[original + i]] = 0;
        } 
        return;
    };

    if constexpr(!ByRow) {
        for (size_t c = 0; c < NC; ++c) {
            choose(cols, rows, c);
        }
    } else {
        for (size_t r = 0; r < NR; ++r) {
            choose(rows, cols, r);
        }
    }

    return;
}

#endif
