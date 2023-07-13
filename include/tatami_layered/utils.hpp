#ifndef TATAMI_LAYERED_UTILS_HPP
#define TATAMI_LAYERED_UTILS_HPP

#include <limits>
#include <vector>
#include <numeric>

#include "tatami/tatami.hpp"

namespace tatami_layered {

enum class Category : unsigned char { U8, U16, U32 };

template<typename T>
Category categorize(T v) {
    constexpr uint8_t max8 = std::numeric_limits<uint8_t>::max();
    constexpr uint16_t max16 = std::numeric_limits<uint16_t>::max();
    constexpr T maxv = std::numeric_limits<T>::max();

    if constexpr(maxv <= max8) {
        return Category::U8;
    } else {
        if (v <= max8) {
            return Category::U8;
        }
    }

    if constexpr(maxv <= max16) {
        return Category::U16;
    } else {
        if (v <= max16) {
            return Category::U16;
        }
    }

    return Category::U32;
}

template<typename Int_, typename Index_, typename ColIndex_>
struct Holder {
    Holder() : ptr(1) {}
    std::vector<ColIndex_> index;
    std::vector<Int_> value;
    std::vector<size_t> ptr;

    void fill() {
        index.resize(ptr.back());
        value.resize(ptr.back());
    }
};

template<typename IndexIn_, typename ColIndex_> 
void allocate_rows(
    const std::vector<std::vector<Category> >& max_per_chunk,
    const std::vector<std::vector<IndexIn_> >& num_per_chunk,
    std::vector<std::vector<IndexIn_> >& identities8,
    std::vector<std::vector<IndexIn_> >& identities16,
    std::vector<std::vector<IndexIn_> >& identities32,
    std::vector<Holder<uint8_t, IndexIn_, ColIndex_> >& store8,
    std::vector<Holder<uint16_t, IndexIn_, ColIndex_> >& store16,
    std::vector<Holder<uint32_t, IndexIn_, ColIndex_> >& store32,
    std::vector<std::vector<Category> >& assigned_category,
    std::vector<std::vector<IndexIn_> >& assigned_position)
{
    size_t num_chunks = max_per_chunk.size();

    for (size_t chunk = 0; chunk < num_chunks; ++chunk) {
        IndexIn_ counter8 = 0, counter16 = 0, counter32 = 0;
        const auto& current_max = max_per_chunk[chunk];
        const auto& current_num = num_per_chunk[chunk];
        IndexIn_ NR = current_max.size();

        auto& asscat = assigned_category[chunk];
        asscat.resize(NR);
        auto& asspos = assigned_position[chunk];
        asspos.resize(NR);

        for (IndexIn_ r = 0; r < NR; ++r) {
            auto cat = current_max[r];
            auto num = current_num[r];
            IndexIn_ counter;

            switch(current_max[r]) {
                case Category::U8:
                    store8[chunk].ptr.push_back(store8[chunk].ptr.back() + num);
                    counter = counter8++;
                    identities8[chunk].push_back(r);
                    break;

                case Category::U16:
                    store16[chunk].ptr.push_back(store16[chunk].ptr.back() + num);
                    counter = counter16++;
                    identities16[chunk].push_back(r);
                    break;

                case Category::U32:
                    store32[chunk].ptr.push_back(store32[chunk].ptr.back() + num);
                    counter = counter32++;
                    identities32[chunk].push_back(r);
                    break;
            }

            asscat[r] = cat;
            asspos[r] = counter;
        }

        store8[chunk].fill();
        store16[chunk].fill();
        store32[chunk].fill();
    }
}

template<typename IndexIn_, typename ColIndex_> 
size_t get_sparse_ptr(
    const std::vector<Holder<uint8_t, IndexIn_, ColIndex_> >& store8,
    const std::vector<Holder<uint16_t, IndexIn_, ColIndex_> >& store16,
    const std::vector<Holder<uint32_t, IndexIn_, ColIndex_> >& store32,
    const std::vector<std::vector<Category> >& assigned_category,
    const std::vector<std::vector<IndexIn_> >& assigned_position,
    size_t chunk,
    IndexIn_ row)
{
    auto arow = assigned_position[chunk][row];
    switch (assigned_category[chunk][row]) {
        case Category::U8:
            return store8[chunk].ptr[arow];
        case Category::U16:
            return store16[chunk].ptr[arow];
        case Category::U32:
            return store32[chunk].ptr[arow];
    }
    return 0;
}

template<typename IndexIn_, typename ColIndex_, typename ValueIn_>
void fill_sparse_value(
    std::vector<Holder<uint8_t, IndexIn_, ColIndex_> >& store8,
    std::vector<Holder<uint16_t, IndexIn_, ColIndex_> >& store16,
    std::vector<Holder<uint32_t, IndexIn_, ColIndex_> >& store32,
    Category cat,
    size_t chunk, 
    IndexIn_ col, 
    ValueIn_ val,
    size_t output_position) 
{
    switch (cat) {
        case Category::U8:
            store8[chunk].value[output_position] = val;
            store8[chunk].index[output_position] = col;
            break;

        case Category::U16:
            store16[chunk].value[output_position] = val;
            store16[chunk].index[output_position] = col;
            break;

        case Category::U32:
            store32[chunk].value[output_position] = val;
            store32[chunk].index[output_position] = col;
            break;
    }
}

template<typename ValueOut_, typename IndexOut_, typename IndexIn_, typename ColIndex_> 
std::shared_ptr<tatami::Matrix<ValueOut_, IndexOut_> > consolidate_matrices(
    const std::vector<std::vector<IndexIn_> >& identities8,
    const std::vector<std::vector<IndexIn_> >& identities16,
    const std::vector<std::vector<IndexIn_> >& identities32,
    std::vector<Holder<uint8_t, IndexIn_, ColIndex_> > store8,
    std::vector<Holder<uint16_t, IndexIn_, ColIndex_> > store16,
    std::vector<Holder<uint32_t, IndexIn_, ColIndex_> > store32,
    IndexIn_ NR,
    IndexIn_ chunk_size,
    size_t leftovers)
{
    std::vector<std::shared_ptr<tatami::Matrix<ValueOut_, IndexOut_> > > col_combined;
    size_t num_chunks = identities8.size();
    col_combined.reserve(num_chunks);

    for (size_t c = 0; c < num_chunks; ++c) {
        std::vector<std::shared_ptr<tatami::Matrix<ValueOut_, IndexOut_> > > row_combined;
        IndexOut_ current_size = (c + 1 == num_chunks ? leftovers : chunk_size);
        row_combined.reserve(3);

        row_combined.emplace_back(new tatami::CompressedSparseRowMatrix<ValueOut_, IndexOut_, std::vector<uint8_t>, std::vector<ColIndex_> >(
            identities8[c].size(), current_size, std::move(store8[c].value), std::move(store8[c].index), std::move(store8[c].ptr), false
        ));
        row_combined.emplace_back(new tatami::CompressedSparseRowMatrix<ValueOut_, IndexOut_, std::vector<uint16_t>, std::vector<ColIndex_> >(
            identities16[c].size(), current_size, std::move(store16[c].value), std::move(store16[c].index), std::move(store16[c].ptr), false
        ));
        row_combined.emplace_back(new tatami::CompressedSparseRowMatrix<ValueOut_, IndexOut_, std::vector<uint32_t>, std::vector<ColIndex_> >(
            identities32[c].size(), current_size, std::move(store32[c].value), std::move(store32[c].index), std::move(store32[c].ptr), false
        ));

        std::vector<IndexIn_> reordered(NR);
        IndexIn_ counter = 0;
        for (auto& i : identities8[c]) {
            reordered[i] = counter++;
        }
        for (auto& i : identities16[c]) {
            reordered[i] = counter++;
        }
        for (auto& i : identities32[c]) {
            reordered[i] = counter++;
        }

        col_combined.push_back(tatami::make_DelayedSubset<0>(tatami::make_DelayedBind<0>(std::move(row_combined)), std::move(reordered)));
    }

    return tatami::make_DelayedBind<1>(std::move(col_combined));
}

}

#endif
