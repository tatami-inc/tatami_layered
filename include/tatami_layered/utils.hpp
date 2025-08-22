#ifndef TATAMI_LAYERED_UTILS_HPP
#define TATAMI_LAYERED_UTILS_HPP

#include <limits>
#include <vector>
#include <numeric>
#include <cstdint>
#include <cstddef>
#include <type_traits>

#include "tatami/tatami.hpp"
#include "sanisizer/sanisizer.hpp"

namespace tatami_layered {

enum class Category : unsigned char { U8, U16, U32 };

template<typename Value_>
Category categorize(const Value_ v) {
    if (v < 0) {
        throw std::runtime_error("cannot categorize negative value in a layered matrix");
    }

    constexpr std::uint32_t max32 = std::numeric_limits<std::uint32_t>::max();
    if constexpr(std::is_floating_point<Value_>::value) {
        if (v > static_cast<double>(max32)) {
            throw std::runtime_error("floating-point value is outside of the range of a 32-bit unsigned integer");
        }
        return categorize<std::uint32_t>(v);

    } else {
        constexpr std::uint8_t max8 = std::numeric_limits<std::uint8_t>::max();
        if (sanisizer::is_less_than_or_equal(v, max8)) {
            return Category::U8;
        }

        constexpr std::uint16_t max16 = std::numeric_limits<std::uint16_t>::max();
        if (sanisizer::is_less_than_or_equal(v, max16)) {
            return Category::U16;
        }

        if (sanisizer::is_greater_than(v, max32)) {
            throw std::runtime_error("integer value is outside of the range of a 32-bit unsigned integer");
        }

        return Category::U32;
    }
}

template<typename Int_, typename Index_, typename ColIndex_>
struct Holder {
    Holder() : ptr(1) {}
    std::vector<ColIndex_> index;
    std::vector<Int_> value;
    std::vector<std::size_t> ptr;

    void fill() {
        sanisizer::resize(index, ptr.back());
        sanisizer::resize(value, ptr.back());
    }
};

template<typename Input_>
std::remove_cv_t<std::remove_reference_t<Input_> > I(Input_ x) {
    return x;
}

template<typename IndexIn_, typename ColIndex_> 
void allocate_rows(
    const std::vector<std::vector<Category> >& max_per_chunk,
    const std::vector<std::vector<IndexIn_> >& num_per_chunk,
    std::vector<std::vector<IndexIn_> >& identities8,
    std::vector<std::vector<IndexIn_> >& identities16,
    std::vector<std::vector<IndexIn_> >& identities32,
    std::vector<Holder<std::uint8_t, IndexIn_, ColIndex_> >& store8,
    std::vector<Holder<std::uint16_t, IndexIn_, ColIndex_> >& store16,
    std::vector<Holder<std::uint32_t, IndexIn_, ColIndex_> >& store32,
    std::vector<std::vector<Category> >& assigned_category,
    std::vector<std::vector<IndexIn_> >& assigned_position)
{
    const IndexIn_ num_chunks = max_per_chunk.size();
    for (decltype(I(num_chunks)) chunk = 0; chunk < num_chunks; ++chunk) {
        IndexIn_ counter8 = 0, counter16 = 0, counter32 = 0;
        const auto& current_max = max_per_chunk[chunk];
        const auto& current_num = num_per_chunk[chunk];
        const IndexIn_ NR = current_max.size();

        auto& asscat = assigned_category[chunk];
        tatami::resize_container_to_Index_size(asscat, NR);
        auto& asspos = assigned_position[chunk];
        tatami::resize_container_to_Index_size(asspos, NR);

        for (decltype(I(NR)) r = 0; r < NR; ++r) {
            const auto cat = current_max[r];
            const auto num = current_num[r];
            IndexIn_ counter;

            switch(current_max[r]) {
                case Category::U8:
                    store8[chunk].ptr.push_back(sanisizer::sum<std::size_t>(store8[chunk].ptr.back(), num));
                    counter = counter8++;
                    identities8[chunk].push_back(r);
                    break;

                case Category::U16:
                    store16[chunk].ptr.push_back(sanisizer::sum<std::size_t>(store16[chunk].ptr.back(), num));
                    counter = counter16++;
                    identities16[chunk].push_back(r);
                    break;

                case Category::U32:
                    store32[chunk].ptr.push_back(sanisizer::sum<std::size_t>(store32[chunk].ptr.back(), num));
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
std::size_t get_sparse_ptr(
    const std::vector<Holder< std::uint8_t, IndexIn_, ColIndex_> >& store8,
    const std::vector<Holder<std::uint16_t, IndexIn_, ColIndex_> >& store16,
    const std::vector<Holder<std::uint32_t, IndexIn_, ColIndex_> >& store32,
    const std::vector<std::vector<Category> >& assigned_category,
    const std::vector<std::vector<IndexIn_> >& assigned_position,
    const IndexIn_ chunk,
    const IndexIn_ row)
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
    std::vector<Holder< std::uint8_t, IndexIn_, ColIndex_> >& store8,
    std::vector<Holder<std::uint16_t, IndexIn_, ColIndex_> >& store16,
    std::vector<Holder<std::uint32_t, IndexIn_, ColIndex_> >& store32,
    const Category cat,
    const IndexIn_ chunk, 
    const IndexIn_ col, 
    const ValueIn_ val,
    const std::size_t output_position) 
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
    std::vector<Holder< std::uint8_t, IndexIn_, ColIndex_> > store8,
    std::vector<Holder<std::uint16_t, IndexIn_, ColIndex_> > store16,
    std::vector<Holder<std::uint32_t, IndexIn_, ColIndex_> > store32,
    const IndexIn_ NR,
    const IndexIn_ chunk_size,
    const IndexIn_ leftovers)
{
    std::vector<std::shared_ptr<const tatami::Matrix<ValueOut_, IndexOut_> > > col_combined;
    IndexIn_ num_chunks = identities8.size();
    col_combined.reserve(num_chunks);

    for (decltype(I(num_chunks)) c = 0; c < num_chunks; ++c) {
        std::vector<std::shared_ptr<const tatami::Matrix<ValueOut_, IndexOut_> > > row_combined;
        IndexOut_ current_size = (c + 1 == num_chunks ? leftovers : chunk_size);
        row_combined.reserve(3);

        std::vector<IndexIn_> reordered(NR);
        IndexIn_ counter = 0;

        if (!(identities8[c].empty())) {
            row_combined.emplace_back(new tatami::CompressedSparseRowMatrix<ValueOut_, IndexOut_, std::vector<std::uint8_t>, std::vector<ColIndex_>, std::vector<std::size_t> >(
                identities8[c].size(), current_size, std::move(store8[c].value), std::move(store8[c].index), std::move(store8[c].ptr), false
            ));
            for (auto& i : identities8[c]) {
                reordered[i] = counter++;
            }
        }

        if (!(identities16[c].empty())) {
            row_combined.emplace_back(new tatami::CompressedSparseRowMatrix<ValueOut_, IndexOut_, std::vector<std::uint16_t>, std::vector<ColIndex_>, std::vector<std::size_t> >(
                identities16[c].size(), current_size, std::move(store16[c].value), std::move(store16[c].index), std::move(store16[c].ptr), false
            ));
            for (auto& i : identities16[c]) {
                reordered[i] = counter++;
            }
        }

        if (!(identities32[c].empty()) || row_combined.empty()) { // make sure that at least one matrix is fed into the DelayedBind.
            row_combined.emplace_back(new tatami::CompressedSparseRowMatrix<ValueOut_, IndexOut_, std::vector<uint32_t>, std::vector<ColIndex_>, std::vector<std::size_t> >(
                identities32[c].size(), current_size, std::move(store32[c].value), std::move(store32[c].index), std::move(store32[c].ptr), false
            ));
            for (auto& i : identities32[c]) {
                reordered[i] = counter++;
            }
        }

        if (row_combined.size() == 1) {
            col_combined.push_back(std::move(row_combined.front()));
        } else {
            col_combined.push_back(
                tatami::make_DelayedSubset<ValueOut_, IndexOut_>(
                    std::make_shared<const tatami::DelayedBind<ValueOut_, IndexOut_> >(std::move(row_combined), true),
                    std::move(reordered),
                    true
                )
            );
        }
    }

    return std::make_shared<tatami::DelayedBind<ValueOut_, IndexOut_> >(std::move(col_combined), false);
}

template<typename Output_, typename ColumnIndex_, typename Input_>
Output_ check_chunk_size(const Input_ chunk_size) {
    if (chunk_size <= 0) {
        throw std::runtime_error("chunk size should be positive");
    }
    constexpr auto max_index = std::numeric_limits<ColumnIndex_>::max();
    if (sanisizer::is_greater_than(chunk_size - 1, max_index)) { // same as chunk_size > max_index + 1 but this way avoids any risk of integer overflow.
        return sanisizer::sum<Output_>(max_index, 1);
    } else {
        return sanisizer::cast<Output_>(chunk_size);
    }
}

}

#endif
