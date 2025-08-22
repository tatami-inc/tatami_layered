#ifndef TATAMI_CONVERT_TO_LAYERED_SPARSE_HPP
#define TATAMI_CONVERT_TO_LAYERED_SPARSE_HPP

#include <cstdint>
#include <vector>
#include <memory>
#include <limits>
#include <algorithm>

#include "tatami/tatami.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils.hpp"

/**
 * @file convert_to_layered_sparse.hpp
 * @brief Create a layered sparse matrix for non-negative integers.
 */

namespace tatami_layered {

/**
 * @cond
 */
template<typename ColIndex_, typename ValueOut_ = double, typename IndexOut_ = int, typename ValueIn_, typename IndexIn_>
std::shared_ptr<tatami::Matrix<ValueOut_, IndexOut_> > convert_by_row(const tatami::Matrix<ValueIn_, IndexIn_>& mat, IndexIn_ chunk_size, int nthreads) {
    auto NR = mat.nrow(), NC = mat.ncol();
    IndexIn_ leftovers = NC % chunk_size;
    IndexIn_ nchunks = sanisizer::max(1, NC / chunk_size + (leftovers != 0));

    auto store8  = tatami::create_container_of_Index_size<std::vector<Holder< std::uint8_t, IndexOut_, ColIndex_> > >(nchunks);
    auto store16 = tatami::create_container_of_Index_size<std::vector<Holder<std::uint16_t, IndexOut_, ColIndex_> > >(nchunks);
    auto store32 = tatami::create_container_of_Index_size<std::vector<Holder<std::uint32_t, IndexOut_, ColIndex_> > >(nchunks);

    auto identities8  = tatami::create_container_of_Index_size<std::vector<std::vector<IndexOut_> > >(nchunks);
    auto identities16 = tatami::create_container_of_Index_size<std::vector<std::vector<IndexOut_> > >(nchunks);
    auto identities32 = tatami::create_container_of_Index_size<std::vector<std::vector<IndexOut_> > >(nchunks);

    auto assigned_position = tatami::create_container_of_Index_size<std::vector<std::vector<IndexOut_> > >(nchunks);
    auto assigned_category = tatami::create_container_of_Index_size<std::vector<std::vector<Category> > >(nchunks);

    // First pass to define the allocations.
    {
        auto max_per_chunk = tatami::create_container_of_Index_size<std::vector<std::vector<Category> > >(nchunks);
        for (auto& x : max_per_chunk) {
            tatami::resize_container_to_Index_size(x, NR);
        }

        auto num_per_chunk = tatami::create_container_of_Index_size<std::vector<std::vector<IndexIn_> > >(nchunks);
        for (auto& x : num_per_chunk) {
            tatami::resize_container_to_Index_size(x, NR);
        }

        if (mat.sparse()) {
            tatami::parallelize([&](int, IndexIn_ start, IndexIn_ length) -> void {
                auto ext = tatami::consecutive_extractor<true>(mat, true, start, length, [&]{
                    tatami::Options opt;
                    opt.sparse_ordered_index = false;
                    return opt;
                }());
                auto dbuffer = tatami::create_container_of_Index_size<std::vector<ValueIn_> >(NC);
                auto ibuffer = tatami::create_container_of_Index_size<std::vector<IndexIn_> >(NC);

                for (IndexIn_ r = start, end = start + length; r < end; ++r) {
                    auto range = ext->fetch(r, dbuffer.data(), ibuffer.data());
                    for (IndexIn_ i = 0; i < range.number; ++i) {
                        if (range.value[i]) {
                            auto chunk = range.index[i] / chunk_size;
                            auto cat = categorize(range.value[i]);
                            max_per_chunk[chunk][r] = std::max(max_per_chunk[chunk][r], cat);
                            ++num_per_chunk[chunk][r];
                        }
                    }
                }
            }, NR, nthreads);

        } else {
            tatami::parallelize([&](int, IndexIn_ start, IndexIn_ length) -> void {
                auto ext = tatami::consecutive_extractor<false>(mat, true, start, length);
                auto dbuffer = tatami::create_container_of_Index_size<std::vector<ValueIn_> >(NC);

                for (IndexIn_ r = start, end = start + length; r < end; ++r) {
                    auto ptr = ext->fetch(r, dbuffer.data());
                    for (IndexIn_ c = 0; c < NC; ++c) {
                        if (ptr[c]) {
                            auto chunk = c / chunk_size;
                            auto cat = categorize(ptr[c]);
                            max_per_chunk[chunk][r] = std::max(max_per_chunk[chunk][r], cat);
                            ++num_per_chunk[chunk][r];
                        }
                    }
                }
            }, NR, nthreads);
        }

        allocate_rows(
            max_per_chunk, 
            num_per_chunk, 
            identities8, 
            identities16, 
            identities32, 
            store8, 
            store16, 
            store32, 
            assigned_category, 
            assigned_position
        );
    }

    // Second pass to actually fill the vectors.
    {
        tatami::parallelize([&](int, IndexIn_ start, IndexIn_ length) -> void {
            auto output_positions = tatami::create_container_of_Index_size<std::vector<std::size_t> >(nchunks);
            auto dbuffer = tatami::create_container_of_Index_size<std::vector<ValueIn_> >(NC);

            if (mat.sparse()) {
                auto ibuffer = tatami::create_container_of_Index_size<std::vector<IndexIn_> >(NC);
                auto ext = tatami::consecutive_extractor<true>(mat, true, start, length);

                for (IndexIn_ r = start, end = start + length; r < end; ++r) {
                    for (decltype(nchunks) chunk = 0; chunk < nchunks; ++chunk) {
                        output_positions[chunk] = get_sparse_ptr(store8, store16, store32, assigned_category, assigned_position, chunk, r);
                    }

                    auto range = ext->fetch(r, dbuffer.data(), ibuffer.data());
                    for (IndexIn_ i = 0; i < range.number; ++i) {
                        if (range.value[i]) {
                            IndexIn_ chunk = range.index[i] / chunk_size;
                            IndexIn_ col = range.index[i] % chunk_size;
                            fill_sparse_value(store8, store16, store32, assigned_category[chunk][r], chunk, col, range.value[i], output_positions[chunk]++);
                        }
                    }
                }

            } else {
                auto ext = tatami::consecutive_extractor<false>(mat, true, start, length);

                for (IndexIn_ r = start, end = start + length; r < end; ++r) {
                    for (decltype(nchunks) chunk = 0; chunk < nchunks; ++chunk) {
                        output_positions[chunk] = get_sparse_ptr(store8, store16, store32, assigned_category, assigned_position, chunk, r);
                    }

                    auto ptr = ext->fetch(r, dbuffer.data());
                    for (IndexIn_ c = 0; c < NC; ++c) {
                        if (ptr[c]) {
                            IndexIn_ chunk = c / chunk_size;
                            IndexIn_ col = c % chunk_size;
                            fill_sparse_value(store8, store16, store32, assigned_category[chunk][r], chunk, col, ptr[c], output_positions[chunk]++);
                        }
                    }
                }
            }

        }, NR, nthreads);
    }

    return consolidate_matrices<ValueOut_, IndexOut_>(
        identities8, 
        identities16, 
        identities32, 
        std::move(store8), 
        std::move(store16), 
        std::move(store32),
        NR,
        chunk_size,
        leftovers
    );
}

template<typename ColIndex_, typename ValueOut_ = double, typename IndexOut_ = int, typename ValueIn_, typename IndexIn_>
std::shared_ptr<tatami::Matrix<ValueOut_, IndexOut_> > convert_by_column(const tatami::Matrix<ValueIn_, IndexIn_>& mat, IndexIn_ chunk_size, int nthreads) {
    auto NR = mat.nrow(), NC = mat.ncol();
    IndexIn_ leftovers = NC % chunk_size;
    IndexIn_ nchunks = sanisizer::max(1, NC / chunk_size + (leftovers != 0));

    auto store8  = tatami::create_container_of_Index_size<std::vector<Holder< std::uint8_t, IndexOut_, ColIndex_> > >(nchunks);
    auto store16 = tatami::create_container_of_Index_size<std::vector<Holder<std::uint16_t, IndexOut_, ColIndex_> > >(nchunks);
    auto store32 = tatami::create_container_of_Index_size<std::vector<Holder<std::uint32_t, IndexOut_, ColIndex_> > >(nchunks);

    auto identities8  = tatami::create_container_of_Index_size<std::vector<std::vector<IndexOut_> > >(nchunks);
    auto identities16 = tatami::create_container_of_Index_size<std::vector<std::vector<IndexOut_> > >(nchunks);
    auto identities32 = tatami::create_container_of_Index_size<std::vector<std::vector<IndexOut_> > >(nchunks);

    auto assigned_position = tatami::create_container_of_Index_size<std::vector<std::vector<IndexOut_> > >(nchunks);
    auto assigned_category = tatami::create_container_of_Index_size<std::vector<std::vector<Category> > >(nchunks);

    // First pass to define the allocations.
    {
        auto max_per_chunk_threaded = sanisizer::create<std::vector<std::vector<std::vector<Category> > > >(nthreads);
        for (auto& max_per_chunk : max_per_chunk_threaded) { 
            tatami::resize_container_to_Index_size<std::vector<std::vector<Category> > >(max_per_chunk, nchunks);
            for (auto& x : max_per_chunk) {
                tatami::resize_container_to_Index_size(x, NR);
            }
        }

        auto num_per_chunk_threaded = sanisizer::create<std::vector<std::vector<std::vector<IndexIn_> > > >(nthreads);
        for (auto& num_per_chunk : num_per_chunk_threaded) { 
            tatami::resize_container_to_Index_size<std::vector<std::vector<IndexIn_> > >(num_per_chunk, nchunks);
            for (auto& x : num_per_chunk) {
                tatami::resize_container_to_Index_size(x, NR);
            }
        }

        if (mat.sparse()) {
            tatami::parallelize([&](int t, IndexIn_ start, IndexIn_ length) -> void {
                auto ext = tatami::consecutive_extractor<true>(&mat, false, start, length, [&]{
                    tatami::Options opt;
                    opt.sparse_ordered_index = false;
                    return opt;
                }());
                auto dbuffer = tatami::create_container_of_Index_size<std::vector<ValueIn_> >(NR);
                auto ibuffer = tatami::create_container_of_Index_size<std::vector<IndexIn_> >(NR);

                auto& max_per_chunk = max_per_chunk_threaded[t];
                auto& num_per_chunk = num_per_chunk_threaded[t];

                for (IndexIn_ c = start, end = start + length; c < end; ++c) {
                    auto range = ext->fetch(c, dbuffer.data(), ibuffer.data());
                    auto chunk = c / chunk_size;
                    auto& max_vec = max_per_chunk[chunk];
                    auto& num_vec = num_per_chunk[chunk];

                    for (IndexIn_ i = 0; i < range.number; ++i) {
                        if (range.value[i]) {
                            auto cat = categorize(range.value[i]);
                            auto r = range.index[i];
                            max_vec[r] = std::max(max_vec[r], cat);
                            ++num_vec[r];
                        }
                    }
                }
            }, NC, nthreads);

        } else {
            tatami::parallelize([&](int t, IndexIn_ start, IndexIn_ length) -> void {
                auto ext = tatami::consecutive_extractor<false>(&mat, false, start, length);
                auto dbuffer = tatami::create_container_of_Index_size<std::vector<ValueIn_> >(NR);

                auto& max_per_chunk = max_per_chunk_threaded[t];
                auto& num_per_chunk = num_per_chunk_threaded[t];

                for (IndexIn_ c = start, end = start + length; c < end; ++c) {
                    auto ptr = ext->fetch(c, dbuffer.data());
                    auto chunk = c / chunk_size;
                    auto& max_vec = max_per_chunk[chunk];
                    auto& num_vec = num_per_chunk[chunk];

                    for (IndexIn_ r = 0; r < NR; ++r) {
                        if (ptr[r]) {
                            auto cat = categorize(ptr[r]);
                            max_vec[r] = std::max(max_vec[r], cat);
                            ++num_vec[r];
                        }
                    }
                }
            }, NC, nthreads);
        }

        auto max_per_chunk = tatami::create_container_of_Index_size<std::vector<std::vector<Category> > >(nchunks);
        auto num_per_chunk = tatami::create_container_of_Index_size<std::vector<std::vector<IndexIn_> > >(nchunks);

        for (decltype(nchunks) chunk = 0; chunk < nchunks; ++chunk) {
            // Assume we have at least one thread!
            max_per_chunk[chunk].swap(max_per_chunk_threaded[0][chunk]);
            num_per_chunk[chunk].swap(num_per_chunk_threaded[0][chunk]);

            for (int t = 1; t < nthreads; ++t) {
                for (IndexIn_ r = 0; r < NR; ++r) {
                    max_per_chunk[chunk][r] = std::max(max_per_chunk[chunk][r], max_per_chunk_threaded[t][chunk][r]);
                    num_per_chunk[chunk][r] += num_per_chunk_threaded[t][chunk][r];
                }
            }
        }

        allocate_rows(
            max_per_chunk, 
            num_per_chunk, 
            identities8, 
            identities16, 
            identities32, 
            store8, 
            store16, 
            store32, 
            assigned_category, 
            assigned_position
        );
    }

    // Second pass to actually fill the vectors.
    {
        tatami::parallelize([&](int, IndexIn_ start, IndexIn_ length) -> void {
            auto output_positions = tatami::create_container_of_Index_size<std::vector<std::vector<std::size_t> > >(nchunks);
            for (decltype(nchunks) chunk = 0; chunk < nchunks; ++chunk) {
                output_positions[chunk].resize(length);
                for (IndexIn_ r = 0; r < length; ++r) {
                    output_positions[chunk][r] = get_sparse_ptr(store8, store16, store32, assigned_category, assigned_position, chunk, r + start);
                }
            }

            auto dbuffer = tatami::create_container_of_Index_size<std::vector<ValueIn_> >(length);

            if (mat.sparse()) {
                auto ibuffer = tatami::create_container_of_Index_size<std::vector<IndexIn_> >(length);
                auto ext = tatami::consecutive_extractor<true>(mat, false, static_cast<IndexIn_>(0), NC, start, length);

                for (IndexIn_ c = 0; c < NC; ++c) {
                    auto range = ext->fetch(c, dbuffer.data(), ibuffer.data());
                    auto chunk = c / chunk_size;
                    IndexIn_ col = c % chunk_size;
                    auto& outpos = output_positions[chunk];

                    for (IndexIn_ i = 0; i < range.number; ++i) {
                        if (range.value[i]) {
                            auto r = range.index[i];
                            fill_sparse_value(store8, store16, store32, assigned_category[chunk][r], chunk, col, range.value[i], outpos[r - start]++);
                        }
                    }
                }

            } else {
                auto ext = tatami::consecutive_extractor<false>(mat, false, static_cast<IndexIn_>(0), NC, start, length);

                for (IndexIn_ c = 0; c < NC; ++c) {
                    auto ptr = ext->fetch(c, dbuffer.data());
                    auto chunk = c / chunk_size;
                    IndexIn_ col = c % chunk_size;
                    auto& outpos = output_positions[chunk];

                    for (IndexIn_ r = 0; r < NR; ++r) {
                        if (ptr[r]) {
                            fill_sparse_value(store8, store16, store32, assigned_category[chunk][r], chunk, col, ptr[r], outpos[r - start]++);
                        }
                    }
                }
            }

        }, NR, nthreads);
    }

    return consolidate_matrices<ValueOut_, IndexOut_>(
        identities8, 
        identities16, 
        identities32, 
        std::move(store8), 
        std::move(store16), 
        std::move(store32),
        NR,
        chunk_size,
        leftovers
    );
}
/**
 * @endcond
 */

/**
 * @brief Options for `convert_to_layered_sparse()`.
 */
struct ConvertToLayeredSparseOptions {
    /**
     * Chunk size to use for partitioning columns.
     * This should be a positive integer.
     */
    std::size_t chunk_size = sanisizer::cap<std::size_t>(65536);

    /**
     * Number of threads to use.
     * This should be a positive integer.
     */
    int num_threads = 1;
};

/**
 * @param mat A `tatami::Matrix` object containing non-negative integers.
 * @param options Further options.
 *
 * @return A `tatami::Matrix` object containing a layered sparse matrix.
 *
 * @tparam ValueOut_ Type of data value for the output `tatami::Matrix` interface.
 * @tparam IndexOut_ Integer type for the row/column indices of the output.
 * @tparam ColumnIndex_ Integer type for the stored column indices.
 * @tparam ValueIn_ Type of data value for the input.
 * @tparam IndexIn_ Integer type for the row/column indices of the input.
 *
 * This function converts an existing sparse matrix of non-negative integers into a layered sparse matrix.
 * The aim is to reduce memory usage by storing each row's data in the smallest unsigned integer type that can hold them.
 * To create a layered sparse matrix:
 *
 * 1. We split the input matrix into chunks of `options.chunk_size` contiguous columns.
 * 2. Within each chunk, we identify the maximum integer for each row.
 * 3. Data for each row are stored in one of three sparse matrices using 8, 16, or 32-bit unsigned integers as the data type, depending on the row's maximum value.
 * 4. The three sparse matrices are combined together with `tatami::DelayedBind`, with some use of `tatami::DelayedSubset` to restore the original order of rows.
 * 5. Each chunk is then combined together with `tatami::DelayedBind`.
 *
 * We improve the chances of being able to use small types by splitting the matrix columns into chunks.
 * This ensures that a few large values in a particular row only cause promotion to a larger integer type for the chunks in which they occur.
 *
 * Setting `ColumnIndex_` to the smallest type that can hold `options.chunk_size - 1` can be used to further reduce memory usage.
 * If `ColumnIndex_` is not able to hold `options.chunk_size - 1`, the chunk size is automatically set to the largest value that can be represented by `ColumnIndex_` plus 1.
 * For example, if `ColumnIndex_` was set to an unsigned 8-bit integer, `chunk_size` would be automatically reduced to 256.
 */
template<typename ValueOut_ = double, typename IndexOut_ = int, typename ColumnIndex_ = std::uint16_t, typename ValueIn_, typename IndexIn_>
std::shared_ptr<tatami::Matrix<ValueOut_, IndexOut_> > convert_to_layered_sparse(const tatami::Matrix<ValueIn_, IndexIn_>& mat, const ConvertToLayeredSparseOptions& options) {
    IndexIn_ chunk_size = check_chunk_size<IndexIn_, ColumnIndex_>(options.chunk_size);
    if (mat.prefer_rows()) {
        return convert_by_row<ColumnIndex_, ValueOut_, IndexOut_>(mat, chunk_size, options.num_threads);
    } else {
        return convert_by_column<ColumnIndex_, ValueOut_, IndexOut_>(mat, chunk_size, options.num_threads);
    }
}

/**
 * @cond
 */
// Provided for back-compatibility.
template<typename ValueOut_ = double, typename IndexOut_ = int, typename ColumnIndex_ = std::uint16_t, typename ValueIn_, typename IndexIn_>
std::shared_ptr<tatami::Matrix<ValueOut_, IndexOut_> > convert_to_layered_sparse(const tatami::Matrix<ValueIn_, IndexIn_>& mat, IndexIn_ chunk_size = 65536, int num_threads = 1) {
    return convert_to_layered_sparse(mat, [&]{
        ConvertToLayeredSparseOptions opt;
        opt.chunk_size = chunk_size;
        opt.num_threads = num_threads;
        return opt;
    }());
}

template<typename ValueOut_ = double, typename IndexOut_ = int, typename ColumnIndex_ = std::uint16_t, typename ValueIn_, typename IndexIn_>
std::shared_ptr<tatami::Matrix<ValueOut_, IndexOut_> > convert_to_layered_sparse(const tatami::Matrix<ValueIn_, IndexIn_>* mat, IndexIn_ chunk_size = 65536, int num_threads = 1) {
    return convert_to_layered_sparse<ValueOut_, IndexOut_, ColumnIndex_, ValueIn_, IndexIn_>(*mat, chunk_size, num_threads);
}
/**
 * @endcond
 */

}

#endif
