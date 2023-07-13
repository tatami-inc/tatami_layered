#ifndef TATAMI_CONVERT_TO_LAYERED_SPARSE_HPP
#define TATAMI_CONVERT_TO_LAYERED_SPARSE_HPP

#include "tatami/tatami.hpp"
#include "utils.hpp"

#include <cstdint>
#include <vector>
#include <memory>
#include <limits>
#include <algorithm>

/**
 * @file convert_to_layered_sparse.hpp
 * @brief Create a layered sparse matrix for non-negative integers.
 */

namespace tatami_layered {

/**
 * @cond
 */
template<typename ColIndex_, typename ValueOut_ = double, typename IndexOut_ = int, typename ValueIn_, typename IndexIn_>
std::shared_ptr<tatami::Matrix<ValueOut_, IndexOut_> > convert_by_row(const tatami::Matrix<ValueIn_, IndexIn_>* mat, IndexIn_ chunk_size, int nthreads) {
    auto NR = mat->nrow(), NC = mat->ncol();
    IndexIn_ leftovers = NC % chunk_size;
    size_t nchunks = std::max(static_cast<size_t>(1), static_cast<size_t>(NC) / chunk_size + (leftovers != 0));

    std::vector<Holder<uint8_t, IndexOut_, ColIndex_> > store8(nchunks);
    std::vector<Holder<uint16_t, IndexOut_, ColIndex_> > store16(nchunks);
    std::vector<Holder<uint32_t, IndexOut_, ColIndex_> > store32(nchunks);

    std::vector<std::vector<IndexOut_> > identities8(nchunks), identities16(nchunks), identities32(nchunks);
    std::vector<std::vector<IndexOut_> > assigned_position(nchunks);
    std::vector<std::vector<Category> > assigned_category(nchunks);

    // First pass to define the allocations.
    {
        std::vector<std::vector<Category> > max_per_chunk(nchunks);
        std::vector<std::vector<IndexIn_> > num_per_chunk(nchunks);
        for (auto& x : max_per_chunk) { x.resize(NR); }
        for (auto& x : num_per_chunk) { x.resize(NR); }

        if (mat->sparse()) {
            tatami::parallelize([&](size_t, IndexIn_ start, IndexIn_ length) -> void {
                tatami::Options opt;
                opt.sparse_ordered_index = false;
                auto ext = tatami::consecutive_extractor<true, true>(mat, start, length, opt);
                std::vector<ValueIn_> dbuffer(NC);
                std::vector<IndexIn_> ibuffer(NC);

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
            tatami::parallelize([&](size_t, IndexIn_ start, IndexIn_ length) -> void {
                auto ext = tatami::consecutive_extractor<true, false>(mat, start, length);
                std::vector<ValueIn_> dbuffer(NC);

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
        tatami::parallelize([&](size_t, IndexIn_ start, IndexIn_ length) -> void {
            std::vector<size_t> output_positions(nchunks);
            std::vector<ValueIn_> dbuffer(NC);

            if (mat->sparse()) {
                std::vector<IndexIn_> ibuffer(NC);
                auto ext = tatami::consecutive_extractor<true, true>(mat, start, length);
                for (IndexIn_ r = start, end = start + length; r < end; ++r) {
                    for (size_t chunk = 0; chunk < nchunks; ++chunk) {
                        output_positions[chunk] = get_sparse_ptr(store8, store16, store32, assigned_category, assigned_position, chunk, r);
                    }

                    auto range = ext->fetch(r, dbuffer.data(), ibuffer.data());
                    for (IndexIn_ i = 0; i < range.number; ++i) {
                        if (range.value[i]) {
                            size_t chunk = range.index[i] / chunk_size;
                            IndexIn_ col = range.index[i] % chunk_size;
                            fill_sparse_value(store8, store16, store32, assigned_category[chunk][r], chunk, col, range.value[i], output_positions[chunk]++);
                        }
                    }
                }

            } else {
                auto ext = tatami::consecutive_extractor<true, false>(mat, start, length);
                for (IndexIn_ r = start, end = start + length; r < end; ++r) {
                    for (size_t chunk = 0; chunk < nchunks; ++chunk) {
                        output_positions[chunk] = get_sparse_ptr(store8, store16, store32, assigned_category, assigned_position, chunk, r);
                    }

                    auto ptr = ext->fetch(r, dbuffer.data());
                    for (IndexIn_ c = 0; c < NC; ++c) {
                        if (ptr[c]) {
                            size_t chunk = c / chunk_size;
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
std::shared_ptr<tatami::Matrix<ValueOut_, IndexOut_> > convert_by_column(const tatami::Matrix<ValueIn_, IndexIn_>* mat, IndexIn_ chunk_size, int nthreads) {
    auto NR = mat->nrow(), NC = mat->ncol();
    IndexIn_ leftovers = NC % chunk_size;
    size_t nchunks = std::max(static_cast<size_t>(1), static_cast<size_t>(NC) / chunk_size + (leftovers != 0));
    std::cout << nchunks << std::endl;

    std::vector<Holder<uint8_t, IndexOut_, ColIndex_> > store8(nchunks);
    std::vector<Holder<uint16_t, IndexOut_, ColIndex_> > store16(nchunks);
    std::vector<Holder<uint32_t, IndexOut_, ColIndex_> > store32(nchunks);

    std::vector<std::vector<IndexOut_> > identities8(nchunks), identities16(nchunks), identities32(nchunks);
    std::vector<std::vector<IndexOut_> > assigned_position(nchunks);
    std::vector<std::vector<Category> > assigned_category(nchunks);

    // First pass to define the allocations.
    {
        std::vector<std::vector<std::vector<Category> > > max_per_chunk_threaded(nthreads);
        for (auto& max_per_chunk : max_per_chunk_threaded) { 
            max_per_chunk.resize(nchunks);
            for (auto& x : max_per_chunk) { x.resize(NR); }
        }

        std::vector<std::vector<std::vector<IndexIn_> > > num_per_chunk_threaded(nthreads);
        for (auto& num_per_chunk : num_per_chunk_threaded) { 
            num_per_chunk.resize(nchunks);
            for (auto& x : num_per_chunk) { x.resize(NR); }
        }

        if (mat->sparse()) {
            tatami::parallelize([&](size_t t, IndexIn_ start, IndexIn_ length) -> void {
                tatami::Options opt;
                opt.sparse_ordered_index = false;
                auto ext = tatami::consecutive_extractor<false, true>(mat, start, length, opt);
                std::vector<ValueIn_> dbuffer(NR);
                std::vector<IndexIn_> ibuffer(NR);
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
            tatami::parallelize([&](size_t t, IndexIn_ start, IndexIn_ length) -> void {
                auto ext = tatami::consecutive_extractor<false, false>(mat, start, length);
                std::vector<ValueIn_> dbuffer(NR);
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

        std::vector<std::vector<Category> > max_per_chunk(nchunks);
        std::vector<std::vector<IndexIn_> > num_per_chunk(nchunks);

        for (size_t chunk = 0; chunk < nchunks; ++chunk) {
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
        tatami::parallelize([&](size_t, IndexIn_ start, IndexIn_ length) -> void {
            std::vector<ValueIn_> dbuffer(length);

            std::vector<std::vector<size_t> > output_positions(nchunks);
            for (size_t chunk = 0; chunk < nchunks; ++chunk) {
                output_positions[chunk].resize(length);
                for (IndexIn_ r = 0; r < length; ++r) {
                    output_positions[chunk][r] = get_sparse_ptr(store8, store16, store32, assigned_category, assigned_position, chunk, r + start);
                }
            }

            if (mat->sparse()) {
                std::vector<IndexIn_> ibuffer(length);
                auto ext = tatami::consecutive_extractor<false, true>(mat, 0, NC, start, length);

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
                auto ext = tatami::consecutive_extractor<false, false>(mat, 0, NC, start, length);

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
 * @param mat A `tatami::Matrix` object containing non-negative integers.
 * @param chunk_size Chunk size to use for partitioning columns.
 * @param num_threads Number of threads to use.
 *
 * @return A `tatami::Matrix` object containing a layered sparse matrix.
 *
 * @tparam ValueOut_ Type of data value for the output `tatami::Matrix` interface.
 * @tparam IndexOut_ Integer type for the row/column indices of the output.
 * @tparam ColumnIndex_ Integer type for the stored column indices.
 * @tparam ValueIn_ Type of data value for the input.
 * @tparam IndexIn_ Integer type for the row/column indices of the input.
 *
 * This function converts an existing sparse integer matrix into a layered sparse matrix.
 * The aim is to reduce memory usage by storing each gene's counts in the smallest unsigned integer type that can hold them.
 * We improve the chances of being able to use small types by splitting the matrix columns into contiguous chunks according to `chunk_size`.
 * This ensures that a few large values in a particular row only cause promotion to a larger type for the chunks in which they occur.
 * Setting `ColumnIndex_` can be used to further reduce memory usage, but this must be chosen such that its maximum value is not less than the chunk size (minus 1).
 */
template<typename ValueOut_ = double, typename IndexOut_ = int, typename ColumnIndex_ = uint16_t, typename ValueIn_, typename IndexIn_>
std::shared_ptr<tatami::Matrix<ValueOut_, IndexOut_> > convert_to_layered_sparse(const tatami::Matrix<ValueIn_, IndexIn_>* mat, IndexIn_ chunk_size = 65536, int num_threads = 1) {
    if (static_cast<size_t>(chunk_size) > static_cast<size_t>(std::numeric_limits<ColumnIndex_>::max()) + 1) {
        throw std::runtime_error("'chunk_size' is too large for the specified column index type");
    }

    if (mat->prefer_rows()) {
        return convert_by_row<ColumnIndex_, ValueOut_, IndexOut_>(mat, chunk_size, num_threads);
    } else {
        return convert_by_column<ColumnIndex_, ValueOut_, IndexOut_>(mat, chunk_size, num_threads);
    }
}

}

#endif
