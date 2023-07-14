#ifndef TATAMI_LAYERED_READ_LAYERED_SPARSE_FROM_MATRIX_MARKET_HPP
#define TATAMI_LAYERED_READ_LAYERED_SPARSE_FROM_MATRIX_MARKET_HPP

#include "byteme/byteme.hpp"
#include "eminem/eminem.hpp"

#include <vector>
#include <algorithm>

#include "utils.hpp"

namespace tatami_layered {

/**
 * @cond
 */
template<typename Value_, typename Index_, typename ColumnIndex_, class Creator_>
std::shared_ptr<tatami::Matrix<Value_, Index_> > read_layered_sparse_from_matrix_market(Creator_ create, Index_ chunk_size) {
    size_t NR, NC, nchunks, leftovers;

    std::vector<Holder<uint8_t, Index_, ColumnIndex_> > store8;
    std::vector<Holder<uint16_t, Index_, ColumnIndex_> > store16;
    std::vector<Holder<uint32_t, Index_, ColumnIndex_> > store32;

    std::vector<std::vector<Index_> > identities8, identities16, identities32;
    std::vector<std::vector<Index_> > assigned_position;
    std::vector<std::vector<Category> > assigned_category;

    // First pass, scanning for the max and number.
    {
        auto reader = create();
        eminem::Parser parser(&reader);

        parser.scan_preamble();
        NR = parser.get_nrows();
        NC = parser.get_ncols();
        nchunks = std::max(static_cast<size_t>(1), NC / chunk_size + (leftovers != 0));
        leftovers = NC % chunk_size;

        store8.resize(nchunks);
        store16.resize(nchunks);
        store32.resize(nchunks);
        identities8.resize(nchunks);
        identities16.resize(nchunks);
        identities32.resize(nchunks);
        assigned_position.resize(nchunks);
        assigned_category.resize(nchunks);

        std::vector<std::vector<Category> > max_per_chunk(nchunks);
        std::vector<std::vector<Index_> > num_per_chunk(nchunks);
        for (auto& x : max_per_chunk) { x.resize(NR); }
        for (auto& x : num_per_chunk) { x.resize(NR); }

        auto handler = [&](size_t r, size_t c, int val) -> void {
            auto chunk = (c - 1) / chunk_size;
            --r;
            max_per_chunk[chunk][r] = std::max(max_per_chunk[chunk][r], categorize(val));
            ++num_per_chunk[chunk][r];
        };

        const auto& banner = parser.get_banner();
        if (banner.field == eminem::Field::INTEGER) {
            parser.scan_integer(handler);
        } else if (banner.field == eminem::Field::DOUBLE || banner.field == eminem::Field::REAL) {
            parser.scan_real(handler);
        } else {
            throw std::runtime_error("expected a numeric field in the Matrix Market file");
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
   
    // Now allocating.
    {
        std::vector<std::vector<size_t> > output_positions(nchunks);
        for (size_t chunk = 0; chunk < nchunks; ++chunk) {
            output_positions[chunk].resize(NR);
            for (size_t r = 0; r < NR; ++r) {
                output_positions[chunk][r] = get_sparse_ptr(store8, store16, store32, assigned_category, assigned_position, chunk, static_cast<Index_>(r));
            }
        }

        auto handler = [&](size_t r, size_t c, int val) -> void {
            --c;
            auto chunk = c / chunk_size;
            Index_ offset = c % chunk_size;
            --r;
            fill_sparse_value(store8, store16, store32, assigned_category[chunk][r], chunk, offset, val, output_positions[chunk][r]++);
        };

        auto reader = create();
        eminem::Parser parser(&reader);
        parser.scan_preamble();

        const auto& banner = parser.get_banner();
        if (banner.field == eminem::Field::INTEGER) {
            parser.scan_integer(handler);
        } else if (banner.field == eminem::Field::DOUBLE || banner.field == eminem::Field::REAL) {
            parser.scan_real(handler);
        }

        // Checking that the column indices are sorted properly.
        auto sorter = [&](auto& store) -> void {
            std::vector<std::pair<typename decltype(store[0].index)::value_type, typename decltype(store[0].value)::value_type> > buffer;
            buffer.reserve(chunk_size);

            for (auto& st : store) {
                size_t current_NR = st.ptr.size() - 1;
                for (size_t r = 0; r < current_NR; ++r) {
                    size_t start = st.ptr[r], end = st.ptr[r + 1];

                    if (!std::is_sorted(st.index.begin() + start, st.index.begin() + end)) {
                        buffer.clear();
                        for (size_t i = start; i < end; ++i) {
                            buffer.emplace_back(st.index[i], st.value[i]);
                        }

                        std::sort(buffer.begin(), buffer.end());
                        auto bIt = buffer.begin();
                        for (size_t i = start; i < end; ++i, ++bIt) {
                            st.index[i] = bIt->first;
                            st.value[i] = bIt->second;
                        }
                    }
                }
            }
        };

        sorter(store8);
        sorter(store16);
        sorter(store32);
    }

    return consolidate_matrices<Value_, Index_>(
        identities8, 
        identities16, 
        identities32, 
        std::move(store8), 
        std::move(store16), 
        std::move(store32),
        static_cast<Index_>(NR),
        chunk_size,
        static_cast<Index_>(leftovers)
    );
}
/**
 * @endcond
 */

/**
 * @param filepath Path to a (possibly Gzip-compressed) Matrix Market file.
 * The file should contain integer data in the coordinate format, stored in text without any compression.
 * @param chunk_size Chunk size to use for partitioning columns.
 * @param compression Compression method for the file - no compression (0) or Gzip compression (1).
 * If set to -1, the function will automatically guess the compression based on magic numbers.
 * @param buffer_size Size of the buffer (in bytes) to use when reading from file.
 * 
 * @return A `tatami::Matrix` object containing a layered sparse matrix.
 *
 * @tparam Value_ Type of data value for the output `tatami::Matrix` interface.
 * @tparam Index_ Integer type for the row/column indices of the output.
 * @tparam ColumnIndex_ Integer type for the stored column indices.
 *
 * This function loads a layered sparse integer matrix from a Matrix Market file.
 * The aim is to reduce memory usage by storing each gene's counts in the smallest unsigned integer type that can hold them.
 * See `convert_to_layered_sparse()` for more details.
 */
template<typename Value_ = double, typename Index_ = int, typename ColumnIndex_ = uint16_t>
std::shared_ptr<tatami::Matrix<Value_, Index_> > read_layered_sparse_from_matrix_market_file(
    const char * filepath, 
    Index_ chunk_size = 65536, 
    int compression = -1, 
    size_t buffer_size = 65536)
{
    if (compression != 0) {
#if __has_include("zlib.h")
        if (compression == -1) {
            return read_layered_sparse_from_matrix_market<Value_, Index_, ColumnIndex_>([&]() -> auto { return byteme::SomeFileReader(filepath, buffer_size); }, chunk_size);
        } else if (compression == 1) {
            return read_layered_sparse_from_matrix_market<Value_, Index_, ColumnIndex_>([&]() -> auto { return byteme::GzipFileReader(filepath, buffer_size); }, chunk_size);
        }
#else
        if (compression != -1) {
            throw std::runtime_error("layered sparse matrix reader not compiled with support for non-zero 'compression'");
        }
#endif
    }
    return read_layered_sparse_from_matrix_market<Value_, Index_, ColumnIndex_>([&]() -> auto { return byteme::RawFileReader(filepath, buffer_size); }, chunk_size);
}


/**
 * @param contents Array containing the contents of a (possibly Gzip-compressed) Matrix Market file.
 * The file should contain integer data in the coordinate format, stored in text without any compression.
 * @param length Length of the array.
 * @param chunk_size Chunk size to use for partitioning columns.
 * @param compression Compression method for the file contents - no compression (0) or Gzip/Zlib compression (1).
 * If set to -1, the function will automatically guess the compression based on magic numbers.
 * @param buffer_size Size of the buffer to use for decompression, in bytes.
 * 
 * @return A `tatami::Matrix` object containing a layered sparse matrix.
 *
 * @tparam Value_ Type of data value for the output `tatami::Matrix` interface.
 * @tparam Index_ Integer type for the row/column indices of the output.
 * @tparam ColumnIndex_ Integer type for the stored column indices.
 *
 * This function loads a layered sparse integer matrix from a buffer with the contents of a Matrix Market file.
 * The aim is to reduce memory usage by storing each gene's counts in the smallest unsigned integer type that can hold them.
 * See `convert_to_layered_sparse()` for more details.
 */
template<typename Value_ = double, typename Index_ = int, typename ColumnIndex_ = uint16_t>
std::shared_ptr<tatami::Matrix<Value_, Index_> > read_layered_sparse_from_matrix_market_buffer(
    const unsigned char* contents, 
    size_t length, 
    Index_ chunk_size = 65536, 
    int compression = -1, 
    size_t buffer_size = 65536) 
{
    if (compression != 0) {
#if __has_include("zlib.h")
        if (compression == -1) {
            return read_layered_sparse_from_matrix_market<Value_, Index_, ColumnIndex_>([&]() -> auto { return byteme::SomeBufferReader(contents, length, buffer_size); }, chunk_size);
        } else if (compression == 1) {
            return read_layered_sparse_from_matrix_market<Value_, Index_, ColumnIndex_>([&]() -> auto { return byteme::ZlibBufferReader(contents, length, 3, buffer_size); }, chunk_size);
        }
#else
        if (compression != -1) {
            throw std::runtime_error("layered sparse matrix reader not compiled with support for non-zero 'compression'");
        }
#endif
    }
    return read_layered_sparse_from_matrix_market<Value_, Index_, ColumnIndex_>([&]() -> auto { return byteme::RawBufferReader(contents, length); }, chunk_size);

}

}

#endif
