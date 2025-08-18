#ifndef TATAMI_LAYERED_READ_LAYERED_SPARSE_FROM_MATRIX_MARKET_HPP
#define TATAMI_LAYERED_READ_LAYERED_SPARSE_FROM_MATRIX_MARKET_HPP

#include <vector>
#include <algorithm>
#include <cstddef>

#include "byteme/byteme.hpp"
#include "eminem/eminem.hpp"
#include "tatami/tatami.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils.hpp"

/**
 * @file read_layered_sparse_from_matrix_market.hpp
 * @brief Read layered sparse matrices from Matrix Market files.
 */

namespace tatami_layered {

/**
 * @cond
 */
template<typename Value_, typename Index_, typename ColumnIndex_, class Creator_>
std::shared_ptr<tatami::Matrix<Value_, Index_> > read_layered_sparse_from_matrix_market(Creator_ create, Index_ chunk_size, int num_threads) {
    Index_ NR, NC, nchunks, leftovers;

    std::vector<Holder< std::uint8_t, Index_, ColumnIndex_> > store8;
    std::vector<Holder<std::uint16_t, Index_, ColumnIndex_> > store16;
    std::vector<Holder<std::uint32_t, Index_, ColumnIndex_> > store32;

    std::vector<std::vector<Index_> > identities8, identities16, identities32;
    std::vector<std::vector<Index_> > assigned_position;
    std::vector<std::vector<Category> > assigned_category;

    eminem::ParserOptions eopt;
    eopt.num_threads = num_threads;

    // First pass, scanning for the max and number.
    {
        auto reader = create();
        byteme::PerByteSerial<char, byteme::Reader*> pb(&reader);
        eminem::Parser<decltype(&pb), Index_> parser(&pb, eopt);

        parser.scan_preamble();
        NR = parser.get_nrows();
        NC = parser.get_ncols();
        leftovers = NC % chunk_size;
        nchunks = atleastone(NC / chunk_size + (leftovers != 0));

        tatami::resize_container_to_Index_size(store8, nchunks);
        tatami::resize_container_to_Index_size(store16, nchunks);
        tatami::resize_container_to_Index_size(store32, nchunks);
        tatami::resize_container_to_Index_size(identities8, nchunks);
        tatami::resize_container_to_Index_size(identities16, nchunks);
        tatami::resize_container_to_Index_size(identities32, nchunks);
        tatami::resize_container_to_Index_size(assigned_position, nchunks);
        tatami::resize_container_to_Index_size(assigned_category, nchunks);

        auto max_per_chunk = tatami::create_container_of_Index_size<std::vector<std::vector<Category> > >(nchunks);
        tatami::cast_Index_to_container_size<std::vector<Category> >(NR); // check that resize is safe.
        for (auto& x : max_per_chunk) { x.resize(NR); }

        auto num_per_chunk = tatami::create_container_of_Index_size<std::vector<std::vector<Index_> > >(nchunks);
        tatami::cast_Index_to_container_size<std::vector<Index_> >(NR);
        for (auto& x : num_per_chunk) { x.resize(NR); }

        const auto& banner = parser.get_banner();
        if (banner.field == eminem::Field::INTEGER) {
            parser.template scan_integer<std::uint32_t>([&](Index_ r, Index_ c, std::uint32_t val) -> void {
                auto chunk = (c - 1) / chunk_size;
                --r;
                max_per_chunk[chunk][r] = std::max(max_per_chunk[chunk][r], categorize(val));
                ++num_per_chunk[chunk][r];
            });

        } else if (banner.field == eminem::Field::DOUBLE || banner.field == eminem::Field::REAL) {
            parser.scan_real([&](Index_ r, Index_ c, double val) -> void {
                auto chunk = (c - 1) / chunk_size;
                --r;
                max_per_chunk[chunk][r] = std::max(max_per_chunk[chunk][r], categorize(val));
                ++num_per_chunk[chunk][r];
            });

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
        std::vector<std::vector<std::size_t> > output_positions(nchunks);
        tatami::cast_Index_to_container_size<decltype(output_positions.front())>(NR); // check that resize is safe.
        for (decltype(nchunks) chunk = 0; chunk < nchunks; ++chunk) {
            output_positions[chunk].resize(NR);
            for (decltype(NR) r = 0; r < NR; ++r) {
                output_positions[chunk][r] = get_sparse_ptr(store8, store16, store32, assigned_category, assigned_position, chunk, r);
            }
        }

        auto reader = create();
        byteme::PerByteSerial<char, decltype(&reader)> pb(&reader);
        eminem::Parser<decltype(&pb), Index_> parser(&pb, eopt);

        parser.scan_preamble();
        const auto& banner = parser.get_banner();
        if (banner.field == eminem::Field::INTEGER) {
            parser.template scan_integer<std::uint32_t>([&](Index_ r, Index_ c, std::uint32_t val) -> void {
                --c;
                Index_ chunk = c / chunk_size;
                Index_ offset = c % chunk_size;
                --r;
                fill_sparse_value(store8, store16, store32, assigned_category[chunk][r], chunk, offset, val, output_positions[chunk][r]++);
            });

        } else if (banner.field == eminem::Field::DOUBLE || banner.field == eminem::Field::REAL) {
            parser.scan_real([&](Index_ r, Index_ c, double val) -> void {
                --c;
                Index_ chunk = c / chunk_size;
                Index_ offset = c % chunk_size;
                --r;
                fill_sparse_value(store8, store16, store32, assigned_category[chunk][r], chunk, offset, val, output_positions[chunk][r]++);
            });
        }

        // Checking that the column indices are sorted properly.
        auto sorter = [&](auto& store) -> void {
            std::vector<std::pair<typename decltype(store[0].index)::value_type, typename decltype(store[0].value)::value_type> > buffer;
            buffer.reserve(chunk_size);

            for (auto& st : store) {
                auto num_ptr = st.ptr.size();
                for (decltype(num_ptr) r = 1; r < num_ptr; ++r) {
                    auto start = st.ptr[r - 1], end = st.ptr[r];

                    if (!std::is_sorted(st.index.begin() + start, st.index.begin() + end)) {
                        buffer.clear();
                        for (auto i = start; i < end; ++i) {
                            buffer.emplace_back(st.index[i], st.value[i]);
                        }

                        std::sort(buffer.begin(), buffer.end());
                        auto bIt = buffer.begin();
                        for (auto i = start; i < end; ++i, ++bIt) {
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
        NR,
        chunk_size,
        leftovers
    );
}
/**
 * @endcond
 */

/**
 * Options for `read_layered_sparse_from_matrix_market_text_file()` and friends.
 */
struct ReadLayeredSparseFromMatrixMarketOptions {
    /**
     * Chunk size to use for partitioning columns.
     */
    std::size_t chunk_size = sanisizer::cap<std::size_t>(65536);

    /**
     * Size of the buffer (in bytes) to use when reading from file.
     */
    std::size_t buffer_size = sanisizer::cap<std::size_t>(65536);

    /**
     * Number of threads for Matrix Market parsing.
     */
    int num_threads = 1;
};

/**
 * @param filepath Path to an uncompressed Matrix Market text file.
 * @param options Further options.
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
template<typename Value_ = double, typename Index_ = int, typename ColumnIndex_ = std::uint16_t>
std::shared_ptr<tatami::Matrix<Value_, Index_> > read_layered_sparse_from_matrix_market_text_file(const char* filepath, const ReadLayeredSparseFromMatrixMarketOptions& options) {
    return read_layered_sparse_from_matrix_market<Value_, Index_, ColumnIndex_>(
        [&]() -> auto {
            return byteme::RawFileReader(filepath, [&]{
                byteme::RawFileReaderOptions opt;
                opt.buffer_size = options.buffer_size;
                return opt;
            }());
        },
        check_chunk_size<Index_, ColumnIndex_>(options.chunk_size),
        options.num_threads
    );
}

/**
 * @cond
 */
// Back-compatibility.
template<typename Value_ = double, typename Index_ = int, typename ColumnIndex_ = std::uint16_t>
std::shared_ptr<tatami::Matrix<Value_, Index_> > read_layered_sparse_from_matrix_market_text_file(const char* filepath, Index_ chunk_size = 65536, std::size_t buffer_size = 65536) {
    return read_layered_sparse_from_matrix_market_text_file<Value_, Index_, ColumnIndex_>(filepath, [&]{
        ReadLayeredSparseFromMatrixMarketOptions opt;
        opt.chunk_size = chunk_size;
        opt.buffer_size = buffer_size;
        return opt;
    }());
}
/**
 * @endcond
 */

#if __has_include("zlib.h")

/**
 * @param filepath Path to a (possibly Gzip-compressed) Matrix Market file.
 * @param chunk_size Chunk size to use for partitioning columns.
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
template<typename Value_ = double, typename Index_ = int, typename ColumnIndex_ = std::uint16_t>
std::shared_ptr<tatami::Matrix<Value_, Index_> > read_layered_sparse_from_matrix_market_some_file(const char* filepath, const ReadLayeredSparseFromMatrixMarketOptions& options) {
    return read_layered_sparse_from_matrix_market<Value_, Index_, ColumnIndex_>(
        [&]() -> auto {
            return byteme::SomeFileReader(filepath, [&]{
                byteme::SomeFileReaderOptions opt;
                opt.buffer_size = options.buffer_size;
                return opt;
            }());
        },
        check_chunk_size<Index_, ColumnIndex_>(options.chunk_size),
        options.num_threads
    );
}

/**
 * @param filepath Path to a Gzip-compressed Matrix Market file.
 * @param chunk_size Chunk size to use for partitioning columns.
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
template<typename Value_ = double, typename Index_ = int, typename ColumnIndex_ = std::uint16_t>
std::shared_ptr<tatami::Matrix<Value_, Index_> > read_layered_sparse_from_matrix_market_gzip_file(const char* filepath, const ReadLayeredSparseFromMatrixMarketOptions& options) {
    return read_layered_sparse_from_matrix_market<Value_, Index_, ColumnIndex_>(
        [&]() -> auto {
            return byteme::GzipFileReader(filepath, [&]{
                byteme::GzipFileReaderOptions opt;
                opt.buffer_size = options.buffer_size;
                return opt;
            }());
        },
        check_chunk_size<Index_, ColumnIndex_>(options.chunk_size),
        options.num_threads
    );
}

/**
 * @cond
 */
// Back-compatibility.
template<typename Value_ = double, typename Index_ = int, typename ColumnIndex_ = std::uint16_t>
std::shared_ptr<tatami::Matrix<Value_, Index_> > read_layered_sparse_from_matrix_market_some_file(const char* filepath, Index_ chunk_size = 65536, std::size_t buffer_size = 65536) {
    return read_layered_sparse_from_matrix_market_some_file<Value_, Index_, ColumnIndex_>(filepath, [&]{
        ReadLayeredSparseFromMatrixMarketOptions opt;
        opt.chunk_size = chunk_size;
        opt.buffer_size = buffer_size;
        return opt;
    }());
}

template<typename Value_ = double, typename Index_ = int, typename ColumnIndex_ = std::uint16_t>
std::shared_ptr<tatami::Matrix<Value_, Index_> > read_layered_sparse_from_matrix_market_gzip_file(const char* filepath, Index_ chunk_size = 65536, std::size_t buffer_size = 65536) {
    return read_layered_sparse_from_matrix_market_gzip_file<Value_, Index_, ColumnIndex_>(filepath, [&]{
        ReadLayeredSparseFromMatrixMarketOptions opt;
        opt.chunk_size = chunk_size;
        opt.buffer_size = buffer_size;
        return opt;
    }());
}
/**
 * @endcond
 */


#endif

/**
 * @param contents Array containing the contents of an uncompressed Matrix Market text file.
 * @param length Length of the array.
 * @param chunk_size Chunk size to use for partitioning columns.
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
template<typename Value_ = double, typename Index_ = int, typename ColumnIndex_ = std::uint16_t>
std::shared_ptr<tatami::Matrix<Value_, Index_> > read_layered_sparse_from_matrix_market_text_buffer(
    const unsigned char* contents,
    std::size_t length,
    const ReadLayeredSparseFromMatrixMarketOptions& options)
{
    return read_layered_sparse_from_matrix_market<Value_, Index_, ColumnIndex_>(
        [&]() -> auto {
            return byteme::RawBufferReader(contents, length);
        },
        check_chunk_size<Index_, ColumnIndex_>(options.chunk_size),
        options.num_threads
    );
}

/**
 * @cond
 */
// Back-compatibility.
template<typename Value_ = double, typename Index_ = int, typename ColumnIndex_ = std::uint16_t>
std::shared_ptr<tatami::Matrix<Value_, Index_> > read_layered_sparse_from_matrix_market_text_buffer(const unsigned char* contents, std::size_t length, Index_ chunk_size = 65536) {
    return read_layered_sparse_from_matrix_market_text_buffer<Value_, Index_, ColumnIndex_>(contents, length, [&]{
        ReadLayeredSparseFromMatrixMarketOptions opt;
        opt.chunk_size = chunk_size;
        return opt;
    }());
}
/**
 * @endcond
 */


#if __has_include("zlib.h")

/**
 * @param contents Array containing the contents of a (possibly Gzip/Zlib-compressed) Matrix Market file.
 * @param length Length of the array.
 * @param chunk_size Chunk size to use for partitioning columns.
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
template<typename Value_ = double, typename Index_ = int, typename ColumnIndex_ = std::uint16_t>
std::shared_ptr<tatami::Matrix<Value_, Index_> > read_layered_sparse_from_matrix_market_some_buffer(
    const unsigned char* contents,
    std::size_t length,
    const ReadLayeredSparseFromMatrixMarketOptions& options)
{
    return read_layered_sparse_from_matrix_market<Value_, Index_, ColumnIndex_>(
        [&]() -> auto {
            return byteme::SomeBufferReader(contents, length, [&]{
                byteme::SomeBufferReaderOptions opt;
                opt.buffer_size = options.buffer_size;
                return opt;
            }());
        },
        check_chunk_size<Index_, ColumnIndex_>(options.chunk_size),
        options.num_threads
    );
}

/**
 * @param contents Array containing the contents of a Gzip/Zlib-compressed Matrix Market file.
 * @param length Length of the array.
 * @param chunk_size Chunk size to use for partitioning columns.
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
template<typename Value_ = double, typename Index_ = int, typename ColumnIndex_ = std::uint16_t>
std::shared_ptr<tatami::Matrix<Value_, Index_> > read_layered_sparse_from_matrix_market_zlib_buffer(
    const unsigned char* contents,
    std::size_t length,
    const ReadLayeredSparseFromMatrixMarketOptions& options)
{
    return read_layered_sparse_from_matrix_market<Value_, Index_, ColumnIndex_>(
        [&]() -> auto {
            return byteme::ZlibBufferReader(contents, length, [&]{
                byteme::ZlibBufferReaderOptions opt;
                opt.buffer_size = options.buffer_size;
                return opt;
            }());
        },
        check_chunk_size<Index_, ColumnIndex_>(options.chunk_size),
        options.num_threads
    );
}

/**
 * @cond
 */
// Back-compatibility.
template<typename Value_ = double, typename Index_ = int, typename ColumnIndex_ = std::uint16_t>
std::shared_ptr<tatami::Matrix<Value_, Index_> > read_layered_sparse_from_matrix_market_some_buffer(
    const unsigned char* contents,
    std::size_t length,
    Index_ chunk_size = 65536,
    std::size_t buffer_size = 65536)
{
    return read_layered_sparse_from_matrix_market_some_buffer<Value_, Index_, ColumnIndex_>(contents, length, [&]{
        ReadLayeredSparseFromMatrixMarketOptions opt;
        opt.chunk_size = chunk_size;
        opt.buffer_size = buffer_size;
        return opt;
    }());
}

template<typename Value_ = double, typename Index_ = int, typename ColumnIndex_ = std::uint16_t>
std::shared_ptr<tatami::Matrix<Value_, Index_> > read_layered_sparse_from_matrix_market_zlib_buffer(
    const unsigned char* contents,
    std::size_t length,
    Index_ chunk_size = 65536,
    std::size_t buffer_size = 65536)
{
    return read_layered_sparse_from_matrix_market_zlib_buffer<Value_, Index_, ColumnIndex_>(contents, length, [&]{
        ReadLayeredSparseFromMatrixMarketOptions opt;
        opt.chunk_size = chunk_size;
        opt.buffer_size = buffer_size;
        return opt;
    }());
}
/**
 * @endcond
 */


#endif

}

#endif
