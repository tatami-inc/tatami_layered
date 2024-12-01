#include <gtest/gtest.h>

#include "tatami/tatami.hpp"
#include "tatami_test/tatami_test.hpp"
#include "tatami_layered/read_layered_sparse_from_matrix_market.hpp"

#include "temp_file_path.h"
#include "mock_layered_sparse_data.h"

#include <fstream>
#include <sstream>

class ReadLayeredSparseFromMatrixMarketTest : public ::testing::TestWithParam<std::tuple<bool, int> > {
protected:
    size_t NR = 2000, NC = 1000;

    template<class Stream, class U, class V, class W>
    static void write_matrix_market(Stream& stream, size_t nr, size_t nc, const U& vals, const V& rows, const W& cols, bool scramble) {
        stream << "%%MatrixMarket matrix coordinate integer general\n";
        stream << nr << " " << nc << " " << vals.size();

        if (scramble) {
            std::vector<size_t> reorder(vals.size());
            std::iota(reorder.begin(), reorder.end(), 0);
            std::mt19937_64 rng(vals.size() - nr + nc);
            std::shuffle(reorder.begin(), reorder.end(), rng);
            for (auto i : reorder) {
                stream << "\n" << rows[i] + 1 << " " << cols[i] + 1 << " " << vals[i];
            }

        } else {
            for (size_t i = 0; i < vals.size(); ++i) {
                stream << "\n" << rows[i] + 1 << " " << cols[i] + 1 << " " << vals[i];
            }
        }
        stream << "\n";

        return;
    }
};

TEST_P(ReadLayeredSparseFromMatrixMarketTest, File) {
    auto param = GetParam();
    auto scrambled = std::get<0>(param);
    auto compressed = std::get<1>(param);

    std::vector<size_t> rows, cols;
    std::vector<int> vals;
    mock_layered_sparse_data<false>(NR, NC, rows, cols, vals);

    auto indptrs = tatami::compress_sparse_triplets<false>(NR, NC, vals, rows, cols);
    typedef tatami::CompressedSparseColumnMatrix<double, int, decltype(vals), decltype(rows), decltype(indptrs)> SparseMat; 
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new SparseMat(NR, NC, vals, rows, indptrs)); 

    auto path = temp_file_path("tatami-tests-ext-MatrixMarket");
    if (!compressed) {
        std::ofstream file_out(path);
        write_matrix_market(file_out, NR, NC, vals, rows, cols, scrambled);
    } else {
        std::stringstream sstream;
        write_matrix_market(sstream, NR, NC, vals, rows, cols, scrambled);
        auto contents = sstream.str();

        path += ".gz";
        gzFile ohandle = gzopen(path.c_str(), "w");
        gzwrite(ohandle, contents.data(), contents.size());
        gzclose(ohandle);
    }

    std::shared_ptr<tatami::NumericMatrix> out;
    if (compressed == 0) {
        out = tatami_layered::read_layered_sparse_from_matrix_market_text_file(path.c_str());
    } else if (compressed == 1) {
        out = tatami_layered::read_layered_sparse_from_matrix_market_gzip_file(path.c_str());
    } else if (compressed == 2) {
        out = tatami_layered::read_layered_sparse_from_matrix_market_some_file(path.c_str());
    }

    tatami_test::test_simple_row_access(*out, *ref);
    tatami_test::test_simple_column_access(*out, *ref);
}

TEST_P(ReadLayeredSparseFromMatrixMarketTest, Buffer) {
    auto param = GetParam();
    auto scrambled = std::get<0>(param);
    auto compressed = std::get<1>(param);

    std::vector<size_t> rows, cols;
    std::vector<int> vals;
    mock_layered_sparse_data<false>(NR, NC, rows, cols, vals);

    auto indptrs = tatami::compress_sparse_triplets<false>(NR, NC, vals, rows, cols);
    typedef tatami::CompressedSparseColumnMatrix<double, int, decltype(vals), decltype(rows), decltype(indptrs)> SparseMat; 
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new SparseMat(NR, NC, vals, rows, indptrs)); 

    std::stringstream buf_out;
    write_matrix_market(buf_out, NR, NC, vals, rows, cols, scrambled);
    auto contents = buf_out.str();

    if (compressed) {
        auto path = temp_file_path("tatami-tests-ext-MatrixMarket-gzip");
        gzFile ohandle = gzopen(path.c_str(), "w");
        gzwrite(ohandle, contents.data(), contents.size());
        gzclose(ohandle);

        std::ifstream handle(path, std::ios_base::binary);
        handle >> std::noskipws;
        contents.clear();
        contents.insert(contents.end(), std::istream_iterator<unsigned char>{handle}, std::istream_iterator<unsigned char>());
    }

    std::shared_ptr<tatami::NumericMatrix> out;
    auto ptr = reinterpret_cast<const unsigned char*>(contents.data());
    size_t n = contents.size();

    if (compressed == 0) {
        out = tatami_layered::read_layered_sparse_from_matrix_market_text_buffer(ptr, n);
    } else if (compressed == 1) {
        out = tatami_layered::read_layered_sparse_from_matrix_market_zlib_buffer(ptr, n);
    } else if (compressed == 2) {
        out = tatami_layered::read_layered_sparse_from_matrix_market_some_buffer(ptr, n);
    }

    tatami_test::test_simple_row_access(*out, *ref);
    tatami_test::test_simple_column_access(*out, *ref);
}

INSTANTIATE_TEST_SUITE_P(
    ReadLayeredSparseFromMatrixMarket,
    ReadLayeredSparseFromMatrixMarketTest,
    ::testing::Combine(
        ::testing::Values(false, true), // scrambled?
        ::testing::Values(0, 1, 2)  // Gzipped?
    )
);
