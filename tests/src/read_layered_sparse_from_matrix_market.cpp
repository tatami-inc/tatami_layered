#include <gtest/gtest.h>

#include "tatami/tatami.hpp"
#include "tatami_layered/read_layered_sparse_from_matrix_market.hpp"
#include "byteme/temp_file_path.hpp"

#include "mock_layered_sparse_data.h"

#include <fstream>
#include <sstream>

class ReadLayeredSparseFromMatrixMarketTest : public ::testing::TestWithParam<bool> {
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
    auto scrambled = GetParam();

    std::vector<size_t> rows, cols;
    std::vector<int> vals;
    mock_layered_sparse_data<false>(NR, NC, rows, cols, vals);

    auto indptrs = tatami::compress_sparse_triplets<false>(NR, NC, vals, rows, cols);
    typedef tatami::CompressedSparseColumnMatrix<double, int, decltype(vals), decltype(rows), decltype(indptrs)> SparseMat; 
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new SparseMat(NR, NC, vals, rows, indptrs)); 

    auto path = byteme::temp_file_path("tatami-tests-ext-MatrixMarket", ".mtx");
    {
        std::ofstream file_out(path);
        write_matrix_market(file_out, NR, NC, vals, rows, cols, scrambled);
    }

    auto out = tatami_layered::read_layered_sparse_from_matrix_market_file(path.c_str());

    auto owrk = out->dense_row();
    auto rwrk = ref->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto stuff = owrk->fetch(i);
        EXPECT_EQ(stuff, rwrk->fetch(i));
    }
}

TEST_P(ReadLayeredSparseFromMatrixMarketTest, Buffer) {
    auto scrambled = GetParam();

    std::vector<size_t> rows, cols;
    std::vector<int> vals;
    mock_layered_sparse_data<false>(NR, NC, rows, cols, vals);

    auto indptrs = tatami::compress_sparse_triplets<false>(NR, NC, vals, rows, cols);
    typedef tatami::CompressedSparseColumnMatrix<double, int, decltype(vals), decltype(rows), decltype(indptrs)> SparseMat; 
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new SparseMat(NR, NC, vals, rows, indptrs)); 

    std::stringstream buf_out;
    write_matrix_market(buf_out, NR, NC, vals, rows, cols, scrambled);
    auto contents = buf_out.str();

    auto out = tatami_layered::read_layered_sparse_from_matrix_market_buffer(reinterpret_cast<const unsigned char*>(contents.data()), contents.size());

    auto owrk = out->dense_row();
    auto rwrk = ref->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto stuff = owrk->fetch(i);
        EXPECT_EQ(stuff, rwrk->fetch(i));
    }
}

INSTANTIATE_TEST_SUITE_P(
    ReadLayeredSparseFromMatrixMarket,
    ReadLayeredSparseFromMatrixMarketTest,
    ::testing::Values(false, true)
);
