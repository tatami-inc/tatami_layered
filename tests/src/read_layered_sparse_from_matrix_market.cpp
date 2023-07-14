#include <gtest/gtest.h>

#include "tatami/tatami.hpp"
#include "tatami_layered/read_layered_sparse_from_matrix_market.hpp"
#include "byteme/temp_file_path.hpp"

#include "mock_layered_sparse_data.h"

#include <fstream>
#include <sstream>

class ReadLayeredSparseFromMatrixMarket : public ::testing::Test {
protected:
    size_t NR = 2000, NC = 1000;

    template<class Stream, class U, class V, class W>
    static void write_matrix_market(Stream& stream, size_t nr, size_t nc, const U& vals, const V& rows, const W& cols) {
        stream << "%%MatrixMarket matrix coordinate integer general\n";
        stream << nr << " " << nc << " " << vals.size();

        for (size_t i = 0; i < vals.size(); ++i) {
            stream << "\n" << rows[i] + 1 << " " << cols[i] + 1 << " " << vals[i];
        }
        stream << "\n";

        return;
    }
};

TEST_F(ReadLayeredSparseFromMatrixMarket, File) {
    std::vector<size_t> rows, cols;
    std::vector<int> vals;
    mock_layered_sparse_data<false>(NR, NC, rows, cols, vals);

    auto indptrs = tatami::compress_sparse_triplets<false>(NR, NC, vals, rows, cols);
    typedef tatami::CompressedSparseColumnMatrix<double, int, decltype(vals), decltype(rows), decltype(indptrs)> SparseMat; 
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new SparseMat(NR, NC, std::move(vals), std::move(rows), std::move(indptrs))); 

    auto path = byteme::temp_file_path("tatami-tests-ext-MatrixMarket", ".mtx");
    {
        std::ofstream file_out(path);
        write_matrix_market(file_out, NR, NC, vals, rows, cols);
    }

    auto out = tatami_layered::read_layered_sparse_from_matrix_market_file(path.c_str());

    auto owrk = out->dense_row();
    auto rwrk = ref->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto stuff = owrk->fetch(i);
        EXPECT_EQ(stuff, rwrk->fetch(i));
    }
}

TEST_F(ReadLayeredSparseFromMatrixMarket, Buffer) {
    std::vector<size_t> rows, cols;
    std::vector<int> vals;
    mock_layered_sparse_data<false>(NR, NC, rows, cols, vals);

    auto indptrs = tatami::compress_sparse_triplets<false>(NR, NC, vals, rows, cols);
    typedef tatami::CompressedSparseColumnMatrix<double, int, decltype(vals), decltype(rows), decltype(indptrs)> SparseMat; 
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new SparseMat(NR, NC, std::move(vals), std::move(rows), std::move(indptrs))); 

    std::stringstream buf_out;
    write_matrix_market(buf_out, NR, NC, vals, rows, cols);
    auto contents = buf_out.str();

    auto out = tatami_layered::read_layered_sparse_from_matrix_market_file(contents.data(), contents.size());

    auto owrk = out->dense_row();
    auto rwrk = ref->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto stuff = owrk->fetch(i);
        EXPECT_EQ(stuff, rwrk->fetch(i));
    }
}

