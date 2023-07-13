#ifndef TATAMI_LAYERED_READ_LAYERED_SPARSE_FROM_MATRIX_MARKET_HPP
#define TATAMI_LAYERED_READ_LAYERED_SPARSE_FROM_MATRIX_MARKET_HPP

#include "byteme/byteme.hpp"
#include "eminem/eminem.hpp"

namespace tatami_layered {

template<typename Value_, typename Index_, typename ColumnIndex_>
std::shared_ptr<tatami::Matrix<Value_, Index_> > read_layered_sparse_from_matrix_market_file(const char* path, size_t buffer_size, Index_ chunk_size) {
    eminem::SomeFileParser parser(path, buffer_size);


}

}

#endif
