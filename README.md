# tatami helpers for creating layered matrices

![Unit tests](https://github.com/tatami-inc/tatami_layered/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/tatami-inc/tatami_layered/actions/workflows/doxygenate.yaml/badge.svg)
[![codecov](https://codecov.io/gh/tatami-inc/tatami_layered/graph/badge.svg?token=I5w28POARD)](https://codecov.io/gh/tatami-inc/tatami_layered)

## Overview

Layered matrices are a space optimization of sparse matrices containing small positive counts,
where we store rows in different "layers" depending on whether their maximum count is large enough to fit into an unsigned 8-bit, 16-bit or 32-bit integer.
This reduces the memory usage compared to naively storing counts for all rows in the largest integer size across the entire matrix.
It is intended to be used with gene expression data where different genes (rows) can vary widely in their expression.

## Quick start

We can easily convert an existing `tatami::Matrix` to a layered sparse matrix:

```cpp
#include "tatami_layered/tatami_layered.hpp"

auto converted = tatami_layered::convert_to_layered_sparse(*mat);
```

We can also read a layered sparse matrix from a Matrix Market file:

```cpp
auto loaded = tatami_layered::read_layered_sparse_from_matrix_market_text_file(path.c_str());
```

Check out the [documentation](https://tatami-inc.github.io/tatami_layered) for more details.

## Building projects

### CMake with `FetchContent`

If you're using CMake, you just need to add something like this to your `CMakeLists.txt`:

```cmake
include(FetchContent)

FetchContent_Declare(
  tatami_layered
  GIT_REPOSITORY https://github.com/tatami-inc/tatami_layered
  GIT_TAG master # or any version of interest 
)

FetchContent_MakeAvailable(tatami_layered)
```

Then you can link to **tatami_layered** to make the headers available during compilation:

```cmake
# For executables:
target_link_libraries(myexe tatami_layered)

# For libaries
target_link_libraries(mylib INTERFACE tatami_layered)
```

### CMake using `find_package()`

You can install the library by cloning a suitable version of this repository and running the following commands:

```sh
mkdir build && cd build
cmake .. -DTATAMI_LAYERED_TESTS=OFF
cmake --build . --target install
```

Then you can use `find_package()` as usual:

```cmake
find_package(tatami_tatami_layered CONFIG REQUIRED)
target_link_libraries(mylib INTERFACE tatami::tatami_layered)
```

### Manual

If you're not using CMake, the simple approach is to just copy the files - either directly or with Git submodules - and include their path during compilation with, e.g., GCC's `-I`.
This will also require the various dependencies listed in the [`extern/CMakeLists.txt`](extern/CMakeLists.txt) file.
