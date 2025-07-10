//============================================================================
// Name         : test_matrix.cpp
// Author       : Dale Roberts <dale.o.roberts@gmail.com>
// Contributors : 
// Copyright    : Copyright 2017-2025 Geoscience Australia
//
//                Licensed under the Apache License, Version 2.0 (the "License");
//                you may not use this file except in compliance with the License.
//                You may obtain a copy of the License at
//               
//                http ://www.apache.org/licenses/LICENSE-2.0
//               
//                Unless required by applicable law or agreed to in writing, software
//                distributed under the License is distributed on an "AS IS" BASIS,
//                WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//                See the License for the specific language governing permissions and
//                limitations under the License.
//
// Description  : Unit tests
//============================================================================

#define TESTING_MAIN

#include "math/dnamatrix_contiguous.hpp"
#include "testing.hpp"

using namespace dynadjust::math;

TEST_CASE("Constructor initializes matrix correctly", "[matrix_2d]") {
    matrix_2d mat(3, 4);

    REQUIRE(mat.rows() == 3);
    REQUIRE(mat.columns() == 4);
}

TEST_CASE("Copy constructor", "[matrix_2d]") {
    matrix_2d mat1(2, 2);
    mat1.put(0, 0, 1.0);
    mat1.put(0, 1, 2.0);
    mat1.put(1, 0, 3.0);
    mat1.put(1, 1, 4.0);

    matrix_2d mat2(mat1);

    REQUIRE(mat2.rows() == 2);
    REQUIRE(mat2.columns() == 2);
    REQUIRE(mat2.get(0, 0) == 1.0);
    REQUIRE(mat2.get(0, 1) == 2.0);
    REQUIRE(mat2.get(1, 0) == 3.0);
    REQUIRE(mat2.get(1, 1) == 4.0);
}

TEST_CASE("Addition", "[matrix_2d]") {
    matrix_2d mat1(2, 2);
    mat1.put(0, 0, 1.0);
    mat1.put(0, 1, 2.0);
    mat1.put(1, 0, 3.0);
    mat1.put(1, 1, 4.0);

    matrix_2d mat2(2, 2);
    mat2.put(0, 0, 5.0);
    mat2.put(0, 1, 6.0);
    mat2.put(1, 0, 7.0);
    mat2.put(1, 1, 8.0);

    matrix_2d result = mat1.add(mat2);

    REQUIRE(result.rows() == 2);
    REQUIRE(result.columns() == 2);
    REQUIRE(result.get(0, 0) == 6.0);
    REQUIRE(result.get(0, 1) == 8.0);
    REQUIRE(result.get(1, 0) == 10.0);
    REQUIRE(result.get(1, 1) == 12.0);
}

TEST_CASE("Square matrix multiplication", "[matrix_2d]") {
    matrix_2d mat1(2, 2);
    mat1.put(0, 0, 1.0);
    mat1.put(0, 1, 2.0);
    mat1.put(1, 0, 3.0);
    mat1.put(1, 1, 4.0);

    matrix_2d mat2(2, 2);
    mat2.put(0, 0, 1.0);
    mat2.put(0, 1, 2.0);
    mat2.put(1, 0, 3.0);
    mat2.put(1, 1, 4.0);

    matrix_2d result = mat1.multiply("N", mat2, "N");

    REQUIRE(result.rows() == 2);
    REQUIRE(result.columns() == 2);
    REQUIRE(result.get(0, 0) == 7.0);
    REQUIRE(result.get(0, 1) == 10.0);
    REQUIRE(result.get(1, 0) == 15.0);
    REQUIRE(result.get(1, 1) == 22.0);
}

TEST_CASE("Rectanglar matrix multiplication", "[matrix_2d]") {
    matrix_2d mat1(2, 3);
    mat1.put(0, 0, 1.0);
    mat1.put(0, 1, 2.0);
    mat1.put(0, 2, 3.0);
    mat1.put(1, 0, 4.0);
    mat1.put(1, 1, 5.0);
    mat1.put(1, 2, 6.0);

    matrix_2d mat2(3, 2);
    mat2.put(0, 0, 7.0);
    mat2.put(0, 1, 8.0);
    mat2.put(1, 0, 9.0);
    mat2.put(1, 1, 10.0);
    mat2.put(2, 0, 11.0);
    mat2.put(2, 1, 12.0);

    matrix_2d result = mat1.multiply("N", mat2, "N");

    REQUIRE(result.rows() == 2);
    REQUIRE(result.columns() == 2);
    REQUIRE(result.get(0, 0) == 58.0);
    REQUIRE(result.get(0, 1) == 64.0);
    REQUIRE(result.get(1, 0) == 139.0);
    REQUIRE(result.get(1, 1) == 154.0);
}

TEST_CASE("Transpose", "[matrix_2d]") {
    matrix_2d mat(2, 3);
    mat.put(0, 0, 1.0);
    mat.put(0, 1, 2.0);
    mat.put(0, 2, 3.0);
    mat.put(1, 0, 4.0);
    mat.put(1, 1, 5.0);
    mat.put(1, 2, 6.0);

    matrix_2d transposed = mat.transpose();

    REQUIRE(transposed.rows() == 3);
    REQUIRE(transposed.columns() == 2);
    REQUIRE(transposed.get(0, 0) == 1.0);
    REQUIRE(transposed.get(0, 1) == 4.0);
    REQUIRE(transposed.get(1, 0) == 2.0);
    REQUIRE(transposed.get(1, 1) == 5.0);
    REQUIRE(transposed.get(2, 0) == 3.0);
    REQUIRE(transposed.get(2, 1) == 6.0);
}

TEST_CASE("Scale", "[matrix_2d]") {
    matrix_2d mat(2, 2);
    mat.put(0, 0, 1.0);
    mat.put(0, 1, 2.0);
    mat.put(1, 0, 3.0);
    mat.put(1, 1, 4.0);

    matrix_2d scaled = mat.scale(2.0);

    REQUIRE(scaled.rows() == 2);
    REQUIRE(scaled.columns() == 2);
    REQUIRE(scaled.get(0, 0) == 2.0);
    REQUIRE(scaled.get(0, 1) == 4.0);
    REQUIRE(scaled.get(1, 0) == 6.0);
    REQUIRE(scaled.get(1, 1) == 8.0);
}

TEST_CASE("Sweep inverse", "[matrix_2d]") {
    matrix_2d mat(3, 3);
    mat.put(0, 0, 4.0);
    mat.put(0, 1, -1.0);
    mat.put(0, 2, -1.0);
    mat.put(1, 0, -1.0);
    mat.put(1, 1, 3.0);
    mat.put(1, 2, -1.0);
    mat.put(2, 0, -1.0);
    mat.put(2, 1, -1.0);
    mat.put(2, 2, 2.0);

    matrix_2d inverse = mat.sweepinverse();

    REQUIRE(inverse.rows() == 3);
    REQUIRE(inverse.columns() == 3);
    REQUIRE(abs(inverse.get(0, 0) - 0.384615) < 0.0001);
    REQUIRE(abs(inverse.get(0, 1) - 0.230769) < 0.0001);
    REQUIRE(abs(inverse.get(0, 2) - 0.307692) < 0.0001);
    REQUIRE(abs(inverse.get(1, 0) - 0.230769) < 0.0001);
    REQUIRE(abs(inverse.get(1, 1) - 0.538462) < 0.0001);
    REQUIRE(abs(inverse.get(1, 2) - 0.384615) < 0.0001);
    REQUIRE(abs(inverse.get(2, 0) - 0.307692) < 0.0001);
    REQUIRE(abs(inverse.get(2, 1) - 0.384615) < 0.0001);
    REQUIRE(abs(inverse.get(2, 2) - 0.846154) < 0.0001);
}

TEST_CASE("Cholesky inverse", "[matrix_2d]") {
    matrix_2d mat(3, 3);
    mat.put(0, 0, 4.0);
    mat.put(0, 1, -1.0);
    mat.put(0, 2, -1.0);
    mat.put(1, 0, -1.0);
    mat.put(1, 1, 3.0);
    mat.put(1, 2, -1.0);
    mat.put(2, 0, -1.0);
    mat.put(2, 1, -1.0);
    mat.put(2, 2, 2.0);

    matrix_2d inverse = mat.cholesky_inverse();

    REQUIRE(inverse.rows() == 3);
    REQUIRE(inverse.columns() == 3);
    REQUIRE(abs(inverse.get(0, 0) - 0.384615) < 0.0001);
    REQUIRE(abs(inverse.get(0, 1) - 0.230769) < 0.0001);
    REQUIRE(abs(inverse.get(0, 2) - 0.307692) < 0.0001);
    REQUIRE(abs(inverse.get(1, 0) - 0.230769) < 0.0001);
    REQUIRE(abs(inverse.get(1, 1) - 0.538462) < 0.0001);
    REQUIRE(abs(inverse.get(1, 2) - 0.384615) < 0.0001);
    REQUIRE(abs(inverse.get(2, 0) - 0.307692) < 0.0001);
    REQUIRE(abs(inverse.get(2, 1) - 0.384615) < 0.0001);
    REQUIRE(abs(inverse.get(2, 2) - 0.846154) < 0.0001);
}

TEST_CASE("Submatrix", "[matrix_2d]") {
    matrix_2d mat(4, 4);
    for (UINT32 i = 0; i < 4; ++i) {
        for (UINT32 j = 0; j < 4; ++j) {
            mat.put(i, j, static_cast<double>(i * 4 + j));
        }
    }

    matrix_2d submat = mat.submatrix(1, 1, 2, 2);

    REQUIRE(submat.rows() == 2);
    REQUIRE(submat.columns() == 2);
    REQUIRE(submat.get(0, 0) == 5.0);
    REQUIRE(submat.get(0, 1) == 6.0);
    REQUIRE(submat.get(1, 0) == 9.0);
    REQUIRE(submat.get(1, 1) == 10.0);
}

TEST_CASE("Element retrieval", "[matrix_2d]") {
    matrix_2d mat(3, 3);
    for (UINT32 i = 0; i < 3; ++i) {
        for (UINT32 j = 0; j < 3; ++j) {
            mat.put(i, j, static_cast<double>(i * 3 + j));
        }
    }

    REQUIRE(mat.get(0, 0) == 0.0);
    REQUIRE(mat.get(0, 1) == 1.0);
    REQUIRE(mat.get(0, 2) == 2.0);
    REQUIRE(mat.get(1, 0) == 3.0);
    REQUIRE(mat.get(1, 1) == 4.0);
    REQUIRE(mat.get(1, 2) == 5.0);
    REQUIRE(mat.get(2, 0) == 6.0);
    REQUIRE(mat.get(2, 1) == 7.0);
    REQUIRE(mat.get(2, 2) == 8.0);
}

TEST_CASE("Element modification", "[matrix_2d]") {
    matrix_2d mat(3, 3);
    for (UINT32 i = 0; i < 3; ++i) {
        for (UINT32 j = 0; j < 3; ++j) {
            mat.put(i, j, static_cast<double>(i * 3 + j));
        }
    }

    mat.put(1, 1, 99.0);

    REQUIRE(mat.get(1, 1) == 99.0);
}

TEST_CASE("Matrix allocation", "[matrix_2d]") {
    matrix_2d mat;
    mat.allocate(3, 4); // Does not change dimensions?

    REQUIRE(mat.rows() == 0);
    REQUIRE(mat.columns() == 0);
}

TEST_CASE("Matrix redimensioning", "[matrix_2d]") {
    matrix_2d mat(3, 4);
    mat.redim(5, 6);

    REQUIRE(mat.rows() == 5);
    REQUIRE(mat.columns() == 6);
}

TEST_CASE("Matrix shrinking", "[matrix_2d]") {
    matrix_2d mat(5, 6);
    mat.shrink(3, 4);

    REQUIRE(mat.rows() == 2);
    REQUIRE(mat.columns() == 2);
}

TEST_CASE("Clear lower triangle", "[matrix_2d]") {
    matrix_2d mat(3, 3);
    for (UINT32 i = 0; i < 3; ++i) {
        for (UINT32 j = 0; j < 3; ++j) {
            mat.put(i, j, static_cast<double>(i * 3 + j));
        }
    }

    mat.clearlower();

    REQUIRE(mat.get(1, 0) == 0.0);
    REQUIRE(mat.get(2, 0) == 0.0);
    REQUIRE(mat.get(2, 1) == 0.0);
}

TEST_CASE("Fill lower triangle", "[matrix_2d]") {
    matrix_2d mat(3, 3);
    for (UINT32 i = 0; i < 3; ++i) {
        for (UINT32 j = 0; j < 3; ++j) {
            mat.put(i, j, static_cast<double>(i * 3 + j));
        }
    }

    mat.filllower();

    REQUIRE(mat.get(1, 0) == mat.get(0, 1));
    REQUIRE(mat.get(2, 0) == mat.get(0, 2));
    REQUIRE(mat.get(2, 1) == mat.get(1, 2));
}

TEST_CASE("Fill upper triangle", "[matrix_2d]") {
    matrix_2d mat(3, 3);
    for (UINT32 i = 0; i < 3; ++i) {
        for (UINT32 j = 0; j < 3; ++j) {
            mat.put(i, j, static_cast<double>(i * 3 + j));
        }
    }

    mat.fillupper();

    REQUIRE(mat.get(0, 1) == mat.get(1, 0));
    REQUIRE(mat.get(0, 2) == mat.get(2, 0));
    REQUIRE(mat.get(1, 2) == mat.get(2, 1));
}

TEST_CASE("Zeroing matrix", "[matrix_2d]") {
    matrix_2d mat(3, 3);
    for (UINT32 i = 0; i < 3; ++i) {
        for (UINT32 j = 0; j < 3; ++j) {
            mat.put(i, j, static_cast<double>(i * 3 + j));
        }
    }

    mat.zero();

    REQUIRE(mat.get(0, 0) == 0.0);
    REQUIRE(mat.get(0, 1) == 0.0);
    REQUIRE(mat.get(0, 2) == 0.0);
    REQUIRE(mat.get(1, 0) == 0.0);
    REQUIRE(mat.get(1, 1) == 0.0);
    REQUIRE(mat.get(1, 2) == 0.0);
    REQUIRE(mat.get(2, 0) == 0.0);
    REQUIRE(mat.get(2, 1) == 0.0);
    REQUIRE(mat.get(2, 2) == 0.0);
}

TEST_CASE("Zeroing submatrix", "[matrix_2d]") {
    matrix_2d mat(4, 4);
    for (UINT32 i = 0; i < 4; ++i) {
        for (UINT32 j = 0; j < 4; ++j) {
            mat.put(i, j, static_cast<double>(i * 4 + j));
        }
    }

    mat.zero(1, 1, 2, 2);

    REQUIRE(mat.get(1, 1) == 0.0);
    REQUIRE(mat.get(1, 2) == 0.0);
    REQUIRE(mat.get(2, 1) == 0.0);
    REQUIRE(mat.get(2, 2) == 0.0);
}

TEST_CASE("Compute maximum value", "[matrix_2d]") {
    matrix_2d mat(3, 3);
    for (UINT32 i = 0; i < 3; ++i) {
        for (UINT32 j = 0; j < 3; ++j) {
            mat.put(i, j, static_cast<double>(i * 3 + j));
        }
    }

    REQUIRE(mat.compute_maximum_value() == 8.0);
}
