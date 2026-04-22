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

#include <cmath>
#include <sstream>
#include <stdexcept>

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

    matrix_2d result(2, 2);
    result.multiply(mat1, "N", mat2, "N");

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

// ============================================================================
// Bug fix: operator= used || instead of && when checking if rhs fits
// ============================================================================

TEST_CASE("Assignment operator with smaller rhs reuses buffer", "[matrix_2d]") {
    matrix_2d big(4, 4);
    for (UINT32 i = 0; i < 4; ++i)
        for (UINT32 j = 0; j < 4; ++j)
            big.put(i, j, static_cast<double>(i * 4 + j));

    matrix_2d small(2, 2);
    small.put(0, 0, 10.0);
    small.put(0, 1, 20.0);
    small.put(1, 0, 30.0);
    small.put(1, 1, 40.0);

    big = small;

    REQUIRE(big.rows() == 2);
    REQUIRE(big.columns() == 2);
    REQUIRE(big.get(0, 0) == 10.0);
    REQUIRE(big.get(0, 1) == 20.0);
    REQUIRE(big.get(1, 0) == 30.0);
    REQUIRE(big.get(1, 1) == 40.0);
    // mem dimensions should be preserved (buffer reused)
    REQUIRE(big.memRows() == 4);
    REQUIRE(big.memColumns() == 4);
}

TEST_CASE("Assignment operator with larger rhs reallocates", "[matrix_2d]") {
    matrix_2d small(2, 2);
    small.put(0, 0, 1.0);
    small.put(0, 1, 2.0);
    small.put(1, 0, 3.0);
    small.put(1, 1, 4.0);

    matrix_2d big(4, 4);
    for (UINT32 i = 0; i < 4; ++i)
        for (UINT32 j = 0; j < 4; ++j)
            big.put(i, j, static_cast<double>(i * 4 + j));

    small = big;

    REQUIRE(small.rows() == 4);
    REQUIRE(small.columns() == 4);
    for (UINT32 i = 0; i < 4; ++i)
        for (UINT32 j = 0; j < 4; ++j)
            REQUIRE(small.get(i, j) == static_cast<double>(i * 4 + j));
}

TEST_CASE("Assignment operator with mismatched dimensions reallocates", "[matrix_2d]") {
    // The old bug: 2x4 assigned from 3x3 — old code used || so the
    // 2 >= 3 check (false) || 4 >= 3 check (true) would take the
    // reuse path even though rows didn't fit. With && it correctly
    // reallocates.
    matrix_2d lhs(2, 4);
    lhs.put(0, 0, 99.0);

    matrix_2d rhs(3, 3);
    for (UINT32 i = 0; i < 3; ++i)
        for (UINT32 j = 0; j < 3; ++j)
            rhs.put(i, j, static_cast<double>(i * 3 + j + 1));

    lhs = rhs;

    REQUIRE(lhs.rows() == 3);
    REQUIRE(lhs.columns() == 3);
    for (UINT32 i = 0; i < 3; ++i)
        for (UINT32 j = 0; j < 3; ++j)
            REQUIRE(lhs.get(i, j) == static_cast<double>(i * 3 + j + 1));
}

// ============================================================================
// Bug fix: copybuffer fast-path now checks oldmat dimensions too
// ============================================================================

TEST_CASE("copybuffer with different mem dimensions uses per-column copy", "[matrix_2d]") {
    // Create a 4x4 matrix, shrink to 2x2 (mem stays 4x4), then assign
    // from a native 2x2 (mem is 2x2). The copybuffer fast-path must not
    // memcpy the full 4x4 buffer from a 2x2 source.
    matrix_2d big(4, 4);
    for (UINT32 i = 0; i < 4; ++i)
        for (UINT32 j = 0; j < 4; ++j)
            big.put(i, j, 1.0);

    big.shrink(2, 2);  // now rows=2, cols=2, but mem_rows=4, mem_cols=4

    matrix_2d small(2, 2);
    small.put(0, 0, 10.0);
    small.put(0, 1, 20.0);
    small.put(1, 0, 30.0);
    small.put(1, 1, 40.0);

    big = small;

    REQUIRE(big.rows() == 2);
    REQUIRE(big.columns() == 2);
    REQUIRE(big.get(0, 0) == 10.0);
    REQUIRE(big.get(0, 1) == 20.0);
    REQUIRE(big.get(1, 0) == 30.0);
    REQUIRE(big.get(1, 1) == 40.0);
}

// ============================================================================
// Bug fix: LAPACK LDA must be _mem_rows, not _rows
// ============================================================================

TEST_CASE("Cholesky inverse with mem_rows != rows", "[matrix_2d]") {
    // Allocate a 4x4 matrix, shrink logical size to 3x3, then do
    // cholesky_inverse. This exercises the LDA fix (lda=4, n=3).
    matrix_2d mat(4, 4);

    // Fill the 3x3 leading block with a symmetric positive definite matrix
    mat.put(0, 0, 4.0);
    mat.put(0, 1, -1.0);
    mat.put(0, 2, -1.0);
    mat.put(1, 0, -1.0);
    mat.put(1, 1, 3.0);
    mat.put(1, 2, -1.0);
    mat.put(2, 0, -1.0);
    mat.put(2, 1, -1.0);
    mat.put(2, 2, 2.0);

    mat.shrink(1, 1);  // logical 3x3, mem 4x4

    REQUIRE(mat.rows() == 3);
    REQUIRE(mat.columns() == 3);
    REQUIRE(mat.memRows() == 4);
    REQUIRE(mat.memColumns() == 4);

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

TEST_CASE("Cholesky throws MatrixInversionFailure for indefinite matrix", "[matrix_2d]") {
    // Symmetric indefinite matrix — Cholesky must throw MatrixInversionFailure.
    // A = [[2, 1], [1, -1]]  — eigenvalues ~2.41 and ~-1.41
    matrix_2d mat(2, 2);
    mat.put(0, 0, 2.0);
    mat.put(0, 1, 1.0);
    mat.put(1, 0, 1.0);
    mat.put(1, 1, -1.0);

    bool caught = false;
    try {
        mat.cholesky_inverse();
    } catch (const MatrixInversionFailure&) {
        caught = true;
    }
    REQUIRE(caught);
}

TEST_CASE("Cholesky throws MatrixInversionFailure on singular matrix", "[matrix_2d]") {
    // Singular positive semi-definite matrix — dpotrf detects it as not positive definite
    matrix_2d mat(2, 2);
    mat.put(0, 0, 1.0);
    mat.put(0, 1, 2.0);
    mat.put(1, 0, 2.0);
    mat.put(1, 1, 4.0);

    bool caught = false;
    try {
        mat.cholesky_inverse();
    } catch (const MatrixInversionFailure&) {
        caught = true;
    }
    REQUIRE(caught);
}

TEST_CASE("Cholesky inverse with LOWER_IS_CLEARED", "[matrix_2d]") {
    // Upper triangle filled, lower cleared
    matrix_2d mat(3, 3);
    mat.put(0, 0, 4.0);
    mat.put(0, 1, -1.0);
    mat.put(0, 2, -1.0);
    mat.put(1, 1, 3.0);
    mat.put(1, 2, -1.0);
    mat.put(2, 2, 2.0);
    // lower triangle stays zero

    matrix_2d inverse = mat.cholesky_inverse(true);

    REQUIRE(abs(inverse.get(0, 0) - 0.384615) < 0.0001);
    REQUIRE(abs(inverse.get(1, 1) - 0.538462) < 0.0001);
    REQUIRE(abs(inverse.get(2, 2) - 0.846154) < 0.0001);
    // Result should be symmetric
    REQUIRE(abs(inverse.get(0, 1) - inverse.get(1, 0)) < 1e-10);
    REQUIRE(abs(inverse.get(0, 2) - inverse.get(2, 0)) < 1e-10);
    REQUIRE(abs(inverse.get(1, 2) - inverse.get(2, 1)) < 1e-10);
}

// ============================================================================
// Symmetric matrix support
// ============================================================================

TEST_CASE("Symmetric flag default is false", "[matrix_2d][symmetric]") {
    matrix_2d mat(3, 3);
    REQUIRE(mat.is_symmetric() == false);
}

TEST_CASE("Symmetric flag set and get", "[matrix_2d][symmetric]") {
    matrix_2d mat(3, 3);
    mat.set_symmetric(true);
    REQUIRE(mat.is_symmetric() == true);
    mat.set_symmetric(false);
    REQUIRE(mat.is_symmetric() == false);
}

TEST_CASE("Symmetric flag propagates through copy constructor", "[matrix_2d][symmetric]") {
    matrix_2d mat(3, 3);
    mat.put(0, 0, 1.0); mat.put(1, 0, 2.0); mat.put(1, 1, 3.0);
    mat.put(2, 0, 4.0); mat.put(2, 1, 5.0); mat.put(2, 2, 6.0);
    mat.set_symmetric(true);

    matrix_2d copy(mat);
    REQUIRE(copy.is_symmetric() == true);
    REQUIRE(copy.get(0, 1) == 2.0);  // reflects from lower
}

TEST_CASE("Symmetric flag propagates through operator=", "[matrix_2d][symmetric]") {
    matrix_2d mat(3, 3);
    mat.put(0, 0, 1.0); mat.put(1, 0, 2.0); mat.put(1, 1, 3.0);
    mat.set_symmetric(true);

    matrix_2d other(3, 3);
    other = mat;
    REQUIRE(other.is_symmetric() == true);
    REQUIRE(other.get(0, 1) == 2.0);
}

TEST_CASE("Symmetric get reflects upper to lower", "[matrix_2d][symmetric]") {
    // Populate only the lower triangle + diagonal
    matrix_2d mat(4, 4);
    mat.put(0, 0, 1.0);
    mat.put(1, 0, 2.0); mat.put(1, 1, 3.0);
    mat.put(2, 0, 4.0); mat.put(2, 1, 5.0); mat.put(2, 2, 6.0);
    mat.put(3, 0, 7.0); mat.put(3, 1, 8.0); mat.put(3, 2, 9.0); mat.put(3, 3, 10.0);
    mat.set_symmetric(true);

    // Diagonal
    REQUIRE(mat.get(0, 0) == 1.0);
    REQUIRE(mat.get(1, 1) == 3.0);
    REQUIRE(mat.get(2, 2) == 6.0);
    REQUIRE(mat.get(3, 3) == 10.0);

    // Lower triangle (direct)
    REQUIRE(mat.get(1, 0) == 2.0);
    REQUIRE(mat.get(2, 0) == 4.0);
    REQUIRE(mat.get(3, 2) == 9.0);

    // Upper triangle (reflected from lower)
    REQUIRE(mat.get(0, 1) == 2.0);
    REQUIRE(mat.get(0, 2) == 4.0);
    REQUIRE(mat.get(0, 3) == 7.0);
    REQUIRE(mat.get(1, 2) == 5.0);
    REQUIRE(mat.get(1, 3) == 8.0);
    REQUIRE(mat.get(2, 3) == 9.0);
}

TEST_CASE("Symmetric getelementref reflects upper to lower", "[matrix_2d][symmetric]") {
    matrix_2d mat(3, 3);
    mat.put(0, 0, 10.0);
    mat.put(1, 0, 20.0); mat.put(1, 1, 30.0);
    mat.put(2, 0, 40.0); mat.put(2, 1, 50.0); mat.put(2, 2, 60.0);
    mat.set_symmetric(true);

    // const version
    const matrix_2d& cmat = mat;
    REQUIRE(*cmat.getelementref(0, 2) == 40.0);  // reflects (0,2) -> (2,0)
    REQUIRE(*cmat.getelementref(2, 0) == 40.0);  // direct

    // non-const version
    REQUIRE(*mat.getelementref(1, 2) == 50.0);   // reflects (1,2) -> (2,1)
}

TEST_CASE("Symmetric redim clears flag", "[matrix_2d][symmetric]") {
    matrix_2d mat(3, 3);
    mat.set_symmetric(true);
    mat.redim(4, 4);
    REQUIRE(mat.is_symmetric() == false);
}

TEST_CASE("Cholesky inverse with mark_symmetric", "[matrix_2d][symmetric]") {
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

    matrix_2d inverse = mat.cholesky_inverse(false, true);

    REQUIRE(inverse.is_symmetric() == true);
    // Values should match the non-symmetric cholesky inverse
    REQUIRE(abs(inverse.get(0, 0) - 0.384615) < 0.0001);
    REQUIRE(abs(inverse.get(1, 1) - 0.538462) < 0.0001);
    REQUIRE(abs(inverse.get(2, 2) - 0.846154) < 0.0001);
    // Upper triangle via symmetric reflection
    REQUIRE(abs(inverse.get(0, 1) - 0.230769) < 0.0001);
    REQUIRE(abs(inverse.get(0, 2) - 0.307692) < 0.0001);
    REQUIRE(abs(inverse.get(1, 2) - 0.384615) < 0.0001);
    // Symmetry
    REQUIRE(abs(inverse.get(0, 1) - inverse.get(1, 0)) < 1e-10);
    REQUIRE(abs(inverse.get(0, 2) - inverse.get(2, 0)) < 1e-10);
    REQUIRE(abs(inverse.get(1, 2) - inverse.get(2, 1)) < 1e-10);
}

TEST_CASE("Cholesky mark_symmetric matches full inverse", "[matrix_2d][symmetric]") {
    // Larger 5x5 SPD matrix
    matrix_2d full(5, 5);
    matrix_2d sym(5, 5);
    double spd[] = {
        10, 1, 2, 0, 1,
         1, 8, 1, 2, 0,
         2, 1, 7, 1, 1,
         0, 2, 1, 6, 1,
         1, 0, 1, 1, 5
    };
    for (UINT32 i = 0; i < 5; ++i)
        for (UINT32 j = 0; j < 5; ++j) {
            full.put(i, j, spd[i * 5 + j]);
            sym.put(i, j, spd[i * 5 + j]);
        }

    full.cholesky_inverse(false, false);  // full inverse with fillupper
    sym.cholesky_inverse(false, true);    // symmetric inverse

    REQUIRE(sym.is_symmetric() == true);
    REQUIRE(full.is_symmetric() == false);

    // All elements should match
    for (UINT32 i = 0; i < 5; ++i)
        for (UINT32 j = 0; j < 5; ++j)
            REQUIRE(abs(sym.get(i, j) - full.get(i, j)) < 1e-12);
}

TEST_CASE("multiply_sym matches dgemm multiply", "[matrix_2d][symmetric]") {
    // Create a symmetric 4x4 matrix and a 4x2 rhs
    matrix_2d A(4, 4);
    double spd[] = {
        5, 1, 2, 0,
        1, 4, 1, 1,
        2, 1, 6, 1,
        0, 1, 1, 3
    };
    for (UINT32 i = 0; i < 4; ++i)
        for (UINT32 j = 0; j < 4; ++j)
            A.put(i, j, spd[i * 4 + j]);

    matrix_2d B(4, 2);
    B.put(0, 0, 1.0); B.put(0, 1, 5.0);
    B.put(1, 0, 2.0); B.put(1, 1, 6.0);
    B.put(2, 0, 3.0); B.put(2, 1, 7.0);
    B.put(3, 0, 4.0); B.put(3, 1, 8.0);

    // dgemm result
    matrix_2d C_full(4, 2);
    C_full.multiply(A, "N", B, "N");

    // dsymm result — only lower triangle in A_sym
    matrix_2d A_sym(4, 4);
    for (UINT32 i = 0; i < 4; ++i)
        for (UINT32 j = 0; j <= i; ++j)
            A_sym.put(i, j, spd[i * 4 + j]);
    A_sym.set_symmetric(true);

    matrix_2d C_sym(4, 2);
    C_sym.multiply_sym(A_sym, B);

    for (UINT32 i = 0; i < 4; ++i)
        for (UINT32 j = 0; j < 2; ++j)
            REQUIRE(abs(C_sym.get(i, j) - C_full.get(i, j)) < 1e-12);
}

TEST_CASE("multiply_sym with column vector", "[matrix_2d][symmetric]") {
    // This mirrors the Solve() hot path: N_inv * At_Vinv_m
    matrix_2d N(3, 3);
    N.put(0, 0, 4.0);
    N.put(0, 1, -1.0);
    N.put(0, 2, -1.0);
    N.put(1, 0, -1.0);
    N.put(1, 1, 3.0);
    N.put(1, 2, -1.0);
    N.put(2, 0, -1.0);
    N.put(2, 1, -1.0);
    N.put(2, 2, 2.0);

    // Full inverse
    matrix_2d N_full(N);
    N_full.cholesky_inverse(false, false);

    // Symmetric inverse
    matrix_2d N_sym(N);
    N_sym.cholesky_inverse(false, true);

    matrix_2d rhs(3, 1);
    rhs.put(0, 0, 1.0);
    rhs.put(1, 0, 2.0);
    rhs.put(2, 0, 3.0);

    // dgemm
    matrix_2d result_full(3, 1);
    result_full.multiply(N_full, "N", rhs, "N");

    // dsymm
    matrix_2d result_sym(3, 1);
    result_sym.multiply_sym(N_sym, rhs);

    for (UINT32 i = 0; i < 3; ++i)
        REQUIRE(abs(result_sym.get(i, 0) - result_full.get(i, 0)) < 1e-12);
}

TEST_CASE("Symmetric copyelements 3x3 from lower triangle", "[matrix_2d][symmetric]") {
    matrix_2d src(6, 6);
    // Fill lower triangle only
    for (UINT32 i = 0; i < 6; ++i)
        for (UINT32 j = 0; j <= i; ++j)
            src.put(i, j, static_cast<double>(i * 10 + j));
    src.set_symmetric(true);

    matrix_2d dst(3, 3);

    // Copy from lower triangle region (row >= col) — should work directly
    dst.copyelements(0, 0, src, 3, 0, 3, 3);
    REQUIRE(dst.get(0, 0) == 30.0);
    REQUIRE(dst.get(1, 0) == 40.0);
    REQUIRE(dst.get(2, 0) == 50.0);
    REQUIRE(dst.get(0, 1) == 31.0);
    REQUIRE(dst.get(1, 1) == 41.0);
    REQUIRE(dst.get(2, 1) == 51.0);
    REQUIRE(dst.get(0, 2) == 32.0);
    REQUIRE(dst.get(1, 2) == 42.0);
    REQUIRE(dst.get(2, 2) == 52.0);
}

TEST_CASE("Symmetric copyelements 3x3 from upper triangle", "[matrix_2d][symmetric]") {
    matrix_2d src(6, 6);
    // Fill lower triangle only: src(i,j) = i*10+j for j<=i
    for (UINT32 i = 0; i < 6; ++i)
        for (UINT32 j = 0; j <= i; ++j)
            src.put(i, j, static_cast<double>(i * 10 + j));
    src.set_symmetric(true);

    matrix_2d dst(3, 3);

    // Copy from upper triangle region (row_src=0, col_src=3)
    // dst(r,c) = src.get(r, 3+c) which reflects to src(3+c, r) = (3+c)*10 + r
    dst.copyelements(0, 0, src, 0, 3, 3, 3);
    REQUIRE(dst.get(0, 0) == 30.0);  // src.get(0,3) = src(3,0) = 30
    REQUIRE(dst.get(1, 0) == 31.0);  // src.get(1,3) = src(3,1) = 31
    REQUIRE(dst.get(2, 0) == 32.0);  // src.get(2,3) = src(3,2) = 32
    REQUIRE(dst.get(0, 1) == 40.0);  // src.get(0,4) = src(4,0) = 40
    REQUIRE(dst.get(1, 1) == 41.0);  // src.get(1,4) = src(4,1) = 41
    REQUIRE(dst.get(2, 1) == 42.0);  // src.get(2,4) = src(4,2) = 42
    REQUIRE(dst.get(0, 2) == 50.0);  // src.get(0,5) = src(5,0) = 50
    REQUIRE(dst.get(1, 2) == 51.0);  // src.get(1,5) = src(5,1) = 51
    REQUIRE(dst.get(2, 2) == 52.0);  // src.get(2,5) = src(5,2) = 52
}

TEST_CASE("Symmetric blockadd 3x3 from upper triangle", "[matrix_2d][symmetric]") {
    matrix_2d src(6, 6);
    for (UINT32 i = 0; i < 6; ++i)
        for (UINT32 j = 0; j <= i; ++j)
            src.put(i, j, static_cast<double>(i * 10 + j));
    src.set_symmetric(true);

    matrix_2d dst(3, 3);
    // Pre-fill with 1.0 to verify addition
    for (UINT32 i = 0; i < 3; ++i)
        for (UINT32 j = 0; j < 3; ++j)
            dst.put(i, j, 1.0);

    // blockadd from upper triangle (row_src=0, col_src=3)
    dst.blockadd(0, 0, src, 0, 3, 3, 3);

    // dst(r,c) = 1.0 + src.get(r, 3+c)
    REQUIRE(dst.get(0, 0) == 31.0);  // 1 + 30
    REQUIRE(dst.get(1, 0) == 32.0);  // 1 + 31
    REQUIRE(dst.get(2, 0) == 33.0);  // 1 + 32
    REQUIRE(dst.get(0, 1) == 41.0);  // 1 + 40
    REQUIRE(dst.get(1, 1) == 42.0);  // 1 + 41
    REQUIRE(dst.get(2, 1) == 43.0);  // 1 + 42
    REQUIRE(dst.get(0, 2) == 51.0);  // 1 + 50
    REQUIRE(dst.get(1, 2) == 52.0);  // 1 + 51
    REQUIRE(dst.get(2, 2) == 53.0);  // 1 + 52
}

TEST_CASE("Symmetric copyelements generic from upper triangle", "[matrix_2d][symmetric]") {
    matrix_2d src(8, 8);
    for (UINT32 i = 0; i < 8; ++i)
        for (UINT32 j = 0; j <= i; ++j)
            src.put(i, j, static_cast<double>(i * 10 + j));
    src.set_symmetric(true);

    // Copy a 4x4 block from upper triangle (row=0, col=4)
    matrix_2d dst(4, 4);
    dst.copyelements(0, 0, src, 0, 4, 4, 4);

    // Verify via get() which handles reflection
    for (UINT32 i = 0; i < 4; ++i)
        for (UINT32 j = 0; j < 4; ++j)
            REQUIRE(abs(dst.get(i, j) - src.get(i, 4 + j)) < 1e-15);
}

TEST_CASE("End-to-end: cholesky mark_symmetric then multiply_sym", "[matrix_2d][symmetric]") {
    // Simulate the full Solve() path
    matrix_2d N(3, 3);
    N.put(0, 0, 4.0);  N.put(0, 1, -1.0); N.put(0, 2, -1.0);
    N.put(1, 0, -1.0); N.put(1, 1, 3.0);  N.put(1, 2, -1.0);
    N.put(2, 0, -1.0); N.put(2, 1, -1.0); N.put(2, 2, 2.0);

    matrix_2d N_ref(N);

    // Reference: full inverse + dgemm
    N_ref.cholesky_inverse(false, false);

    matrix_2d rhs(3, 1);
    rhs.put(0, 0, 10.0);
    rhs.put(1, 0, 20.0);
    rhs.put(2, 0, 30.0);

    matrix_2d ref_result(3, 1);
    ref_result.multiply(N_ref, "N", rhs, "N");

    // Test: symmetric inverse + dsymm
    N.cholesky_inverse(false, true);
    REQUIRE(N.is_symmetric() == true);

    matrix_2d sym_result(3, 1);
    sym_result.multiply_sym(N, rhs);

    for (UINT32 i = 0; i < 3; ++i)
        REQUIRE(abs(sym_result.get(i, 0) - ref_result.get(i, 0)) < 1e-12);
}

// ============================================================================
// Packed symmetric matrix storage
// ============================================================================

TEST_CASE("Packed allocation has correct buffer size", "[matrix_2d][packed]") {
    matrix_2d mat;
    mat.redim_packed(4);
    REQUIRE(mat.rows() == 4);
    REQUIRE(mat.columns() == 4);
    REQUIRE(mat.is_packed() == true);
    REQUIRE(mat.is_symmetric() == true);
    // n*(n+1)/2 = 10 elements
    REQUIRE(mat.memRows() == 4);
}

TEST_CASE("Packed element access: put and get reflect symmetry", "[matrix_2d][packed]") {
    matrix_2d mat;
    mat.redim_packed(4);

    mat.put(0, 0, 1.0);
    mat.put(1, 0, 2.0); mat.put(1, 1, 3.0);
    mat.put(2, 0, 4.0); mat.put(2, 1, 5.0); mat.put(2, 2, 6.0);
    mat.put(3, 0, 7.0); mat.put(3, 1, 8.0); mat.put(3, 2, 9.0); mat.put(3, 3, 10.0);

    // Diagonal
    REQUIRE(mat.get(0, 0) == 1.0);
    REQUIRE(mat.get(1, 1) == 3.0);
    REQUIRE(mat.get(2, 2) == 6.0);
    REQUIRE(mat.get(3, 3) == 10.0);

    // Lower triangle (direct)
    REQUIRE(mat.get(1, 0) == 2.0);
    REQUIRE(mat.get(2, 0) == 4.0);
    REQUIRE(mat.get(3, 2) == 9.0);

    // Upper triangle (reflected from lower)
    REQUIRE(mat.get(0, 1) == 2.0);
    REQUIRE(mat.get(0, 2) == 4.0);
    REQUIRE(mat.get(0, 3) == 7.0);
    REQUIRE(mat.get(1, 2) == 5.0);
    REQUIRE(mat.get(1, 3) == 8.0);
    REQUIRE(mat.get(2, 3) == 9.0);
}

TEST_CASE("Packed put via upper triangle maps to lower", "[matrix_2d][packed]") {
    matrix_2d mat;
    mat.redim_packed(3);

    mat.put(0, 1, 42.0);
    REQUIRE(mat.get(1, 0) == 42.0);
    REQUIRE(mat.get(0, 1) == 42.0);
}

TEST_CASE("Packed elementadd skips upper triangle", "[matrix_2d][packed]") {
    matrix_2d mat;
    mat.redim_packed(3);

    mat.elementadd(1, 0, 10.0);
    mat.elementadd(0, 1, 5.0);  // skipped — upper triangle
    REQUIRE(mat.get(1, 0) == 10.0);
    REQUIRE(mat.get(0, 1) == 10.0);

    mat.elementadd(1, 0, 3.0);  // lower triangle — accumulates
    REQUIRE(mat.get(1, 0) == 13.0);
}

TEST_CASE("Packed cholesky_inverse matches full inverse", "[matrix_2d][packed]") {
    matrix_2d full(3, 3);
    full.put(0, 0, 4.0);  full.put(0, 1, -1.0); full.put(0, 2, -1.0);
    full.put(1, 0, -1.0); full.put(1, 1, 3.0);  full.put(1, 2, -1.0);
    full.put(2, 0, -1.0); full.put(2, 1, -1.0); full.put(2, 2, 2.0);

    matrix_2d packed;
    packed.redim_packed(3);
    packed.put(0, 0, 4.0);
    packed.put(1, 0, -1.0); packed.put(1, 1, 3.0);
    packed.put(2, 0, -1.0); packed.put(2, 1, -1.0); packed.put(2, 2, 2.0);

    full.cholesky_inverse();
    packed.cholesky_inverse();

    REQUIRE(packed.is_packed() == true);
    for (UINT32 i = 0; i < 3; ++i)
        for (UINT32 j = 0; j < 3; ++j)
            REQUIRE(abs(packed.get(i, j) - full.get(i, j)) < 1e-10);
}

TEST_CASE("Packed cholesky_inverse throws on indefinite matrix", "[matrix_2d][packed]") {
    matrix_2d mat;
    mat.redim_packed(2);
    mat.put(0, 0, 2.0);
    mat.put(1, 0, 1.0);
    mat.put(1, 1, -1.0);

    bool caught = false;
    try { mat.cholesky_inverse(); }
    catch (const MatrixInversionFailure&) { caught = true; }
    REQUIRE(caught);
}

TEST_CASE("Packed multiply_sym with vector matches full result", "[matrix_2d][packed]") {
    matrix_2d full(3, 3);
    full.put(0, 0, 4.0);  full.put(0, 1, -1.0); full.put(0, 2, -1.0);
    full.put(1, 0, -1.0); full.put(1, 1, 3.0);  full.put(1, 2, -1.0);
    full.put(2, 0, -1.0); full.put(2, 1, -1.0); full.put(2, 2, 2.0);

    matrix_2d packed;
    packed.redim_packed(3);
    packed.put(0, 0, 4.0);
    packed.put(1, 0, -1.0); packed.put(1, 1, 3.0);
    packed.put(2, 0, -1.0); packed.put(2, 1, -1.0); packed.put(2, 2, 2.0);

    matrix_2d rhs(3, 1);
    rhs.put(0, 0, 1.0); rhs.put(1, 0, 2.0); rhs.put(2, 0, 3.0);

    matrix_2d result_full(3, 1);
    result_full.multiply(full, "N", rhs, "N");

    matrix_2d result_packed(3, 1);
    result_packed.multiply_sym(packed, rhs);

    for (UINT32 i = 0; i < 3; ++i)
        REQUIRE(abs(result_packed.get(i, 0) - result_full.get(i, 0)) < 1e-12);
}

TEST_CASE("Packed multiply_sym with multi-column rhs", "[matrix_2d][packed]") {
    matrix_2d packed;
    packed.redim_packed(4);
    double spd[] = {5, 1, 2, 0, 1, 4, 1, 1, 2, 1, 6, 1, 0, 1, 1, 3};
    for (UINT32 i = 0; i < 4; ++i)
        for (UINT32 j = 0; j <= i; ++j)
            packed.put(i, j, spd[i * 4 + j]);

    matrix_2d full(4, 4);
    for (UINT32 i = 0; i < 4; ++i)
        for (UINT32 j = 0; j < 4; ++j)
            full.put(i, j, spd[i * 4 + j]);

    matrix_2d B(4, 2);
    B.put(0, 0, 1.0); B.put(0, 1, 5.0);
    B.put(1, 0, 2.0); B.put(1, 1, 6.0);
    B.put(2, 0, 3.0); B.put(2, 1, 7.0);
    B.put(3, 0, 4.0); B.put(3, 1, 8.0);

    matrix_2d C_full(4, 2);
    C_full.multiply(full, "N", B, "N");

    matrix_2d C_packed(4, 2);
    C_packed.multiply_sym(packed, B);

    for (UINT32 i = 0; i < 4; ++i)
        for (UINT32 j = 0; j < 2; ++j)
            REQUIRE(abs(C_packed.get(i, j) - C_full.get(i, j)) < 1e-12);
}

TEST_CASE("Packed scale_symmetric_diagonal", "[matrix_2d][packed]") {
    matrix_2d packed;
    packed.redim_packed(3);
    packed.put(0, 0, 4.0);
    packed.put(1, 0, -1.0); packed.put(1, 1, 3.0);
    packed.put(2, 0, -1.0); packed.put(2, 1, -1.0); packed.put(2, 2, 2.0);

    matrix_2d full(3, 3);
    full.put(0, 0, 4.0);  full.put(0, 1, -1.0); full.put(0, 2, -1.0);
    full.put(1, 0, -1.0); full.put(1, 1, 3.0);  full.put(1, 2, -1.0);
    full.put(2, 0, -1.0); full.put(2, 1, -1.0); full.put(2, 2, 2.0);
    full.set_symmetric(true);

    double diag[] = {2.0, 3.0, 0.5};
    packed.scale_symmetric_diagonal(diag);
    full.scale_symmetric_diagonal(diag);

    for (UINT32 i = 0; i < 3; ++i)
        for (UINT32 j = 0; j < 3; ++j)
            REQUIRE(abs(packed.get(i, j) - full.get(i, j)) < 1e-12);
}

TEST_CASE("Packed blockadd diagonal block", "[matrix_2d][packed]") {
    matrix_2d packed;
    packed.redim_packed(4);

    // Diagonal blocks only write lower triangle via elementadd
    // (mirrors how AddMsrtoNormalsVar works)
    packed.elementadd(1, 1, 10.0);
    packed.elementadd(2, 1, 20.0);
    packed.elementadd(2, 2, 40.0);

    REQUIRE(packed.get(1, 1) == 10.0);
    REQUIRE(packed.get(2, 1) == 20.0);
    REQUIRE(packed.get(1, 2) == 20.0);
    REQUIRE(packed.get(2, 2) == 40.0);
    REQUIRE(packed.get(0, 0) == 0.0);
}

TEST_CASE("Packed submatrix extraction", "[matrix_2d][packed]") {
    matrix_2d packed;
    packed.redim_packed(4);
    packed.put(0, 0, 1.0);
    packed.put(1, 0, 2.0); packed.put(1, 1, 3.0);
    packed.put(2, 0, 4.0); packed.put(2, 1, 5.0); packed.put(2, 2, 6.0);
    packed.put(3, 0, 7.0); packed.put(3, 1, 8.0); packed.put(3, 2, 9.0); packed.put(3, 3, 10.0);

    matrix_2d sub = packed.submatrix(1, 1, 2, 2);
    REQUIRE(sub.rows() == 2);
    REQUIRE(sub.columns() == 2);
    REQUIRE(sub.get(0, 0) == 3.0);
    REQUIRE(sub.get(0, 1) == 5.0);
    REQUIRE(sub.get(1, 0) == 5.0);
    REQUIRE(sub.get(1, 1) == 6.0);
}

TEST_CASE("Packed copy constructor", "[matrix_2d][packed]") {
    matrix_2d orig;
    orig.redim_packed(3);
    orig.put(0, 0, 1.0);
    orig.put(1, 0, 2.0); orig.put(1, 1, 3.0);
    orig.put(2, 0, 4.0); orig.put(2, 1, 5.0); orig.put(2, 2, 6.0);

    matrix_2d copy(orig);
    REQUIRE(copy.is_packed() == true);
    REQUIRE(copy.is_symmetric() == true);
    for (UINT32 i = 0; i < 3; ++i)
        for (UINT32 j = 0; j < 3; ++j)
            REQUIRE(copy.get(i, j) == orig.get(i, j));
}

TEST_CASE("Packed operator= packed-to-packed", "[matrix_2d][packed]") {
    matrix_2d a;
    a.redim_packed(3);
    a.put(0, 0, 10.0); a.put(1, 0, 20.0); a.put(1, 1, 30.0);
    a.put(2, 0, 40.0); a.put(2, 1, 50.0); a.put(2, 2, 60.0);

    matrix_2d b;
    b.redim_packed(2);
    b = a;
    REQUIRE(b.is_packed() == true);
    REQUIRE(b.rows() == 3);
    for (UINT32 i = 0; i < 3; ++i)
        for (UINT32 j = 0; j < 3; ++j)
            REQUIRE(b.get(i, j) == a.get(i, j));
}

TEST_CASE("Packed operator= full-to-packed preserves full", "[matrix_2d][packed]") {
    matrix_2d full(3, 3);
    full.put(0, 0, 1.0); full.put(0, 1, 2.0); full.put(0, 2, 3.0);
    full.put(1, 0, 4.0); full.put(1, 1, 5.0); full.put(1, 2, 6.0);
    full.put(2, 0, 7.0); full.put(2, 1, 8.0); full.put(2, 2, 9.0);

    matrix_2d packed;
    packed.redim_packed(3);
    packed = full;
    // Assigning a non-packed to a packed should convert to non-packed
    REQUIRE(packed.is_packed() == false);
    for (UINT32 i = 0; i < 3; ++i)
        for (UINT32 j = 0; j < 3; ++j)
            REQUIRE(packed.get(i, j) == full.get(i, j));
}

TEST_CASE("Packed redim_packed reuses buffer when shrinking", "[matrix_2d][packed]") {
    matrix_2d mat;
    mat.redim_packed(5);
    mat.put(0, 0, 99.0);

    mat.redim_packed(3);
    REQUIRE(mat.rows() == 3);
    REQUIRE(mat.is_packed() == true);
    REQUIRE(mat.get(0, 0) == 0.0);
}

TEST_CASE("Packed zero", "[matrix_2d][packed]") {
    matrix_2d mat;
    mat.redim_packed(3);
    mat.put(0, 0, 1.0); mat.put(1, 0, 2.0); mat.put(2, 2, 3.0);
    mat.zero();
    for (UINT32 i = 0; i < 3; ++i)
        for (UINT32 j = 0; j < 3; ++j)
            REQUIRE(mat.get(i, j) == 0.0);
}

TEST_CASE("Packed scale", "[matrix_2d][packed]") {
    matrix_2d mat;
    mat.redim_packed(3);
    mat.put(0, 0, 2.0);
    mat.put(1, 0, 3.0); mat.put(1, 1, 4.0);
    mat.put(2, 0, 5.0); mat.put(2, 1, 6.0); mat.put(2, 2, 7.0);
    mat.scale(2.0);

    REQUIRE(mat.get(0, 0) == 4.0);
    REQUIRE(mat.get(1, 0) == 6.0);
    REQUIRE(mat.get(2, 2) == 14.0);
    REQUIRE(mat.get(0, 1) == 6.0); // reflected
}

TEST_CASE("Packed end-to-end: cholesky then multiply_sym", "[matrix_2d][packed]") {
    matrix_2d full(3, 3);
    full.put(0, 0, 4.0);  full.put(0, 1, -1.0); full.put(0, 2, -1.0);
    full.put(1, 0, -1.0); full.put(1, 1, 3.0);  full.put(1, 2, -1.0);
    full.put(2, 0, -1.0); full.put(2, 1, -1.0); full.put(2, 2, 2.0);

    matrix_2d packed;
    packed.redim_packed(3);
    packed.put(0, 0, 4.0);
    packed.put(1, 0, -1.0); packed.put(1, 1, 3.0);
    packed.put(2, 0, -1.0); packed.put(2, 1, -1.0); packed.put(2, 2, 2.0);

    // Reference
    full.cholesky_inverse(false, false);
    matrix_2d rhs(3, 1);
    rhs.put(0, 0, 10.0); rhs.put(1, 0, 20.0); rhs.put(2, 0, 30.0);
    matrix_2d ref_result(3, 1);
    ref_result.multiply(full, "N", rhs, "N");

    // Packed path
    packed.cholesky_inverse();
    matrix_2d packed_result(3, 1);
    packed_result.multiply_sym(packed, rhs);

    for (UINT32 i = 0; i < 3; ++i)
        REQUIRE(abs(packed_result.get(i, 0) - ref_result.get(i, 0)) < 1e-12);
}

TEST_CASE("Packed copyelements from packed source", "[matrix_2d][packed]") {
    matrix_2d packed;
    packed.redim_packed(6);
    for (UINT32 i = 0; i < 6; ++i)
        for (UINT32 j = 0; j <= i; ++j)
            packed.put(i, j, static_cast<double>(i * 10 + j));

    matrix_2d dst(3, 3);
    dst.copyelements(0, 0, packed, 3, 0, 3, 3);

    REQUIRE(dst.get(0, 0) == 30.0);
    REQUIRE(dst.get(1, 0) == 40.0);
    REQUIRE(dst.get(2, 0) == 50.0);
    REQUIRE(dst.get(0, 1) == 31.0);
    REQUIRE(dst.get(1, 1) == 41.0);
    REQUIRE(dst.get(2, 1) == 51.0);
}

TEST_CASE("Packed 5x5 cholesky_inverse matches full", "[matrix_2d][packed]") {
    double spd[] = {
        10, 1, 2, 0, 1,
         1, 8, 1, 2, 0,
         2, 1, 7, 1, 1,
         0, 2, 1, 6, 1,
         1, 0, 1, 1, 5
    };

    matrix_2d full(5, 5);
    matrix_2d packed;
    packed.redim_packed(5);

    for (UINT32 i = 0; i < 5; ++i)
        for (UINT32 j = 0; j < 5; ++j) {
            full.put(i, j, spd[i * 5 + j]);
            if (j <= i)
                packed.put(i, j, spd[i * 5 + j]);
        }

    full.cholesky_inverse(false, false);
    packed.cholesky_inverse();

    for (UINT32 i = 0; i < 5; ++i)
        for (UINT32 j = 0; j < 5; ++j)
            REQUIRE(abs(packed.get(i, j) - full.get(i, j)) < 1e-12);
}

// ================================================================
// In-place mmap buffer tests (AttachMappedFileRegion / DetachMappedFileRegion)
// ================================================================

TEST_CASE("In-place mmap round-trip: full matrix", "[matrix_2d][mmap]") {
    // Create a 4x3 full matrix
    matrix_2d orig(4, 3);
    for (UINT32 r = 0; r < 4; ++r)
        for (UINT32 c = 0; c < 3; ++c)
            orig.put(r, c, (r + 1) * 10.0 + c);
    orig.compute_maximum_value();

    // Allocate a buffer to simulate an mmap region
    std::size_t region_size = orig.get_size();
    std::vector<char> region(region_size, 0);
    void* addr = region.data();

    // Write the matrix to the "mmap" region
    orig.WriteMappedFileRegion(addr);

    // Attach in-place — should point _buffer at the data in the region
    matrix_2d attached;
    attached.AttachMappedFileRegion(addr);

    REQUIRE(!attached.owns_buffer());
    REQUIRE(attached.rows() == 4);
    REQUIRE(attached.columns() == 3);
    REQUIRE(attached.matrixType() == mtx_full);

    // Verify all elements match
    for (UINT32 r = 0; r < 4; ++r)
        for (UINT32 c = 0; c < 3; ++c)
            REQUIRE(attached.get(r, c) == orig.get(r, c));

    // Modify in-place (writes directly to region)
    attached.put(2, 1, 999.0);
    REQUIRE(attached.get(2, 1) == 999.0);

    // Detach
    attached.DetachMappedFileRegion(addr);
    REQUIRE(attached.empty());
    REQUIRE(attached.owns_buffer());

    // Read back with copy-based Read to verify data persisted
    matrix_2d readback;
    readback.ReadMappedFileRegion(addr);
    REQUIRE(readback.get(2, 1) == 999.0);
    REQUIRE(readback.get(0, 0) == 10.0);
}

TEST_CASE("In-place mmap round-trip: packed lower-triangular", "[matrix_2d][mmap][packed]") {
    // Create a 4x4 packed lower-triangular matrix
    matrix_2d orig;
    orig.redim_packed(4);
    for (UINT32 i = 0; i < 4; ++i)
        for (UINT32 j = 0; j <= i; ++j)
            orig.put(i, j, (i + 1) * 10.0 + j);
    orig.compute_maximum_value();

    std::size_t region_size = orig.get_size();
    std::vector<char> region(region_size, 0);
    void* addr = region.data();

    orig.WriteMappedFileRegion(addr);

    matrix_2d attached;
    attached.AttachMappedFileRegion(addr);

    REQUIRE(!attached.owns_buffer());
    REQUIRE(attached.is_packed());
    REQUIRE(attached.is_symmetric());
    REQUIRE(attached.rows() == 4);

    // Verify all elements
    for (UINT32 i = 0; i < 4; ++i)
        for (UINT32 j = 0; j <= i; ++j)
            REQUIRE(attached.get(i, j) == orig.get(i, j));

    // Symmetry: upper = lower
    REQUIRE(attached.get(0, 1) == attached.get(1, 0));

    // Modify in-place
    attached.put(3, 0, -7.5);

    // WriteMappedFileRegion on in-place buffer should only write footer
    attached.WriteMappedFileRegion(addr);

    // deallocate detaches without freeing
    attached.deallocate();
    REQUIRE(attached.empty());
    REQUIRE(attached.owns_buffer());

    // Read back to verify
    matrix_2d readback;
    readback.ReadMappedFileRegion(addr);
    REQUIRE(readback.get(3, 0) == -7.5);
    REQUIRE(readback.get(1, 0) == 20.0);
}

TEST_CASE("In-place mmap: deallocate does not crash", "[matrix_2d][mmap]") {
    matrix_2d orig(3, 3);
    for (UINT32 i = 0; i < 9; ++i)
        orig.put(i / 3, i % 3, i + 1.0);

    std::size_t region_size = orig.get_size();
    std::vector<char> region(region_size, 0);
    void* addr = region.data();
    orig.WriteMappedFileRegion(addr);

    matrix_2d attached;
    attached.AttachMappedFileRegion(addr);
    REQUIRE(!attached.owns_buffer());

    // deallocate should not crash (no delete[] on mmap pointer)
    attached.deallocate();
    REQUIRE(attached.empty());
    REQUIRE(attached.owns_buffer());

    // Re-attach and verify data still intact
    attached.AttachMappedFileRegion(addr);
    REQUIRE(attached.get(1, 1) == 5.0);
    attached.deallocate();
}

TEST_CASE("In-place mmap: copy constructor makes owned copy", "[matrix_2d][mmap]") {
    matrix_2d orig(2, 2);
    orig.put(0, 0, 1.0); orig.put(0, 1, 2.0);
    orig.put(1, 0, 3.0); orig.put(1, 1, 4.0);

    std::size_t region_size = orig.get_size();
    std::vector<char> region(region_size, 0);
    void* addr = region.data();
    orig.WriteMappedFileRegion(addr);

    matrix_2d attached;
    attached.AttachMappedFileRegion(addr);

    // Copy constructor should create an owned copy
    matrix_2d copy(attached);
    REQUIRE(copy.owns_buffer());
    REQUIRE(copy.get(1, 1) == 4.0);

    // Modifying copy should not affect mmap
    copy.put(1, 1, 99.0);
    REQUIRE(attached.get(1, 1) == 4.0);

    attached.deallocate();
}

TEST_CASE("In-place mmap: operator= into owned buffer", "[matrix_2d][mmap]") {
    matrix_2d orig(3, 2);
    for (UINT32 r = 0; r < 3; ++r)
        for (UINT32 c = 0; c < 2; ++c)
            orig.put(r, c, r * 2.0 + c);

    std::size_t region_size = orig.get_size();
    std::vector<char> region(region_size, 0);
    void* addr = region.data();
    orig.WriteMappedFileRegion(addr);

    matrix_2d attached;
    attached.AttachMappedFileRegion(addr);

    // Assign to a pre-allocated owned matrix of same size
    matrix_2d dest(3, 2);
    dest = attached;
    REQUIRE(dest.owns_buffer());
    REQUIRE(dest.get(2, 1) == 5.0);

    attached.deallocate();
}

TEST_CASE("In-place mmap: region size alignment", "[matrix_2d][mmap]") {
    // Verify that get_size() returns 8-byte aligned sizes for proper mmap alignment
    matrix_2d full(7, 5);
    REQUIRE(full.get_size() % 8 == 0);

    matrix_2d packed;
    packed.redim_packed(11);
    REQUIRE(packed.get_size() % 8 == 0);

    matrix_2d col_vec(100, 1);
    REQUIRE(col_vec.get_size() % 8 == 0);
}

TEST_CASE("operator<< binary size matches get_size for mmap regions", "[matrix_2d][mmap]") {
    // Verify that operator<< writes exactly get_size() bytes in binary mode.
    // A mismatch here causes mmap region offsets to be wrong for all blocks
    // after the first, leading to corrupted header reads.

    auto check_stream_size = [](matrix_2d& m, const std::string& label) {
        std::ostringstream oss;
        oss.iword(0) = binary;
        oss << m;
        std::size_t stream_bytes = oss.str().size();
        std::size_t region_bytes = m.get_size();
        REQUIRE(stream_bytes == region_bytes);
    };

    // Full matrix (column vector — the case that crashed)
    matrix_2d col_vec(50, 1);
    for (UINT32 r = 0; r < 50; ++r)
        col_vec.put(r, 0, static_cast<double>(r));
    check_stream_size(col_vec, "column vector");

    // Full matrix
    matrix_2d full(7, 5);
    for (UINT32 r = 0; r < 7; ++r)
        for (UINT32 c = 0; c < 5; ++c)
            full.put(r, c, static_cast<double>(r * 10 + c));
    check_stream_size(full, "full 7x5");

    // Packed lower-triangular
    matrix_2d packed;
    packed.redim_packed(6);
    for (UINT32 r = 0; r < 6; ++r)
        for (UINT32 c = 0; c <= r; ++c)
            packed.put(r, c, static_cast<double>(r * 10 + c));
    check_stream_size(packed, "packed 6x6");

    // Full matrix with slack (redim + shrink, simulating staged adjustment)
    matrix_2d with_slack;
    with_slack.redim(100, 1);
    with_slack.shrink(30, 0);  // _rows=70, _mem_rows=100
    for (UINT32 r = 0; r < 70; ++r)
        with_slack.put(r, 0, static_cast<double>(r));
    check_stream_size(with_slack, "column vector with slack");
}
