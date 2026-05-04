//============================================================================
// Name         : dnamatrix_contiguous.hpp
// Author       : Roger Fraser
// Contributors : Dale Roberts <dale.o.roberts@gmail.com>
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
// Description  : DynAdjust Matrix library
//============================================================================

#ifndef DNAMATRIX_CONTIGUOUS_H_
#define DNAMATRIX_CONTIGUOUS_H_

/// \cond
#include <cassert>
#include <cstring>
/// \endcond

#include <include/config/dnatypes.hpp>
#include <include/config/dnaversion.hpp>
#include <include/config/dnaexports.hpp>
#include <include/exception/dnaexception.hpp>
#include <include/functions/dnatemplatecalcfuncs.hpp>
#include <include/memory/dnamemory_handler.hpp>

#ifdef _MSDEBUG
#include <include/ide/trace.hpp>
#endif

#if defined(__APPLE__)
// Apple Accelerate Framework

#ifndef ACCELERATE_NEW_LAPACK
#define ACCELERATE_NEW_LAPACK
#endif

#include <Accelerate/Accelerate.h>

#ifdef USE_ILP64

#ifndef ACCELERATE_LAPACK_ILP64
#define ACCELERATE_LAPACK_ILP64
#endif

#define LAPACK_SYMBOL_PREFIX
#define LAPACK_FORTRAN_SUFFIX
#define LAPACK_SYMBOL_SUFFIX $NEWLAPACK$ILP64
#define BLAS_SYMBOL_PREFIX cblas_
#define BLAS_FORTRAN_SUFFIX
#define BLAS_SYMBOL_SUFFIX $NEWLAPACK$ILP64
typedef long lapack_int;

#else

#define LAPACK_SYMBOL_PREFIX
#define LAPACK_FORTRAN_SUFFIX
#define LAPACK_SYMBOL_SUFFIX $NEWLAPACK
#define BLAS_SYMBOL_PREFIX cblas_
#define BLAS_FORTRAN_SUFFIX
#define BLAS_SYMBOL_SUFFIX $NEWLAPACK
typedef int lapack_int;

#endif

// End - Apple Accelerate Framework

#elif defined(USE_MKL) || defined(__MKL__)
// Intel MKL
// #pragma message("Using Intel MKL for LAPACK/BLAS")

#include <mkl.h>

#ifdef USE_ILP64

// Force Intel MKL to use ILP64 interface
#ifndef MKL_ILP64
#define MKL_ILP64
#endif

#define LAPACK_SYMBOL_PREFIX
#define LAPACK_FORTRAN_SUFFIX _
#define LAPACK_SYMBOL_SUFFIX
#define BLAS_SYMBOL_PREFIX cblas_
#define BLAS_FORTRAN_SUFFIX
#define BLAS_SYMBOL_SUFFIX

#define lapack_int MKL_INT

#else

#define LAPACK_SYMBOL_PREFIX
#define LAPACK_FORTRAN_SUFFIX _
#define LAPACK_SYMBOL_SUFFIX
#define BLAS_SYMBOL_PREFIX cblas_
#define BLAS_FORTRAN_SUFFIX
#define BLAS_SYMBOL_SUFFIX

#define lapack_int MKL_INT

#endif

// End - Intel MKL

#else
// Default LAPACK/BLAS
// #pragma message("Using default LAPACK/BLAS")

#include <cblas.h>

#ifdef USE_ILP64

#define LAPACK_SYMBOL_PREFIX
#define LAPACK_FORTRAN_SUFFIX _
#define LAPACK_SYMBOL_SUFFIX 64_
#define BLAS_SYMBOL_PREFIX cblas_
#define BLAS_FORTRAN_SUFFIX
#define BLAS_SYMBOL_SUFFIX 64_
typedef long lapack_int;

#else

#define LAPACK_SYMBOL_PREFIX
#define LAPACK_FORTRAN_SUFFIX _
#define LAPACK_SYMBOL_SUFFIX
#define BLAS_SYMBOL_PREFIX cblas_
#define BLAS_FORTRAN_SUFFIX
#define BLAS_SYMBOL_SUFFIX
typedef int lapack_int;

#endif

// End - Default LAPACK/BLAS
#endif

#ifdef USE_ILP64
// Ensure that lapack_int is 64 bits for ILP64
static_assert(sizeof(lapack_int) == 8, "ILP64 interface requires 64-bit integers");
#else
// Ensure that lapack_int is 32 bits for LP64
static_assert(sizeof(lapack_int) == 4, "LP64 interface requires 32-bit integers");
#endif

#define DNAMATRIX_INDEX(no_rows, no_cols, row, column) column* no_rows + row
#define DNAMATRIX_ELEMENT(A, no_rows, no_cols, row, column) A[DNAMATRIX_INDEX(no_rows, no_cols, row, column)]

#define LAPACK_FUNC_CONCAT(name, prefix, suffix, suffix2) prefix##name##suffix##suffix2
#define LAPACK_FUNC_EXPAND(name, prefix, suffix, suffix2) LAPACK_FUNC_CONCAT(name, prefix, suffix, suffix2)
#define LAPACK_FUNC(name) LAPACK_FUNC_EXPAND(name, LAPACK_SYMBOL_PREFIX, LAPACK_FORTRAN_SUFFIX, LAPACK_SYMBOL_SUFFIX)

#define BLAS_FUNC_CONCAT(name, prefix, suffix, suffix2) prefix##name##suffix##suffix2
#define BLAS_FUNC_EXPAND(name, prefix, suffix, suffix2) BLAS_FUNC_CONCAT(name, prefix, suffix, suffix2)
#define BLAS_FUNC(name) BLAS_FUNC_EXPAND(name, BLAS_SYMBOL_PREFIX, BLAS_FORTRAN_SUFFIX, BLAS_SYMBOL_SUFFIX)

#ifndef USE_MKL
extern "C" {
void LAPACK_FUNC(dpotrf)(const char* uplo, const lapack_int* n, double* a, const lapack_int* lda, lapack_int* info);
void LAPACK_FUNC(dpotri)(const char* uplo, const lapack_int* n, double* a, const lapack_int* lda, lapack_int* info);
void LAPACK_FUNC(dpptrf)(const char* uplo, const lapack_int* n, double* ap, lapack_int* info);
void LAPACK_FUNC(dpptri)(const char* uplo, const lapack_int* n, double* ap, lapack_int* info);
void LAPACK_FUNC(dsytrf)(const char* uplo, const lapack_int* n, double* a, const lapack_int* lda,
                         lapack_int* ipiv, double* work, const lapack_int* lwork, lapack_int* info);
void LAPACK_FUNC(dsytri)(const char* uplo, const lapack_int* n, double* a, const lapack_int* lda,
                         const lapack_int* ipiv, double* work, lapack_int* info);
void LAPACK_FUNC(dpotrs)(const char* uplo, const lapack_int* n, const lapack_int* nrhs,
                         const double* a, const lapack_int* lda, double* b, const lapack_int* ldb, lapack_int* info);
void BLAS_FUNC(dgemm)(const enum CBLAS_ORDER ORDER, const enum CBLAS_TRANSPOSE TRANSA,
                      const enum CBLAS_TRANSPOSE TRANSB, const lapack_int M, const lapack_int N, const lapack_int K,
                      const double ALPHA, const double* A, const lapack_int LDA, const double* B, const lapack_int LDB,
                      const double BETA, double* C, const lapack_int LDC);
void BLAS_FUNC(dsymm)(const enum CBLAS_ORDER ORDER, const enum CBLAS_SIDE SIDE, const enum CBLAS_UPLO UPLO,
                      const lapack_int M, const lapack_int N,
                      const double ALPHA, const double* A, const lapack_int LDA, const double* B, const lapack_int LDB,
                      const double BETA, double* C, const lapack_int LDC);
void BLAS_FUNC(dspmv)(const enum CBLAS_ORDER ORDER, const enum CBLAS_UPLO UPLO,
                      const lapack_int N, const double ALPHA, const double* AP,
                      const double* X, const lapack_int INCX,
                      const double BETA, double* Y, const lapack_int INCY);
}
#endif

using namespace dynadjust::memory;
using namespace dynadjust::exception;

namespace dynadjust {
namespace math {

class MatrixInversionFailure : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Set/get the maximum BLAS thread count used by matrix operations.
void set_max_blas_threads(int n);
int get_max_blas_threads();

class matrix_2d;
typedef std::vector<matrix_2d> v_mat_2d, *pv_mat_2d;
typedef v_mat_2d::iterator _it_v_mat_2d;
typedef v_mat_2d::const_iterator _it_v_mat_2d_const;
typedef std::vector<v_mat_2d> vv_mat_2d;

template <typename T> std::size_t byteSize(const UINT32 elements = 1) { return elements * sizeof(T); }

class matrix_2d : public new_handler_support<matrix_2d> {
  public:
    // Constructors/deconstructors
    matrix_2d();
    matrix_2d(const UINT32& rows, const UINT32& columns); // explicit constructor
    matrix_2d(const UINT32& rows, const UINT32& columns, const double data[], const std::size_t& data_size,
              const UINT32& matrix_type = mtx_full);
    matrix_2d(const matrix_2d&); // copy constructor
    ~matrix_2d();                // destructor

    inline bool empty() { return _buffer == nullptr; }

    std::size_t get_size();

    ///////////////////////////////////////////////////////////////////////
    // Get
    inline UINT32 memRows() const { return _mem_rows; }
    inline UINT32 memColumns() const { return _mem_cols; }
    inline UINT32 rows() const { return _rows; }
    inline UINT32 columns() const { return _cols; }
    inline double* getbuffer() const { return _buffer; }

    inline double dense_get(const UINT32& row, const UINT32& column) const {
        assert(_buffer != nullptr);
        assert(!_packed && "dense_get(): matrix must not use packed storage");
        assert(!_symmetric && "dense_get(): matrix must not use symmetric element redirection");
        assert(row < _rows && "dense_get(): row out of bounds");
        assert(column < _cols && "dense_get(): column out of bounds");
        return _buffer[static_cast<std::size_t>(column) * _mem_rows + row];
    }

    inline void dense_put(const UINT32& row, const UINT32& column, const double& value) {
        assert(_buffer != nullptr);
        assert(!_packed && "dense_put(): matrix must not use packed storage");
        assert(!_symmetric && "dense_put(): matrix must not use symmetric element redirection");
        assert(row < _rows && "dense_put(): row out of bounds");
        assert(column < _cols && "dense_put(): column out of bounds");
        _buffer[static_cast<std::size_t>(column) * _mem_rows + row] = value;
    }

    inline void dense_add(const UINT32& row, const UINT32& column, const double& increment) {
        assert(_buffer != nullptr);
        assert(!_packed && "dense_add(): matrix must not use packed storage");
        assert(!_symmetric && "dense_add(): matrix must not use symmetric element redirection");
        assert(row < _rows && "dense_add(): row out of bounds");
        assert(column < _cols && "dense_add(): column out of bounds");
        _buffer[static_cast<std::size_t>(column) * _mem_rows + row] += increment;
    }

    inline double* dense_ptr(const UINT32& row, const UINT32& column) const {
        assert(_buffer != nullptr);
        assert(!_packed && "dense_ptr(): matrix must not use packed storage");
        assert(!_symmetric && "dense_ptr(): matrix must not use symmetric element redirection");
        assert(row < _rows && "dense_ptr(): row out of bounds");
        assert(column < _cols && "dense_ptr(): column out of bounds");
        return _buffer + static_cast<std::size_t>(column) * _mem_rows + row;
    }

    // element retrieval
    // see DNAMATRIX_ROW_WISE
    inline double& get(const UINT32& row, const UINT32& column) const {
        assert(_buffer != nullptr);
        assert(row < _mem_rows && "get(): row out of bounds");
        assert(column < _mem_cols && "get(): column out of bounds");
        if (_packed) {
            UINT32 i = row, j = column;
            if (i < j) { i = column; j = row; }
            return _buffer[packed_index(_rows, i, j)];
        }
        if (_symmetric && row < column)
            return DNAMATRIX_ELEMENT(_buffer, _mem_rows, _mem_cols, column, row);
        return DNAMATRIX_ELEMENT(_buffer, _mem_rows, _mem_cols, row, column);
    }
    inline double* getbuffer(const UINT32& row, const UINT32& column) const {
        assert(_buffer != nullptr);
        assert(!_packed && "getbuffer(row,col): not valid for packed storage");
        assert(row < _mem_rows && "getbuffer(): row out of bounds");
        assert(column < _mem_cols && "getbuffer(): column out of bounds");
        return _buffer + DNAMATRIX_INDEX(_mem_rows, _mem_cols, row, column);
    }

    void submatrix(const UINT32& row_begin, const UINT32& col_begin, matrix_2d* dest, const UINT32& rows,
                   const UINT32& columns) const;
    matrix_2d
    submatrix(const UINT32& row_begin, const UINT32& col_begin, const UINT32& rows, const UINT32& columns) const;

    inline double maxvalue() const { return get(_maxvalRow, _maxvalCol); }
    inline UINT32 maxvalueRow() const { return _maxvalRow; }
    inline UINT32 maxvalueCol() const { return _maxvalCol; }

    inline double* getelementref(const UINT32& row, const UINT32& column) const {
        assert(_buffer != nullptr);
        assert(row < _mem_rows && "getelementref(): row out of bounds");
        assert(column < _mem_cols && "getelementref(): column out of bounds");
        if (_packed) {
            UINT32 i = row, j = column;
            if (i < j) { i = column; j = row; }
            return &_buffer[packed_index(_rows, i, j)];
        }
        if (_symmetric && row < column)
            return &(DNAMATRIX_ELEMENT(_buffer, _mem_rows, _mem_cols, column, row));
        return &(DNAMATRIX_ELEMENT(_buffer, _mem_rows, _mem_cols, row, column));
    }
    inline double* getelementref(const UINT32& row, const UINT32& column) {
        assert(_buffer != nullptr);
        assert(row < _mem_rows && "getelementref(): row out of bounds");
        assert(column < _mem_cols && "getelementref(): column out of bounds");
        if (_packed) {
            UINT32 i = row, j = column;
            if (i < j) { i = column; j = row; }
            return &_buffer[packed_index(_rows, i, j)];
        }
        if (_symmetric && row < column)
            return &(DNAMATRIX_ELEMENT(_buffer, _mem_rows, _mem_cols, column, row));
        return &(DNAMATRIX_ELEMENT(_buffer, _mem_rows, _mem_cols, row, column));
    }

    inline void mem_rows(const UINT32& r) { _mem_rows = r; }
    inline void mem_columns(const UINT32& c) { _mem_cols = c; }
    inline void rows(const UINT32& r) { _rows = r; }
    inline void columns(const UINT32& c) { _cols = c; }
    inline void maxvalueRow(const UINT32& r) { _maxvalRow = r; }
    inline void maxvalueCol(const UINT32& c) { _maxvalCol = c; }

    inline void put(const UINT32& row, const UINT32& column, const double& value) {
        assert(_buffer != nullptr);
        assert(row < _mem_rows && "put(): row out of bounds");
        assert(column < _mem_cols && "put(): column out of bounds");
        if (_packed) {
            UINT32 i = row, j = column;
            if (i < j) { i = column; j = row; }
            _buffer[packed_index(_rows, i, j)] = value;
            return;
        }
        DNAMATRIX_ELEMENT(_buffer, _mem_rows, _mem_cols, row, column) = value;
    }

    inline UINT32 matrixType() const { return _matrixType; }
    inline void matrixType(const UINT32 t) { _matrixType = t; }

    inline bool is_symmetric() const { return _symmetric; }
    inline void set_symmetric(bool s) {
        assert((!s || _rows == _cols) && "set_symmetric(true): matrix must be square");
        _symmetric = s;
    }

    inline bool is_packed() const { return _packed; }

    static inline std::size_t packed_index(UINT32 n, UINT32 i, UINT32 j) {
        return static_cast<std::size_t>(j) * n - static_cast<std::size_t>(j) * (j - 1) / 2 + (i - j);
    }

    static inline std::size_t packed_size(UINT32 n) {
        return static_cast<std::size_t>(n) * (n + 1) / 2;
    }

    // Matrix functions
    inline void copyelements(const UINT32& row_dest, const UINT32& column_dest, const matrix_2d& src,
                             const UINT32& row_src, const UINT32& column_src,
                             const UINT32& rows, const UINT32& columns) {
        assert(_buffer != nullptr && src._buffer != nullptr);
        assert(row_dest + rows <= _mem_rows && "copyelements(): dest row overflow");
        assert(column_dest + columns <= _mem_cols && "copyelements(): dest col overflow");
        assert(row_src + rows <= src._mem_rows && "copyelements(): src row overflow");
        assert(column_src + columns <= src._mem_cols && "copyelements(): src col overflow");
        if (rows == 3 && columns == 3 && !src._symmetric && !src._packed && !_packed) {
            for (UINT32 c = 0; c < 3; ++c) {
                double* dst = _buffer + static_cast<std::size_t>(column_dest + c) * _mem_rows + row_dest;
                const double* s = src._buffer + static_cast<std::size_t>(column_src + c) * src._mem_rows + row_src;
                dst[0] = s[0]; dst[1] = s[1]; dst[2] = s[2];
            }
            return;
        }
        copyelements_generic(row_dest, column_dest, src, row_src, column_src, rows, columns);
    }
    void copyelements_generic(const UINT32& row_dest, const UINT32& column_dest, const matrix_2d& src,
                              const UINT32& row_src, const UINT32& column_src,
                              const UINT32& rows, const UINT32& columns);
    inline void copyelements(const UINT32& row_dest, const UINT32& column_dest, const matrix_2d* src,
                             const UINT32& row_src, const UINT32& column_src,
                             const UINT32& rows, const UINT32& columns) {
        copyelements(row_dest, column_dest, *src, row_src, column_src, rows, columns);
    }

    inline void elementadd(const UINT32& row, const UINT32& column, const double& increment) {
        assert(row < _rows && column < _cols && "elementadd(): out of bounds");
        if (_packed && row < column) return;
        *getelementref(row, column) += increment;
    }

    inline void lower_add(const UINT32& row, const UINT32& column, const double& increment) {
        assert(_buffer != nullptr);
        assert(row < _rows && column < _cols && "lower_add(): out of bounds");
        if (_packed) {
            if (row < column) return;
            _buffer[packed_index(_rows, row, column)] += increment;
            return;
        }
        if (_symmetric && row < column) {
            _buffer[static_cast<std::size_t>(row) * _mem_rows + column] += increment;
            return;
        }
        _buffer[static_cast<std::size_t>(column) * _mem_rows + row] += increment;
    }

    inline void elementsubtract(const UINT32& row, const UINT32& column, const double& decrement) {
        assert(row < _rows && column < _cols && "elementsubtract(): out of bounds");
        if (_packed && row < column) return;
        *getelementref(row, column) -= decrement;
    }

    inline void elementmultiply(const UINT32& row, const UINT32& column, const double& scale) {
        assert(row < _rows && column < _cols && "elementmultiply(): out of bounds");
        *getelementref(row, column) *= scale;
    }

    inline void blockadd(const UINT32& row_dest, const UINT32& col_dest, const matrix_2d& mat_src,
                         const UINT32& row_src, const UINT32& col_src, const UINT32& rows, const UINT32& cols) {
        assert(_buffer != nullptr && mat_src._buffer != nullptr);
        assert(row_dest + rows <= _mem_rows && "blockadd(): dest row overflow");
        assert(col_dest + cols <= _mem_cols && "blockadd(): dest col overflow");
        assert(row_src + rows <= mat_src._mem_rows && "blockadd(): src row overflow");
        assert(col_src + cols <= mat_src._mem_cols && "blockadd(): src col overflow");
        if (rows == 3 && cols == 3 && !mat_src._symmetric && !mat_src._packed && !_packed) {
            for (UINT32 c = 0; c < 3; ++c) {
                double* dst = _buffer + static_cast<std::size_t>(col_dest + c) * _mem_rows + row_dest;
                const double* src = mat_src._buffer + static_cast<std::size_t>(col_src + c) * mat_src._mem_rows + row_src;
                dst[0] += src[0]; dst[1] += src[1]; dst[2] += src[2];
            }
            return;
        }
        blockadd_generic(row_dest, col_dest, mat_src, row_src, col_src, rows, cols);
    }
    void blockadd_generic(const UINT32& row_dest, const UINT32& col_dest, const matrix_2d& mat_src,
                          const UINT32& row_src, const UINT32& col_src, const UINT32& rows, const UINT32& cols);

    inline void blockTadd(const UINT32& row_dest, const UINT32& col_dest, const matrix_2d& mat_src,
                          const UINT32& row_src, const UINT32& col_src, const UINT32& rows, const UINT32& cols) {
        assert(_buffer != nullptr && mat_src._buffer != nullptr);
        assert(row_dest + rows <= _mem_rows && "blockTadd(): dest row overflow");
        assert(col_dest + cols <= _mem_cols && "blockTadd(): dest col overflow");
        assert(row_src + cols <= mat_src._mem_rows && "blockTadd(): src row overflow (transposed)");
        assert(col_src + rows <= mat_src._mem_cols && "blockTadd(): src col overflow (transposed)");
        if (rows == 3 && cols == 3 && !_packed) {
            for (UINT32 c = 0; c < 3; ++c) {
                double* dst = _buffer + static_cast<std::size_t>(col_dest + c) * _mem_rows + row_dest;
                const std::size_t smr = mat_src._mem_rows;
                const double* src_r0 = mat_src._buffer + static_cast<std::size_t>(row_src) * smr + col_src + c;
                dst[0] += src_r0[0];
                dst[1] += src_r0[smr];
                dst[2] += src_r0[2 * smr];
            }
            return;
        }
        blockTadd_generic(row_dest, col_dest, mat_src, row_src, col_src, rows, cols);
    }
    void blockTadd_generic(const UINT32& row_dest, const UINT32& col_dest, const matrix_2d& mat_src,
                           const UINT32& row_src, const UINT32& col_src, const UINT32& rows, const UINT32& cols);
    void blocksubtract(const UINT32& row_dest, const UINT32& col_dest, const matrix_2d& mat_src, const UINT32& row_src,
                       const UINT32& col_src, const UINT32& rows, const UINT32& cols);

    matrix_2d add(const matrix_2d& rhs);
    matrix_2d add(const matrix_2d& lhs, const matrix_2d& rhs);

    matrix_2d multiply(const char* lhs_trans, const matrix_2d& rhs, const char* rhs_trans); // multiplication
    matrix_2d multiply(const matrix_2d& lhs, const char* lhs_trans, const matrix_2d& rhs,
                       const char* rhs_trans); // multiplication
    // C = sym_lhs * rhs using dsymm (sym_lhs must be symmetric, lower triangle populated)
    matrix_2d multiply_sym(const matrix_2d& sym_lhs, const matrix_2d& rhs);

    matrix_2d sweepinverse();                                  // Sweep inverse (good for rotation matrices)
    matrix_2d cholesky_inverse(bool LOWER_IS_CLEARED = false, bool mark_symmetric = false); // Cholesky inverse
    matrix_2d cholesky_factor(bool LOWER_IS_CLEARED = false);  // Cholesky factor only (dpotrf), no inverse
    void cholesky_solve(matrix_2d& rhs, bool LOWER_IS_CLEARED = false);  // Solve from pre-factored matrix (dpotrs)
    double dot(const matrix_2d& other) const;  // Inner product of two column vectors

    matrix_2d transpose(const matrix_2d&); // Transpose
    matrix_2d transpose();                 //  ''
    matrix_2d scale(const double& scalar); // scale
    void scale_symmetric_diagonal(const double* diag);

    // overloaded operators
    // equality
    bool operator==(const matrix_2d& rhs) const {
        if (_mem_cols != rhs._mem_cols)
            return false;
        if (_mem_rows != rhs._mem_rows)
            return false;
        if (_cols != rhs._cols)
            return false;
        if (_rows != rhs._rows)
            return false;
        if (*_buffer != *rhs._buffer)
            return false;
        if (_maxvalCol != rhs._maxvalCol)
            return false;
        if (_maxvalRow != rhs._maxvalRow)
            return false;
        if (_matrixType != rhs._matrixType)
            return false;
        return true;
    }

    // equality
    bool operator!=(const matrix_2d& rhs) const {
        if (_mem_cols == rhs._mem_cols)
            return false;
        if (_mem_rows == rhs._mem_rows)
            return false;
        if (_cols == rhs._cols)
            return false;
        if (_rows == rhs._rows)
            return false;
        if (*_buffer == *rhs._buffer)
            return false;
        if (_maxvalCol == rhs._maxvalCol)
            return false;
        if (_maxvalRow == rhs._maxvalRow)
            return false;
        if (_matrixType == rhs._matrixType)
            return false;
        return true;
    }

    matrix_2d& operator=(const matrix_2d& rhs);
    matrix_2d operator*(const double& rhs) const;
    //
    // Initialisation / manipulation
    void allocate();
    void allocate(const UINT32& rows, const UINT32& columns);
    void setsize(const UINT32& rows,
                 const UINT32& columns); // sets matrix size to rows * columns only (buffer not allocated any memory)
    void redim(const UINT32& rows, const UINT32& columns); // redimensions matrix to rows * columns
    void redim_packed(const UINT32& n);                    // redimensions to packed symmetric n×n
    void replace(const UINT32& rowstart, const UINT32& columnstart, const matrix_2d& newmat);
    void replace(const UINT32& rowstart, const UINT32& columnstart, const UINT32& rows, const UINT32& columns,
                 const matrix_2d& newmat);
    void shrink(const UINT32& rows, const UINT32& columns);
    void grow(const UINT32& rows, const UINT32& columns);
    void clearlower(); // sets lower tri elements to zero
    void clearupper(); // sets upper tri elements to zero
    void filllower();  // copies upper tri to lower
    void fillupper();  // copies lower tri to upper
    void zero();       // sets all elements to zero
    void zero(const UINT32& row_begin, const UINT32& col_begin, const UINT32& rows, const UINT32& columns);
    double compute_maximum_value();

    // Printing
    friend std::ostream& operator<<(std::ostream& os, const matrix_2d& rhs);

    // Reading from memory mapped file
    void ReadMappedFileRegion(void* addr);

    // Writing to memory mapped file
    void WriteMappedFileRegion(void* addr);

    // In-place mmap: attach _buffer directly to mmap data region (no memcpy)
    void AttachMappedFileRegion(void* addr);

    // In-place mmap: write footer back, detach _buffer
    void DetachMappedFileRegion(void* addr);

    inline bool owns_buffer() const { return _owns_buffer; }

    // debug
#ifdef _MSDEBUG
    void trace(const std::string& comment, const std::string& format) const;
    void trace(const std::string& comment, const std::string& submat_comment, const std::string& format,
               const UINT32& row_begin, const UINT32& col_begin, const UINT32& rows, const UINT32& columns) const;
#endif

    void deallocate();

  private:
    inline std::size_t buffersize() const {
        if (_packed)
            return packed_size(_mem_rows) * sizeof(double);
        return static_cast<std::size_t>(_mem_rows) * static_cast<std::size_t>(_mem_cols) * sizeof(double);
    }
    void buy(const UINT32& rows, const UINT32& columns, double** mem_space);
    void copybuffer(const UINT32& rows, const UINT32& columns, const matrix_2d& oldmat);
    void copybuffer(const UINT32& rowstart, const UINT32& columnstart, const UINT32& rows, const UINT32& columns,
                    const matrix_2d& mat);
    // void copybuffer(const UINT32& rows, const UINT32& columns, double**	buffer);

    void sweep(UINT32 k1, UINT32 k2);

    UINT32 _mem_cols; // actual buffer size (cols)
    UINT32 _mem_rows; // actual buffer size (rows)
    UINT32 _cols;     // number of actual cols
    UINT32 _rows;     // number of actual rows
    double* _buffer;  // matrix buffer elements
    bool _owns_buffer; // true if _buffer is heap-allocated (false when attached to mmap)

    UINT32 _maxvalCol; // col of max value
    UINT32 _maxvalRow; // row of max value

    UINT32 _matrixType; // full, upper/lower, sparse
    bool _symmetric;    // only lower triangle populated (upper is mirror)
    bool _packed;       // packed column-major lower-triangle storage (n*(n+1)/2 elements)
};

} // namespace math
} // namespace dynadjust

#endif // DNAMATRIX_CONTIGUOUS_H_
