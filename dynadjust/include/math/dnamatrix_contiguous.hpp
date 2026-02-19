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
void LAPACK_FUNC(dsytrf)(const char* uplo, const lapack_int* n, double* a, const lapack_int* lda,
                         lapack_int* ipiv, double* work, const lapack_int* lwork, lapack_int* info);
void LAPACK_FUNC(dsytri)(const char* uplo, const lapack_int* n, double* a, const lapack_int* lda,
                         const lapack_int* ipiv, double* work, lapack_int* info);
void BLAS_FUNC(dgemm)(const enum CBLAS_ORDER ORDER, const enum CBLAS_TRANSPOSE TRANSA,
                      const enum CBLAS_TRANSPOSE TRANSB, const lapack_int M, const lapack_int N, const lapack_int K,
                      const double ALPHA, const double* A, const lapack_int LDA, const double* B, const lapack_int LDB,
                      const double BETA, double* C, const lapack_int LDC);
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

    // element retrieval
    // see DNAMATRIX_ROW_WISE
    inline double& get(const UINT32& row, const UINT32& column) const {
        return DNAMATRIX_ELEMENT(_buffer, _mem_rows, _mem_cols, row, column);
    }
    inline double* getbuffer(const UINT32& row, const UINT32& column) const {
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
        return &(DNAMATRIX_ELEMENT(_buffer, _mem_rows, _mem_cols, row, column));
    }
    inline double* getelementref(const UINT32& row, const UINT32& column) {
        return &(DNAMATRIX_ELEMENT(_buffer, _mem_rows, _mem_cols, row, column));
    }

    inline void mem_rows(const UINT32& r) { _mem_rows = r; }
    inline void mem_columns(const UINT32& c) { _mem_cols = c; }
    inline void rows(const UINT32& r) { _rows = r; }
    inline void columns(const UINT32& c) { _cols = c; }
    inline void maxvalueRow(const UINT32& r) { _maxvalRow = r; }
    inline void maxvalueCol(const UINT32& c) { _maxvalCol = c; }

    inline void put(const UINT32& row, const UINT32& column, const double& value) {
        DNAMATRIX_ELEMENT(_buffer, _mem_rows, _mem_cols, row, column) = value;
    }

    inline UINT32 matrixType() const { return _matrixType; }
    inline void matrixType(const UINT32 t) { _matrixType = t; }

    // Matrix functions
    void copyelements(const UINT32& row_dest, const UINT32& column_dest, const matrix_2d& src, const UINT32& row_src,
                      const UINT32& column_src, const UINT32& rows, const UINT32& columns);
    void copyelements(const UINT32& row_dest, const UINT32& column_dest, const matrix_2d* src, const UINT32& row_src,
                      const UINT32& column_src, const UINT32& rows, const UINT32& columns);

    inline void elementadd(const UINT32& row, const UINT32& column, const double& increment) {
        *getelementref(row, column) += increment;
    }

    inline void elementsubtract(const UINT32& row, const UINT32& column, const double& decrement) {
        *getelementref(row, column) -= decrement;
    }

    inline void elementmultiply(const UINT32& row, const UINT32& column, const double& scale) {
        *getelementref(row, column) *= scale;
    }

    void blockadd(const UINT32& row_dest, const UINT32& col_dest, const matrix_2d& mat_src, const UINT32& row_src,
                  const UINT32& col_src, const UINT32& rows, const UINT32& cols);
    void blockTadd(const UINT32& row_dest, const UINT32& col_dest, const matrix_2d& mat_src, const UINT32& row_src,
                   const UINT32& col_src, const UINT32& rows, const UINT32& cols);
    void blocksubtract(const UINT32& row_dest, const UINT32& col_dest, const matrix_2d& mat_src, const UINT32& row_src,
                       const UINT32& col_src, const UINT32& rows, const UINT32& cols);

    matrix_2d add(const matrix_2d& rhs);
    matrix_2d add(const matrix_2d& lhs, const matrix_2d& rhs);

    matrix_2d multiply(const char* lhs_trans, const matrix_2d& rhs, const char* rhs_trans); // multiplication
    matrix_2d multiply(const matrix_2d& lhs, const char* lhs_trans, const matrix_2d& rhs,
                       const char* rhs_trans); // multiplication

    matrix_2d sweepinverse();                                  // Sweep inverse (good for rotation matrices)
    matrix_2d cholesky_inverse(bool LOWER_IS_CLEARED = false); // Cholesky inverse

    matrix_2d transpose(const matrix_2d&); // Transpose
    matrix_2d transpose();                 //  ''
    matrix_2d scale(const double& scalar); // scale

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

    // debug
#ifdef _MSDEBUG
    void trace(const std::string& comment, const std::string& format) const;
    void trace(const std::string& comment, const std::string& submat_comment, const std::string& format,
               const UINT32& row_begin, const UINT32& col_begin, const UINT32& rows, const UINT32& columns) const;
#endif

  private:
    inline std::size_t buffersize() const {
        return static_cast<std::size_t>(_mem_rows) * static_cast<std::size_t>(_mem_cols) * sizeof(double);
    }

    void deallocate();
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

    UINT32 _maxvalCol; // col of max value
    UINT32 _maxvalRow; // row of max value

    UINT32 _matrixType; // full, upper/lower, sparse
};

} // namespace math
} // namespace dynadjust

#endif // DNAMATRIX_CONTIGUOUS_H_
