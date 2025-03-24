#ifndef DNAMATRIX_CONTIGUOUS_H_
#define DNAMATRIX_CONTIGUOUS_H_

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

#if defined(_MSC_VER)
#if defined(LIST_INCLUDES_ON_BUILD)
#pragma message("  " __FILE__)
#endif
#endif

#ifdef _MSDEBUG
#include <include/ide/trace.hpp>
#endif

#if defined(__APPLE__)

/* Use the new LAPACK interface on Apple by default */
#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>

/* Define these missing constants on Apple */
#define LAPACK_COL_MAJOR 102
#define LAPACK_ROW_MAJOR 101

/* The ILP64 interface can be enabled by defining ACCELERATE_LAPACK_ILP64 */
#if defined(ACCELERATE_LAPACK_ILP64)
typedef int64_t lapack_int;
#else
typedef int lapack_int;
#endif

/* We need to define the LAPACKE interface for Apple  as it currently does not
 * exist */

inline int LAPACKE_dpotrf(int layout, char uplo, int n, double* a, int lda) {
    int info = 0;
    if (layout != CblasColMajor) {
        return -1;
    }
    dpotrf_(&uplo, &n, a, &lda, &info);
    return info;
}

inline int LAPACKE_dpotri(int layout, char uplo, int n, double* a, int lda) {
    int info = 0;
    if (layout != CblasColMajor) {
        return -1;
    }
    dpotri_(&uplo, &n, a, &lda, &info);
    return info;
}

#elif defined(MKL_ILP64) // Intel MKL with ILP64

#include <mkl.h>
#include <mkl_lapacke.h>
typedef MKL_INT64 lapack_int;

#elif defined(MKL_LP64) // Intel MKL with LP64

#include <mkl.h>
#include <mkl_lapacke.h>
typedef MKL_INT lapack_int;

#else // Fallback to standard LAPACKE

#include <cblas.h>
#include <lapacke.h>
typedef int lapack_int;

#endif

#include <cstring>
#include <include/config/dnatypes.hpp>
#include <include/functions/dnatemplatecalcfuncs.hpp>
#include <include/memory/dnamemory_handler.hpp>
#include <iostream>

// #include <include/config/dnaconsts.hpp>
// #include <include/config/dnaversion.hpp>

#include <include/exception/dnaexception.hpp>

using namespace dynadjust::memory;
using namespace dynadjust::exception;

namespace dynadjust {
namespace math {

// Forward declarations
class matrix_2d;
typedef std::vector<matrix_2d> v_mat_2d, *pv_mat_2d;
typedef v_mat_2d::iterator _it_v_mat_2d;
typedef v_mat_2d::const_iterator _it_v_mat_2d_const;
typedef std::vector<v_mat_2d> vv_mat_2d;

/*
  By default uses column-major order
*/
#define DNAMATRIX_INDEX(num_rows, num_cols, row, column)                       \
    ((column) * (num_rows) + (row))

/*
  Retrieves an element based on DNAMATRIX_INDEX based indexing
*/
#define DNAMATRIX_ELEMENT(A, no_rows, no_cols, row, column)                    \
    A[DNAMATRIX_INDEX(no_rows, no_cols, row, column)]

/*
  byteSize template: calculates memory size for 'elements' of type T
*/
template <typename T> std::size_t byteSize(const index_t elements = 1) {
    return elements * sizeof(T);
}

/*
  A matrix implementation using column-major storage.
*/
class matrix_2d : public new_handler_support<matrix_2d> {
  public:
    /* Default constructor: initialises an empty matrix */
    matrix_2d();

    /* Constructor: creates matrix with given rows and columns */
    matrix_2d(const index_t& rows, const index_t& columns);

    /*
       Constructor: creates matrix from data array, can handle full/lower forms
       matrix_type default = mtx_full
    */
    matrix_2d(const index_t& rows, const index_t& columns, const double data[],
              const index_t& data_size, const index_t& matrix_type = mtx_full);

    /* Copy constructor: makes a deep copy of newmat */
    matrix_2d(const matrix_2d&);

    /* Destructor: deallocates matrix buffer */
    ~matrix_2d();

    /* Checks if matrix buffer is null */
    inline bool empty() { return _buffer == nullptr; }

    /* Returns size of matrix in bytes, depends on matrixType */
    std::size_t get_size();

    /* Returns total allocated memory rows */
    inline index_t memRows() const { return _mem_rows; }

    /* Returns total allocated memory columns */
    inline index_t memColumns() const { return _mem_cols; }

    /* Returns actual rows in matrix */
    inline index_t rows() const { return _rows; }

    /* Returns actual columns in matrix */
    inline index_t columns() const { return _cols; }

    /* Returns pointer to underlying buffer */
    inline double* getbuffer() const { return _buffer; }

    /* Returns reference to element at row r, column c */
    inline double& get(const index_t& r, const index_t& c) const {
        return DNAMATRIX_ELEMENT(_buffer, _mem_rows, _mem_cols, r, c);
    }

    /* Returns pointer to element at row r, column c */
    inline double* getbuffer(const index_t& r, const index_t& c) const {
        return _buffer + DNAMATRIX_INDEX(_mem_rows, _mem_cols, r, c);
    }

    /*
       Extracts a submatrix of size rows x columns
       starting at (row_begin, col_begin) into 'dest'
    */
    void submatrix(const index_t& row_begin, const index_t& col_begin,
                   matrix_2d* dest, const index_t& rows,
                   const index_t& columns) const;

    /* Returns a submatrix of size rows x columns from (row_begin, col_begin) */
    matrix_2d submatrix(const index_t& row_begin, const index_t& col_begin,
                        const index_t& rows, const index_t& columns) const;

    /* Returns the current maximum value in the matrix */
    inline double maxvalue() const { return get(_maxvalRow, _maxvalCol); }

    /* Returns the row index of the current maximum value */
    inline index_t maxvalueRow() const { return _maxvalRow; }

    /* Returns the column index of the current maximum value */
    inline index_t maxvalueCol() const { return _maxvalCol; }

    /* Returns pointer to the element reference at row r, column c */
    inline double* getelementref(const index_t& r, const index_t& c) const {
        return &(DNAMATRIX_ELEMENT(_buffer, _mem_rows, _mem_cols, r, c));
    }

    /* Returns pointer to the element reference at row r, column c (non-const) */
    inline double* getelementref(const index_t& r, const index_t& c) {
        return &(DNAMATRIX_ELEMENT(_buffer, _mem_rows, _mem_cols, r, c));
    }

    /* Sets the total allocated rows (memory) */
    inline void mem_rows(const index_t& r) { _mem_rows = r; }

    /* Sets the total allocated columns (memory) */
    inline void mem_columns(const index_t& c) { _mem_cols = c; }

    /* Sets the actual number of rows */
    inline void rows(const index_t& r) { _rows = r; }

    /* Sets the actual number of columns */
    inline void columns(const index_t& c) { _cols = c; }

    /* Sets the row index for the maximum value */
    inline void maxvalueRow(const index_t& r) { _maxvalRow = r; }

    /* Sets the column index for the maximum value */
    inline void maxvalueCol(const index_t& c) { _maxvalCol = c; }

    /* Places val into row r, column c */
    inline void put(const index_t& r, const index_t& c, const double& val) {
        DNAMATRIX_ELEMENT(_buffer, _mem_rows, _mem_cols, r, c) = val;
    }

    /* Returns current matrix type */
    inline index_t matrixType() const { return _matrixType; }

    /* Sets the current matrix type */
    inline void matrixType(const index_t t) { _matrixType = t; }

    /*
       Copies rows x columns of src from (row_src, column_src)
       into this matrix at (row_dest, column_dest)
    */
    void copyelements(const index_t& row_dest, const index_t& column_dest,
                      const matrix_2d& src, const index_t& row_src,
                      const index_t& column_src, const index_t& rows,
                      const index_t& columns);

    /* Overload: same as above, but src is a pointer */
    void copyelements(const index_t& row_dest, const index_t& column_dest,
                      const matrix_2d* src, const index_t& row_src,
                      const index_t& column_src, const index_t& rows,
                      const index_t& columns);

    /* Adds 'inc' to the element at (r, c) */
    inline void
    elementadd(const index_t& r, const index_t& c, const double& inc) {
        *getelementref(r, c) += inc;
    }

    /* Subtracts 'dec' from the element at (r, c) */
    inline void
    elementsubtract(const index_t& r, const index_t& c, const double& dec) {
        *getelementref(r, c) -= dec;
    }

    /* Multiplies the element at (r, c) by 'scale' */
    inline void
    elementmultiply(const index_t& r, const index_t& c, const double& scale) {
        *getelementref(r, c) *= scale;
    }

    /*
       Adds a sub-block of mat_src at (row_src, col_src)
       to this matrix at (row_dest, col_dest)
    */
    void
    blockadd(const index_t& row_dest, const index_t& col_dest,
             const matrix_2d& mat_src, const index_t& row_src,
             const index_t& col_src, const index_t& rows, const index_t& cols);

    /*
       Adds the transpose of a sub-block of mat_src
       to this matrix at (row_dest, col_dest)
    */
    void
    blockTadd(const index_t& row_dest, const index_t& col_dest,
              const matrix_2d& mat_src, const index_t& row_src,
              const index_t& col_src, const index_t& rows, const index_t& cols);

    /*
       Subtracts a sub-block of mat_src at (row_src, col_src)
       from this matrix at (row_dest, col_dest)
    */
    void blocksubtract(const index_t& row_dest, const index_t& col_dest,
                       const matrix_2d& mat_src, const index_t& row_src,
                       const index_t& col_src, const index_t& rows,
                       const index_t& cols);

    /*
       In-place add of rhs to this matrix (dimensions must match)
    */
    matrix_2d add(const matrix_2d& rhs);

    /*
       Static add: sets this = lhs + rhs, returns result
       (dimensions must match)
    */
    matrix_2d add(const matrix_2d& lhs, const matrix_2d& rhs);

    /*
       Matrix multiplication using MKL (or fallback):
       multiplies *this by rhs
       lhs_trans or rhs_trans can be "T" or "N"
    */
    matrix_2d multiply_mkl(const char* lhs_trans, const matrix_2d& rhs,
                           const char* rhs_trans);

    /*
       Matrix multiplication using MKL (or fallback):
       stores lhs*rhs in *this
    */
    matrix_2d multiply_mkl(const matrix_2d& lhs, const char* lhs_trans,
                           const matrix_2d& rhs, const char* rhs_trans);

    /*
       Performs a sweep operation on the entire matrix in-place
       to find the inverse (requires square matrix)
    */
    matrix_2d sweepinverse();

    /*
       In-place Cholesky inversion using MKL (or fallback)
       if LOWER_IS_CLEARED=true, lower part is zeroed before
    */
    matrix_2d cholesky_inverse(bool LOWER_IS_CLEARED = false);

    /*
       In-place transpose: sets *this = transpose(matA)
       (Requires matching dimensions for 'matA' and *this)
    */
    matrix_2d transpose(const matrix_2d&);

    /* Returns a new transposed copy of *this */
    matrix_2d transpose();

    /* Scales matrix elements by 'scalar' in-place */
    matrix_2d scale(const double& scalar);

    /* Equality operator */
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

    /* Inequality operator */
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

    /* Copy assignment operator */
    matrix_2d operator=(const matrix_2d& rhs);

    /* Scalar multiplication: returns a new matrix_2d(*this * rhs) */
    matrix_2d operator*(const double& rhs) const;

    /* Allocates memory for _rows, _cols (in-place) */
    void allocate();

    /* Allocates memory for given rows, columns (in-place) */
    void allocate(const index_t& rows, const index_t& columns);

    /*
       Sets new matrix dimension (rows, columns) and
       deallocates old buffer
    */
    void setsize(const index_t& rows, const index_t& columns);

    /*
       Re-dimensions the matrix, possibly creating new buffer
       if current memory is insufficient
    */
    void redim(const index_t& rows, const index_t& columns);

    /*
       Replaces a sub-block at (rowstart, columnstart)
       with contents of newmat
    */
    void replace(const index_t& rowstart, const index_t& columnstart,
                 const matrix_2d& newmat);

    /*
       Overload: replaces (rows x columns) from newmat
       into *this at (rowstart, columnstart)
    */
    void replace(const index_t& rowstart, const index_t& columnstart,
                 const index_t& rows, const index_t& columns,
                 const matrix_2d& newmat);

    /* Shrinks matrix dimensions by (rows, columns) without reallocation */
    void shrink(const index_t& rows, const index_t& columns);

    /* Grows matrix dimensions by (rows, columns) without reallocation */
    void grow(const index_t& rows, const index_t& columns);

    /* Sets lower triangle to zero */
    void clearlower();

    /* Copies upper to lower triangle (makes symmetrical) */
    void filllower();

    /* Copies lower to upper triangle (makes symmetrical) */
    void fillupper();

    /* Sets entire matrix to zero */
    void zero();

    /* Zeros a sub-block of size (rows, columns) from (row_begin, col_begin) */
    void zero(const index_t& row_begin, const index_t& col_begin,
              const index_t& rows, const index_t& columns);

    /* Finds maximum value in matrix and updates _maxvalRow/_maxvalCol */
    double compute_maximum_value();

    /* Overloaded << operator for matrix output */
    friend std::ostream& operator<<(std::ostream& os, const matrix_2d& rhs);

    /* Reads matrix data from mapped file region */
    void ReadMappedFileRegion(void* addr);

    /* Writes matrix data to mapped file region */
    void WriteMappedFileRegion(void* addr);

#ifdef _MSDEBUG
    /* Debug trace of entire matrix */
    void trace(const std::string& comment, const std::string& format) const;

    /* Debug trace of submatrix */
    void trace(const std::string& comment, const std::string& submat_comment,
               const std::string& format, const index_t& row_begin,
               const index_t& col_begin, const index_t& rows,
               const index_t& columns) const;
#endif

  private:
    /* Returns size in bytes of buffer */
    inline std::size_t buffersize() const {
        return _mem_rows * _mem_cols * sizeof(double);
    }

    /* Frees the matrix buffer */
    void deallocate();

    /* Allocates buffer for rows x columns into mem_space */
    void buy(const index_t& rows, const index_t& columns, double** mem_space);

    /* Copies oldmat into this matrix (size must match) */
    void copybuffer(const index_t& rows, const index_t& columns,
                    const matrix_2d& oldmat);

    /* Copies rows x columns of mat into this matrix at rowstart, columnstart */
    void copybuffer(const index_t& rowstart, const index_t& columnstart,
                    const index_t& rows, const index_t& columns,
                    const matrix_2d& mat);

    /* Sweep function for matrix inversion */
    void sweep(index_t k1, index_t k2);

    index_t _cols;             // actual columns
    index_t _rows;             // actual rows
    index_t _mem_cols;         // memory columns
    index_t _mem_rows;         // memory rows
    index_t _maxvalCol;        // column index of maximum value
    index_t _maxvalRow;        // row index of maximum value
    index_t _matrixType;       // matrix type (full, lower, sparse) - not used?
    double* _buffer = nullptr; // pointer to actual matrix memory
};

} // namespace math
} // namespace dynadjust

#endif // DNAMATRIX_CONTIGUOUS_H_
