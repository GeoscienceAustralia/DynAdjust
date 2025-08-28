//============================================================================
// Name         : dnamatrix_contiguous.cpp
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

#include <include/math/dnamatrix_contiguous.hpp>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <limits>
#include <include/ide/trace.hpp>

namespace dynadjust {
namespace math {

std::ostream& operator<<(std::ostream& os, const matrix_2d& rhs) {
    if (os.iword(0) == binary) {
        // Binary output

        // matrix type
        os.write(reinterpret_cast<const char*>(&rhs._matrixType), sizeof(UINT32));

        // output rows and columns
        os.write(reinterpret_cast<const char*>(&rhs._rows), sizeof(UINT32));
        os.write(reinterpret_cast<const char*>(&rhs._cols), sizeof(UINT32));

        os.write(reinterpret_cast<const char*>(&rhs._mem_rows), sizeof(UINT32));
        os.write(reinterpret_cast<const char*>(&rhs._mem_cols), sizeof(UINT32));

        UINT32 c, r;

        switch (rhs._matrixType) {
        case mtx_lower:
            // output lower triangular part of a square matrix
            if (rhs._mem_rows != rhs._mem_cols)
                throw std::runtime_error("matrix_2d operator<< (): Matrix is not square.");

            // print each column
            for (c = 0; c < rhs._mem_cols; ++c)
                os.write(reinterpret_cast<const char*>(rhs.getelementref(c, c)), (rhs._mem_rows - c) * sizeof(double));

            break;
        case mtx_sparse: break;
        case mtx_full:
        default:
            // Output full matrix data
            for (r = 0; r < rhs._mem_rows; ++r)
                for (c = 0; c < rhs._mem_cols; ++c)
                    os.write(reinterpret_cast<const char*>(rhs.getelementref(r, c)), sizeof(double));
            break;
        }

        // output max value info
        os.write(reinterpret_cast<const char*>(&rhs._maxvalRow), sizeof(UINT32));
        os.write(reinterpret_cast<const char*>(&rhs._maxvalCol), sizeof(UINT32));
    } else {
        // ASCII output

        os << rhs._matrixType << " " << rhs._rows << " " << rhs._cols << " " << rhs._mem_rows << " " << rhs._mem_cols
           << std::endl;

        for (UINT32 c, r = 0; r < rhs._mem_rows; ++r) {
            for (c = 0; c < rhs._mem_cols; ++c)
                os << std::scientific << std::setprecision(16) << rhs.get(r, c) << " ";
            os << std::endl;
        }
        os << rhs._maxvalRow << " " << rhs._maxvalCol << std::endl;
        os << std::endl;
    }
    return os;
}

UINT32 __row__;
UINT32 __col__;
// string _method_;

void out_of_memory_handler() {
    std::size_t mem = static_cast<std::size_t>(__row__) * static_cast<std::size_t>(__col__) * sizeof(double);

    std::stringstream ss;
    ss << "Insufficient memory available to create a " << __row__ << " x " << __col__ << " matrix (" << std::fixed
       << std::setprecision(2);

    if (mem < MEGABYTE_SIZE)
        ss << (mem / KILOBYTE_SIZE) << " KB).";
    else if (mem < GIGABYTE_SIZE)
        ss << (mem / MEGABYTE_SIZE) << " MB).";
    else // if (mem >= GIGABYTE_SIZE)
        ss << (mem / GIGABYTE_SIZE) << " GB).";

    throw NetMemoryException(ss.str());
}

matrix_2d::matrix_2d()
    : _mem_cols(0), _mem_rows(0), _cols(0), _rows(0), _buffer(0), _maxvalCol(0), _maxvalRow(0), _matrixType(mtx_full) {
    std::set_new_handler(out_of_memory_handler);

    // if this class were to be modified to use templates, each
    // instance could be tested for an invalid data type as follows
    //
    // if (strcmp(typeid(a(1,1)).name(), "double") != 0 &&
    //	strcmp(typeid(a(1,1)).name(), "float") != 0 ) {
    //	throw std::runtime_error("Not a floating point type");
}

matrix_2d::matrix_2d(const UINT32& rows, const UINT32& columns)
    : _mem_cols(columns), _mem_rows(rows), _cols(columns), _rows(rows), _buffer(0), _maxvalCol(0), _maxvalRow(0),
      _matrixType(mtx_full) {
    std::set_new_handler(out_of_memory_handler);

    allocate(_rows, _cols);
}

matrix_2d::matrix_2d(const UINT32& rows, const UINT32& columns, const double data[], const std::size_t& data_size,
                     const UINT32& matrix_type)
    : _mem_cols(columns), _mem_rows(rows), _cols(columns), _rows(rows), _buffer(0), _maxvalCol(0), _maxvalRow(0),
      _matrixType(matrix_type) {
    std::set_new_handler(out_of_memory_handler);

    std::stringstream ss;
    std::size_t upperMatrixElements(sumOfConsecutiveIntegers(rows));
    UINT32 j;

    const double* dataptr = &data[0];

    switch (matrix_type) {
    case mtx_lower:
        // Lower triangular part of a square matrix
        if (upperMatrixElements != data_size) {
            ss << "Data size must be equivalent to upper matrix element count for " << rows << " x " << columns << ".";
            throw std::runtime_error(ss.str());
        }

        // Create memory and store the data
        allocate(_rows, _cols);

        for (j = 0; j < columns; ++j) {
            memcpy(getelementref(j, j), dataptr, (static_cast<std::size_t>(rows) - j) * sizeof(double));
            dataptr += (static_cast<std::size_t>(rows) - j);
        }

        fillupper();
        break;

    case mtx_sparse:
        ss << "matrix_2d(): A sparse matrix cannot be initialised with a double array.";
        throw std::runtime_error(ss.str());
        break;

    case mtx_full:
    default:
        // Full matrix
        if (data_size != static_cast<std::size_t>(rows) * static_cast<std::size_t>(columns)) {
            ss << "Data size must be equivalent to matrix dimensions (" << rows << " x " << columns << ").";
            throw std::runtime_error(ss.str());
        }

        // Create memory and store the data
        allocate(_rows, _cols);
        memcpy(_buffer, data, data_size * sizeof(double));

        break;
    }
}

matrix_2d::matrix_2d(const matrix_2d& newmat)
    : _mem_cols(newmat.memColumns()), _mem_rows(newmat.memRows()), _cols(newmat.columns()), _rows(newmat.rows()),
      _buffer(0), _maxvalCol(newmat.maxvalueCol()), _maxvalRow(newmat.maxvalueRow()), _matrixType(newmat.matrixType()) {
    std::set_new_handler(out_of_memory_handler);

    allocate(_mem_rows, _mem_cols);

    const double* ptr = newmat.getbuffer();

    // copy buffer
    memcpy(_buffer, ptr, newmat.buffersize());
}

matrix_2d::~matrix_2d() {
    // Default destructor
    deallocate();
}

std::size_t matrix_2d::get_size() {
    size_t size =
        (7 * sizeof(UINT32)); // UINT32 _matrixType, _mem_cols, _mem_rows, _cols, _rows, _maxvalRow, _maxvalCol

    switch (_matrixType) {
    case mtx_lower: size += sumOfConsecutiveIntegers(_mem_rows) * sizeof(double); break;
    case mtx_sparse: break;
    case mtx_full:
    default: size += buffersize();
    }
    return size;
}

// Read data from memory mapped file
void matrix_2d::ReadMappedFileRegion(void* addr) {
    // IMPORTANT
    // The following read statements must correspond
    // with that which is written in operator>> below.

    PUINT32 data_U = reinterpret_cast<PUINT32>(addr);
    _matrixType = *data_U++;
    _rows = *data_U++;
    _cols = *data_U++;

    switch (_matrixType) {
    case mtx_sparse:
        // _mem_cols and _mem_rows equal _cols and _rows
        _mem_rows = _rows;
        _mem_cols = _cols;
        break;
    case mtx_lower:
    case mtx_full:
    default:
        _mem_rows = *data_U++;
        _mem_cols = *data_U++;
        break;
    }

    allocate(_mem_rows, _mem_cols);

    double* data_d;
    int* data_i;

    UINT32 c, r, i;
    int ci;

    switch (_matrixType) {
    case mtx_sparse:
        data_i = reinterpret_cast<int*>(data_U);
        for (r = 0; r < _rows; ++r) {
            // A row corresponding to stations 1, 2 and 3
            for (i = 0; i < 3; ++i) {
                // get column index of first element
                ci = *data_i++;

                // get elements
                data_d = reinterpret_cast<double*>(data_i);

                if (ci < 0) {
                    data_d += 3;
                    data_i = reinterpret_cast<int*>(data_d);
                    continue;
                }

                memcpy(getelementref(r, ci), data_d, sizeof(double)); // xValue
                data_d++;
                memcpy(getelementref(r, ci + 1), data_d, sizeof(double)); // yValue
                data_d++;
                memcpy(getelementref(r, ci + 2), data_d, sizeof(double)); // zValue
                data_d++;

                data_i = reinterpret_cast<int*>(data_d);
            }
        }
        return;
        break;
    case mtx_lower:
        data_d = reinterpret_cast<double*>(data_U);

        // read each column
        for (c = 0; c < _mem_cols; ++c) {
            memcpy(getelementref(c, c), data_d, (_mem_rows - c) * sizeof(double));
            data_d += (_mem_rows - c);
        }

        fillupper();
        break;
    case mtx_full:
    default:
        data_d = reinterpret_cast<double*>(data_U);

        // get contiguous block from memory
        memcpy(_buffer, data_d, buffersize());
        // skip to UINT32 elements
        data_d += _mem_rows * _mem_cols;
        break;
    }

    // Get UINT pointer
    data_U = reinterpret_cast<UINT32*>(data_d);

    _maxvalRow = *data_U++;
    _maxvalCol = *data_U;
}

// Write data to memory mapped file
void matrix_2d::WriteMappedFileRegion(void* addr) {
    // IMPORTANT
    // The following write statements must correspond
    // with that which is written in operator<< above.

    PUINT32 data_U = reinterpret_cast<UINT32*>(addr);
    *data_U++ = _matrixType;
    *data_U++ = _rows;
    *data_U++ = _cols;

    switch (_matrixType) {
    case mtx_sparse:
        // _mem_cols and _mem_rows aren't written
        // because for sparse matrices they are the
        // same size as _cols and _rows
        break;
    case mtx_lower:
    case mtx_full:
    default:
        *data_U++ = _mem_rows;
        *data_U++ = _mem_cols;
        break;
    }

    double* data_d;
    int* data_i;

    UINT32 c, r, i;
    int ci;

    switch (_matrixType) {
    case mtx_sparse:
        data_i = reinterpret_cast<int*>(data_U);
        for (r = 0; r < _rows; ++r) {
            // A row corresponding to stations 1, 2 and 3
            for (i = 0; i < 3; ++i) {
                // get column index of first element
                ci = *data_i++;

                // write elements
                data_d = reinterpret_cast<double*>(data_i);

                if (ci < 0) {
                    data_d += 3;
                    data_i = reinterpret_cast<int*>(data_d);
                    continue;
                }

                memcpy(data_d, getelementref(r, ci), sizeof(double)); // xValue
                data_d++;
                memcpy(data_d, getelementref(r, ci + 1), sizeof(double)); // yValue
                data_d++;
                memcpy(data_d, getelementref(r, ci + 2), sizeof(double)); // zValue
                data_d++;

                data_i = reinterpret_cast<int*>(data_d);
            }
        }
        return;
        break;
    case mtx_lower:
        data_d = reinterpret_cast<double*>(data_U);

        // print each column
        for (c = 0; c < _mem_cols; ++c) {
            memcpy(data_d, getbuffer(c, c), (_mem_rows - c) * sizeof(double));
            data_d += (_mem_rows - c);
        }
        break;
    case mtx_full:
    default:
        data_d = reinterpret_cast<double*>(data_U);

        // write contiguous block to memory
        memcpy(data_d, _buffer, _mem_rows * _mem_cols * sizeof(double));
        // skip to UINT32 elements
        data_d += _mem_rows * _mem_cols;
        break;
    }

    // Get UINT pointer
    data_U = reinterpret_cast<UINT32*>(data_d);

    *data_U++ = _maxvalRow;
    *data_U = _maxvalCol;
}

void matrix_2d::allocate() { allocate(_mem_rows, _mem_cols); }

// creates memory for desired "memory size", not matrix dimensions
void matrix_2d::allocate(const UINT32& rows, const UINT32& columns) {
    //_method_ = "allocate";

    deallocate();

    // an exception will be thrown by out_of_memory_handler
    // if memory cannot be allocated
    buy(rows, columns, &_buffer);
}

// creates memory for desired "memory size", not matrix dimensions
void matrix_2d::buy(const UINT32& rows, const UINT32& columns, double** mem_space) {
    //_method_ = "buy";

    // set globals for new_memory_handler function
    __row__ = rows;
    __col__ = columns;

    // an exception will be thrown by out_of_memory_handler
    // if memory cannot be allocated
    std::size_t total_size = static_cast<std::size_t>(rows) * static_cast<std::size_t>(columns);
    (*mem_space) = new double[total_size];

    if ((*mem_space) == nullptr) {
        std::stringstream ss;
        ss << "Insufficient memory for a " << rows << " x " << columns << " matrix.";
        throw NetMemoryException(ss.str());
    }
    
    // Initialize memory to zero to prevent uninitialized values
    std::memset((*mem_space), 0, total_size * sizeof(double));
}

void matrix_2d::deallocate() {
    if (_buffer != nullptr) {
        delete[] _buffer;
        _buffer = nullptr;
    }
}

matrix_2d matrix_2d::submatrix(const UINT32& row_begin, const UINT32& col_begin, const UINT32& rows,
                               const UINT32& columns) const {
    matrix_2d b(rows, columns);
    submatrix(row_begin, col_begin, &b, rows, columns);
    return b;
}

void matrix_2d::submatrix(const UINT32& row_begin, const UINT32& col_begin, matrix_2d* dest, const UINT32& subrows,
                          const UINT32& subcolumns) const {
    if (row_begin >= _rows) {
        std::stringstream ss;
        ss << row_begin << ", " << col_begin << " lies outside the range of the matrix (" << _rows << ", " << _cols
           << ").";
        throw std::runtime_error(ss.str());
    }

    if (col_begin >= _cols) {
        std::stringstream ss;
        ss << row_begin << ", " << col_begin << " lies outside the range of the matrix (" << _rows << ", " << _cols
           << ").";
        throw std::runtime_error(ss.str());
    }

    if (subrows > dest->rows()) {
        std::stringstream ss;
        ss << subrows << ", " << subcolumns << " exceeds the size of the matrix (" << dest->rows() << ", "
           << dest->columns() << ").";
        throw std::runtime_error(ss.str());
    }

    if (subcolumns > dest->columns()) {
        std::stringstream ss;
        ss << subrows << ", " << subcolumns << " exceeds the size of the matrix (" << dest->rows() << ", "
           << dest->columns() << ").";
        throw std::runtime_error(ss.str());
    }

    if (row_begin + subrows > _rows) {
        std::stringstream ss;
        ss << row_begin + subrows << ", " << col_begin + subcolumns << " lies outside the range of the matrix ("
           << _rows << ", " << _cols << ").";
        throw std::runtime_error(ss.str());
    }

    if (col_begin + subcolumns > _cols) {
        std::stringstream ss;
        ss << row_begin + subrows << ", " << col_begin + subcolumns << " lies outside the range of the matrix ("
           << _rows << ", " << _cols << ").";
        throw std::runtime_error(ss.str());
    }

    UINT32 i, j, m(0), n(0), row_end(row_begin + subrows), col_end(col_begin + subcolumns);
    for (i = row_begin; i < row_end; ++i) {
        for (j = col_begin; j < col_end; ++j) {
            dest->put(m, n, get(i, j));
            n++;
        }
        n = 0;
        m++;
    }
}

void matrix_2d::redim(const UINT32& rows, const UINT32& columns) {
    // if new matrix size is smaller than or equal to the previous
    // matrix size, then simply change dimensions and return
    if (rows <= _mem_rows && columns <= _mem_cols) {
        _rows = rows;
        _cols = columns;
        return;
    }

    //_method_ = "redim";

    // double* new_buffer;
    // buy(rows, columns, &new_buffer);
    // copybuffer(_rows, _cols, &new_buffer);
    // deallocate();
    //_buffer = new_buffer;

    std::set_new_handler(out_of_memory_handler);

    deallocate();
    allocate(rows, columns);

    _rows = _mem_rows = rows;
    _cols = _mem_cols = columns;
}

void matrix_2d::shrink(const UINT32& rows, const UINT32& columns) {
    if (rows > _rows || columns > _cols) {
        std::stringstream ss;
        ss << " " << std::endl;
        if (rows >= _rows)
            ss << "    Cannot shrink by " << rows << " rows on a matrix of " << _rows << " rows. " << std::endl;
        if (columns >= _cols)
            ss << "    Cannot shrink by " << columns << " columns on a matrix of " << _cols << " columns.";
        throw std::runtime_error(ss.str());
    }

    _rows -= rows;
    _cols -= columns;
}

void matrix_2d::grow(const UINT32& rows, const UINT32& columns) {
    if ((rows + _rows) > _mem_rows || (columns + _cols) > _mem_cols) {
        std::stringstream ss;
        ss << " " << std::endl;
        if (rows >= _rows)
            ss << "    Cannot grow matrix by " << rows << " rows: growth exceeds row memory limit (" << _mem_rows
               << ").";
        if (columns >= _cols)
            ss << "    Cannot grow matrix by " << columns << " columns: growth exceeds column memory limit ("
               << _mem_cols << ").";
        throw std::runtime_error(ss.str());
    }

    _rows += rows;
    _cols += columns;
}

void matrix_2d::setsize(const UINT32& rows, const UINT32& columns) {
    deallocate();
    _rows = _mem_rows = rows;
    _cols = _mem_cols = columns;
}

void matrix_2d::replace(const UINT32& rowstart, const UINT32& columnstart, const matrix_2d& newmat) {
    copybuffer(rowstart, columnstart, newmat.rows(), newmat.columns(), newmat);
}

void matrix_2d::copybuffer(const UINT32& rows, const UINT32& columns, const matrix_2d& oldmat) {
    if (rows == _mem_rows && columns == _mem_cols) {
        memcpy(_buffer, oldmat.getbuffer(), buffersize());
        return;
    }

    UINT32 column;
    for (column = 0; column < columns; column++) {
        memcpy(getelementref(0, column), oldmat.getbuffer(0, column), static_cast<std::size_t>(rows) * sizeof(double));
    }
}

void matrix_2d::copybuffer(const UINT32& rowstart, const UINT32& columnstart, const UINT32& rows, const UINT32& columns,
                           const matrix_2d& mat) {
    UINT32 rowend(rowstart + rows), columnend(columnstart + columns);
    if (rowend > _rows || columnend > _cols) {
        std::stringstream ss;
        ss << " " << std::endl;
        if (rowend >= _rows)
            ss << "    Row index " << rowend << " exceeds the matrix row count (" << _rows << "). " << std::endl;
        if (columnend >= _cols)
            ss << "    Column index " << columnend << " exceeds the matrix column count (" << _cols << ").";
        throw std::runtime_error(ss.str());
    }

    UINT32 column(0), c(0);
    for (column = columnstart; column < columnend; ++column, ++c) {
        memcpy(getelementref(rowstart, column), mat.getbuffer(0, c), static_cast<std::size_t>(rows) * sizeof(double));
    }
}

void matrix_2d::copyelements(const UINT32& row_dest, const UINT32& column_dest, const matrix_2d& src,
                             const UINT32& row_src, const UINT32& column_src, const UINT32& rows,
                             const UINT32& columns) {
    UINT32 cd(0), cs(0), colend_dest(column_dest + columns);
    for (cd = column_dest, cs = column_src; cd < colend_dest; ++cd, ++cs)
        memcpy(getelementref(row_dest, cd), src.getbuffer(row_src, cs),
               static_cast<std::size_t>(rows) * sizeof(double));
}

void matrix_2d::copyelements(const UINT32& row_dest, const UINT32& column_dest, const matrix_2d* src,
                             const UINT32& row_src, const UINT32& column_src, const UINT32& rows,
                             const UINT32& columns) {
    copyelements(row_dest, column_dest, *src, row_src, column_src, rows, columns);
}

void matrix_2d::sweep(UINT32 k1, UINT32 k2) {
    double eps(1.0e-8), d;
    UINT32 i, j, k, it;

    if (k2 < k1) {
        k = k1;
        k1 = k2;
        k2 = k;
    }
    //	n = a.nrows();
    for (k = k1; k < k2; k++)      //	for (k = k1; k <= k2; k++)
    {                              //	{
        if (fabs(get(k, k)) < eps) //		if ( fabs( a(k, k) ) < eps)
        {
            for (it = 0; it < _rows; it++) //			for (it = 1; it <= n; it++)
            {
                put(it, k, 0.);
                put(k, it, 0.); //				a(it, k) = a(k, it) = 0.0;
            }
        } else {                                  //		else {
            d = 1.0 / get(k, k);                  //			d = 1.0 / a(k, k);
            put(k, k, d);                         //			a(k, k) = d;
            for (i = 0; i < _rows; i++)           //			for (i = 1; i <= n; i++)
                if (i != k)                       //				if (i != k)
                    *getelementref(i, k) *= -d;   //					a(i, k) *= (T) - d;
            for (j = 0; j < _rows; j++)           //			for (j = 1; j <= n; j++)
                if (j != k)                       //				if (j != k)
                    *getelementref(k, j) *= d;    //					a(k, j) *= (T) d;
            for (i = 0; i < _rows; i++) {         //			for (i = 1; i <= n; i++) {
                if (i != k) {                     //				if (i != k) {
                    for (j = 0; j < _rows; j++) { //					for (j = 1; j <= n; j++) {
                        if (j != k)               //						if (j != k)
                            *getelementref(i, j) +=
                                get(i, k) * get(k, j) / d; //							a(i, j) += a(i, k) *a(k, j) / d;
                    } // end for j													//					} // end for j
                } // end for i != k													//				} // end for i != k
            } // end for i															//			} // end for i
        } // end else																//		} // end else
    } // end for k																	//	} // end for k
}

matrix_2d matrix_2d::sweepinverse() {
    if (_rows != _cols)
        throw std::runtime_error("sweepinverse: Matrix is not square.");

    sweep(0, _rows);
    return *this;
}

matrix_2d matrix_2d::cholesky_inverse(bool LOWER_IS_CLEARED /*=false*/) {
    if (_rows < 1)
        return *this;

    if (_rows != _cols)
        throw std::runtime_error("cholesky_inverse(): Matrix is not square.");

    char uplo(LOWER_TRIANGLE);

    // Which triangle is filled - upper or lower?
    if (LOWER_IS_CLEARED)
        uplo = UPPER_TRIANGLE;

    lapack_int info, n = _rows;

    // Perform Cholesky factorisation
    LAPACK_FUNC(dpotrf)(&uplo, &n, _buffer, &n, &info);

    if (info != 0)
        throw std::runtime_error("choleskyinverse_mkl(): Cholesky factorisation failed.");

    // Perform Cholesky inverse
    LAPACK_FUNC(dpotri)(&uplo, &n, _buffer, &n, &info);

    if (info != 0)
        throw std::runtime_error("choleskyinverse_mkl: Cholesky inversion failed.");

    // Copy empty triangle part
    if (LOWER_IS_CLEARED)
        filllower();
    else
        fillupper();

    return *this;
}

// scale()
//
matrix_2d matrix_2d::scale(const double& scalar) {
    UINT32 i, j;
    for (i = 0; i < _rows; ++i)
        for (j = 0; j < _cols; ++j)
            *getelementref(i, j) *= scalar;
    return *this;
}

void matrix_2d::blockadd(const UINT32& row_dest, const UINT32& col_dest, const matrix_2d& mat_src,
                         const UINT32& row_src, const UINT32& col_src, const UINT32& rows, const UINT32& cols) {
    UINT32 i_dest, j_dest, i_src, j_src;
    UINT32 i_dest_end(row_dest + rows), j_dest_end(col_dest + cols);

    for (i_dest = row_dest, i_src = row_src; i_dest < i_dest_end; ++i_dest, ++i_src)
        for (j_dest = col_dest, j_src = col_src; j_dest < j_dest_end; ++j_dest, ++j_src)
            elementadd(i_dest, j_dest, mat_src.get(i_src, j_src));
}

// Same as blockadd, but adds transpose.  mat_src must be square.
void matrix_2d::blockTadd(const UINT32& row_dest, const UINT32& col_dest, const matrix_2d& mat_src,
                          const UINT32& row_src, const UINT32& col_src, const UINT32& rows, const UINT32& cols) {
    UINT32 i_dest, j_dest, i_src, j_src;
    UINT32 i_dest_end(row_dest + rows), j_dest_end(col_dest + cols);

    for (i_dest = row_dest, i_src = row_src; i_dest < i_dest_end; ++i_dest, ++i_src)
        for (j_dest = col_dest, j_src = col_src; j_dest < j_dest_end; ++j_dest, ++j_src)
            elementadd(i_dest, j_dest, mat_src.get(j_src, i_src));
}

void matrix_2d::blocksubtract(const UINT32& row_dest, const UINT32& col_dest, const matrix_2d& mat_src,
                              const UINT32& row_src, const UINT32& col_src, const UINT32& rows, const UINT32& cols) {
    UINT32 i_dest, j_dest, i_src, j_src;
    UINT32 i_dest_end(row_dest + rows), j_dest_end(col_dest + cols);

    for (i_dest = row_dest, i_src = row_src; i_dest < i_dest_end; ++i_dest, ++i_src)
        for (j_dest = col_dest, j_src = col_src; j_dest < j_dest_end; ++j_dest, ++j_src)
            elementsubtract(i_dest, j_dest, mat_src.get(i_src, j_src));
}

// clearlower()
void matrix_2d::clearlower() {
    // Sets lower triangle elements to zero
    UINT32 col, row;
    for (row = 1, col = 0; col < _mem_cols; ++col, ++row)
        memset(getelementref(row, col), 0, (static_cast<std::size_t>(_mem_rows) - row) * sizeof(double));
}

// filllower()
void matrix_2d::filllower() {
    // copies upper triangle to lower triangle
    UINT32 column, row;
    for (row = 1; row < _rows; row++)
        for (column = 0; column < row; column++)
            put(row, column, get(column, row));
}

// fillupper()
void matrix_2d::fillupper() {
    // copies lower triangle to upper triangle
    UINT32 column, row;
    for (row = 1; row < _rows; row++)
        for (column = 0; column < row; column++)
            put(column, row, get(row, column));
}

// zero()
void matrix_2d::zero() { memset(_buffer, 0, buffersize()); }

// zero()
void matrix_2d::zero(const UINT32& row_begin, const UINT32& col_begin, const UINT32& rows, const UINT32& columns) {

    UINT32 col(0), col_end(col_begin + columns);
    for (col = col_begin; col < col_end; ++col)
        memset(getelementref(row_begin, col), 0, rows * sizeof(double));
}

matrix_2d matrix_2d::operator=(const matrix_2d& rhs) {
    // Overloaded assignment operator
    if (this == &rhs)
        return *this;

    // If rhs data can fit within limits of this matrix, copy
    // and return. Otherwise, allocate new memory
    if (_mem_rows >= rhs.rows() || _mem_cols >= rhs.columns()) {
        // don't change _mem_rows or _mem_cols.  Simply update
        // visible dimensions and copy buffer
        _rows = rhs.rows();
        _cols = rhs.columns();
        copybuffer(_rows, _cols, rhs);

        _maxvalCol = rhs.maxvalueCol(); // col of max value
        _maxvalRow = rhs.maxvalueRow(); // row of max value

        return *this;
    }

    // Okay, rhs is larger, so allocate new memory. Call free
    // memory first before changing row and column dimensions!
    deallocate();
    _mem_rows = rhs.memRows(); // change memory limits
    _mem_cols = rhs.memColumns();
    _rows = rhs.rows(); // change matrix dimensions
    _cols = rhs.columns();
    allocate(_mem_rows, _mem_cols);
    copybuffer(_rows, _cols, rhs);

    _maxvalCol = rhs.maxvalueCol(); // col of max value
    _maxvalRow = rhs.maxvalueRow(); // row of max value

    return *this;
}

// Multiplication operator
matrix_2d matrix_2d::operator*(const double& rhs) const {
    // Answer
    matrix_2d m(_rows, _cols);

    UINT32 row, column;
    for (row = 0; row < _rows; row++)
        for (column = 0; column < _cols; ++column)
            m.put(row, column, get(row, column) * rhs);
    return m;
}

matrix_2d matrix_2d::add(const matrix_2d& rhs) {
    if (_rows != rhs.rows() || _cols != rhs.columns())
        throw std::runtime_error("add: Result matrix dimensions are incompatible.");

    UINT32 row, column;
    for (row = 0; row < _rows; row++) {
        for (column = 0; column < _cols; ++column) {
            *getelementref(row, column) += rhs.get(row, column);
        }
    }
    return *this;
}

// multiplies this matrix by rhs and stores the result in a new matrix
// Uses Intel MKL dgemm
matrix_2d matrix_2d::multiply(const char* lhs_trans, const matrix_2d& rhs, const char* rhs_trans) {
    matrix_2d m(_rows, rhs.columns());

    const double one = 1.0;
    const double zero = 0.0;

    lapack_int lhs_rows(rows()), rhs_cols(rhs.columns());
    lapack_int lhs_cols(columns()), rhs_rows(rhs.rows());
    lapack_int new_mem_rows(memRows());
    lapack_int lhs_mem_rows(memRows());
    lapack_int rhs_mem_rows(memRows());

    if (strcmp(lhs_trans, "T") == 0) {
        lhs_rows = columns(); // transpose
        lhs_cols = rows();    // transpose
    }

    if (strcmp(rhs_trans, "T") == 0) {
        rhs_rows = rhs.columns(); // transpose
        rhs_cols = rhs.rows();    // transpose
    }

    if (lhs_cols != rhs_rows)
        throw std::runtime_error("multiply_mkl(): Matrix dimensions are incompatible.");
    else if (_rows != lhs_rows || _cols != rhs_cols)
        throw std::runtime_error("multiply_mkl(): Result matrix dimensions are incompatible.");

    CBLAS_TRANSPOSE tA = (strcmp(lhs_trans, "T") == 0) ? CblasTrans : CblasNoTrans;
    CBLAS_TRANSPOSE tB = (strcmp(rhs_trans, "T") == 0) ? CblasTrans : CblasNoTrans;

    BLAS_FUNC(dgemm)(CblasColMajor, tA, tB, lhs_rows, rhs_cols, lhs_cols, 1.0, _buffer, _mem_rows, rhs.getbuffer(),
                     rhs.memRows(), 0.0, m.getbuffer(), m.memRows());

    return (*this = m);
}

// Multiplies lhs by rhs and stores the result in this.
// Uses Intel MKL dgemm
matrix_2d
matrix_2d::multiply(const matrix_2d& lhs, const char* lhs_trans, const matrix_2d& rhs, const char* rhs_trans) {
    const double one = 1.0;
    const double zero = 0.0;

    lapack_int lhs_rows(lhs.rows()), rhs_cols(rhs.columns());
    lapack_int lhs_cols(lhs.columns()), rhs_rows(rhs.rows());
    lapack_int new_mem_rows(memRows());
    lapack_int lhs_mem_rows(lhs.memRows());
    lapack_int rhs_mem_rows(rhs.memRows());

    if (strncmp(lhs_trans, "T", 1) == 0) {
        lhs_rows = lhs.columns(); // transpose
        lhs_cols = lhs.rows();    // transpose
    }

    if (strncmp(rhs_trans, "T", 1) == 0) {
        rhs_rows = rhs.columns(); // transpose
        rhs_cols = rhs.rows();    // transpose
    }

    if (lhs_cols != rhs_rows)
        throw std::runtime_error("multiply: Matrix dimensions are incompatible.");
    else if (_rows != lhs_rows || _cols != rhs_cols)
        throw std::runtime_error("multiply(): Result matrix dimensions are incompatible.");

    CBLAS_TRANSPOSE tA = (strncmp(lhs_trans, "T", 1) == 0) ? CblasTrans : CblasNoTrans;
    CBLAS_TRANSPOSE tB = (strncmp(rhs_trans, "T", 1) == 0) ? CblasTrans : CblasNoTrans;

    BLAS_FUNC(dgemm)(CblasColMajor, tA, tB, lhs_rows, rhs_cols, lhs_cols, 1.0, lhs.getbuffer(), lhs.memRows(),
                     rhs.getbuffer(), rhs.memRows(), 0.0, _buffer, _mem_rows);

    return *this;
} // Multiply()

// Transpose()
matrix_2d matrix_2d::transpose(const matrix_2d& matA) {
    if ((matA.columns() != _rows) || (matA.rows() != _cols))
        throw std::runtime_error("transpose: Matrix dimensions are incompatible.");

    UINT32 column, row;
    for (row = 0; row < _rows; row++)
        for (column = 0; column < _cols; column++)
            *getelementref(row, column) = matA.get(column, row);
    return *this;
} // Transpose()

// Transpose()
matrix_2d matrix_2d::transpose() {
    matrix_2d m(_cols, _rows);
    UINT32 column, row;
    for (row = 0; row < _rows; row++)
        for (column = 0; column < _cols; column++)
            m.put(column, row, get(row, column));
    return m;
} // Transpose()

// computes and retains the maximum value in the matrix
double matrix_2d::compute_maximum_value() {
    _maxvalCol = _maxvalRow = 0;
    UINT32 col, row;
    for (row = 0; row < _rows; ++row) {
        for (col = 0; col < _cols; col++) {
            if (fabs(get(row, col)) > fabs(get(_maxvalRow, _maxvalCol))) {
                _maxvalCol = col;
                _maxvalRow = row;
            }
        }
    }
    return get(_maxvalRow, _maxvalCol);
}

#ifdef _MSDEBUG
void matrix_2d::trace(const std::string& comment, const std::string& format) const {
    UINT32 i, j;
    if (comment.empty())
        TRACE("%d %d\n", _rows, _cols);
    else
        TRACE("%s (%d, %d):\n", comment.c_str(), _rows, _cols);
    for (i = 0; i < _rows; ++i) {
        for (j = 0; j < _cols; ++j)
            TRACE(format.c_str(), get(i, j));
        TRACE("\n");
    }
    TRACE("\n");
}
#endif

#ifdef _MSDEBUG
void matrix_2d::trace(const std::string& comment, const std::string& submat_comment, const std::string& format,
                      const UINT32& row_begin, const UINT32& col_begin, const UINT32& rows,
                      const UINT32& columns) const {
    // comparison of unsigned expression < 0 is always false
    if (row_begin >= _rows) {
        TRACE("%d %d lies outside the range of the matrix (%d %d)\n", row_begin, col_begin, _rows, _cols);
        return;
    }
    // comparison of unsigned expression < 0 is always false
    if (col_begin >= _cols) {
        TRACE("%d %d lies outside the range of the matrix (%d %d)\n", row_begin, col_begin, _rows, _cols);
        return;
    }
    if (row_begin + rows > _rows) {
        TRACE("%d %d lies outside the range of the matrix (%d %d)\n", row_begin, col_begin, _rows, _cols);
        return;
    }
    if (col_begin + columns > _cols) {
        TRACE("%d %d lies outside the range of the matrix (%d %d)\n", row_begin, col_begin, _rows, _cols);
        return;
    }

    if (comment.empty())
        TRACE("%d %d, %s submatrix (%d, %d, %d*%d)\n", _rows, _cols, submat_comment.c_str(), row_begin, col_begin, rows,
              columns);
    else
        TRACE("%s (%d, %d), %s submatrix (%d, %d, %d*%d)\n", comment, _rows, _cols, submat_comment.c_str(), row_begin,
              col_begin, rows, columns);

    UINT32 i, j, row_end(row_begin + rows), col_end(col_begin + columns);

    for (i = row_begin; i < row_end; ++i) {
        for (j = col_begin; j < col_end; ++j)
            TRACE(format.c_str(), get(i, j));
        TRACE("\n");
    }
    TRACE("\n");
}
#endif

} // namespace math
} // namespace dynadjust
