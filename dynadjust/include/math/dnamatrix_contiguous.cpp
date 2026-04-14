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

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <include/ide/trace.hpp>
#include <include/math/dnamatrix_contiguous.hpp>
#include <iomanip>
#include <sstream>

namespace dynadjust {
namespace math {

// Maximum BLAS thread count, set via set_max_blas_threads().
static int g_max_blas_threads = 0;

void set_max_blas_threads(int n) { g_max_blas_threads = n; }
int get_max_blas_threads() { return g_max_blas_threads; }

// Compute BLAS thread count scaled proportionally to g_max_blas_threads.
// Tiers: max, max/2, max/4, 1 — selected by matrix dimension n.
static inline int blas_threads_for(lapack_int n) {
    int max_t = g_max_blas_threads;
    if (max_t <= 1) return 1;
    if (n > 8000) return max_t;
    if (n > 4000) return std::max(1, max_t / 2);
    if (n > 2000) return std::max(1, max_t / 4);
    return 1;
}

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

        // Alignment padding — keeps data region at 8-byte aligned offset (24 bytes)
        // Must match WriteMappedFileRegion / ReadMappedFileRegion / AttachMappedFileRegion layout
        const UINT32 pad = 0;
        os.write(reinterpret_cast<const char*>(&pad), sizeof(UINT32));

        UINT32 c, r;

        switch (rhs._matrixType) {
        case mtx_lower:
            if (rhs._mem_rows != rhs._mem_cols)
                throw std::runtime_error("matrix_2d operator<< (): Matrix is not square.");

            if (rhs._packed) {
                os.write(reinterpret_cast<const char*>(rhs._buffer),
                         matrix_2d::packed_size(rhs._mem_rows) * sizeof(double));
            } else {
                for (c = 0; c < rhs._mem_cols; ++c)
                    os.write(reinterpret_cast<const char*>(rhs.getelementref(c, c)), (rhs._mem_rows - c) * sizeof(double));
            }
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
            for (c = 0; c < rhs._mem_cols; ++c) os << std::scientific << std::setprecision(16) << rhs.get(r, c) << " ";
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
    else  // if (mem >= GIGABYTE_SIZE)
        ss << (mem / GIGABYTE_SIZE) << " GB).";

    throw NetMemoryException(ss.str());
}

matrix_2d::matrix_2d()
    : _mem_cols(0), _mem_rows(0), _cols(0), _rows(0), _buffer(0), _owns_buffer(true), _maxvalCol(0), _maxvalRow(0), _matrixType(mtx_full), _symmetric(false), _packed(false) {
    std::set_new_handler(out_of_memory_handler);

    // if this class were to be modified to use templates, each
    // instance could be tested for an invalid data type as follows
    //
    // if (strcmp(typeid(a(1,1)).name(), "double") != 0 &&
    //	strcmp(typeid(a(1,1)).name(), "float") != 0 ) {
    //	throw std::runtime_error("Not a floating point type");
}

matrix_2d::matrix_2d(const UINT32& rows, const UINT32& columns)
    : _mem_cols(columns),
      _mem_rows(rows),
      _cols(columns),
      _rows(rows),
      _buffer(0),
      _owns_buffer(true),
      _maxvalCol(0),
      _maxvalRow(0),
      _matrixType(mtx_full),
      _symmetric(false),
      _packed(false) {
    std::set_new_handler(out_of_memory_handler);

    allocate(_rows, _cols);
}

matrix_2d::matrix_2d(const UINT32& rows, const UINT32& columns, const double data[], const std::size_t& data_size,
                     const UINT32& matrix_type)
    : _mem_cols(columns),
      _mem_rows(rows),
      _cols(columns),
      _rows(rows),
      _buffer(0),
      _owns_buffer(true),
      _maxvalCol(0),
      _maxvalRow(0),
      _matrixType(matrix_type),
      _symmetric(false),
      _packed(false) {
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
    : _mem_cols(newmat.memColumns()),
      _mem_rows(newmat.memRows()),
      _cols(newmat.columns()),
      _rows(newmat.rows()),
      _buffer(0),
      _owns_buffer(true),
      _maxvalCol(newmat.maxvalueCol()),
      _maxvalRow(newmat.maxvalueRow()),
      _matrixType(newmat.matrixType()),
      _symmetric(newmat._symmetric),
      _packed(newmat._packed) {
    std::set_new_handler(out_of_memory_handler);

    if (_packed) {
        std::size_t ps = packed_size(_mem_rows);
        _buffer = static_cast<double*>(std::malloc(ps * sizeof(double)));
        if (!_buffer) throw NetMemoryException("Insufficient memory for packed matrix copy.");
        memcpy(_buffer, newmat.getbuffer(), ps * sizeof(double));
    } else {
        // Allocate without zeroing — memcpy immediately overwrites entire buffer
        std::size_t total_size = static_cast<std::size_t>(_mem_rows) * _mem_cols;
        _buffer = static_cast<double*>(std::malloc(total_size * sizeof(double)));
        if (!_buffer) throw NetMemoryException("Insufficient memory for matrix copy.");
        memcpy(_buffer, newmat.getbuffer(), total_size * sizeof(double));
    }
}

matrix_2d::~matrix_2d() {
    // Default destructor
    deallocate();
}

std::size_t matrix_2d::get_size() {
    // 8 UINT32s: matrixType, rows, cols, mem_rows, mem_cols, _pad, maxvalRow, maxvalCol
    // The padding UINT32 ensures the data region starts at an 8-byte aligned offset (24 bytes)
    // from the region base, enabling in-place mmap buffer attachment.
    size_t size =
        (8 * sizeof(UINT32));

    switch (_matrixType) {
    case mtx_lower: size += sumOfConsecutiveIntegers(_mem_rows) * sizeof(double); break;
    case mtx_sparse: break;
    case mtx_full:
    default:
        if (_packed)
            size += packed_size(_mem_rows) * sizeof(double);
        else
            size += buffersize();
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
        _mem_rows = _rows;
        _mem_cols = _cols;
        break;
    case mtx_lower:
    case mtx_full:
    default:
        _mem_rows = *data_U++;
        _mem_cols = *data_U++;
        ++data_U;  // skip alignment padding UINT32
        break;
    }

    if (_matrixType == mtx_lower && _mem_rows == _mem_cols) {
        deallocate();
        _packed = true;
        _symmetric = true;
        std::size_t ps = packed_size(_mem_rows);
        _buffer = static_cast<double*>(std::calloc(ps, sizeof(double)));
        if (!_buffer) throw NetMemoryException("Insufficient memory for packed matrix read.");
    } else {
        _packed = false;
        allocate(_mem_rows, _mem_cols);
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

                // get elements
                data_d = reinterpret_cast<double*>(data_i);

                if (ci < 0) {
                    data_d += 3;
                    data_i = reinterpret_cast<int*>(data_d);
                    continue;
                }

                memcpy(getelementref(r, ci), data_d, sizeof(double));  // xValue
                data_d++;
                memcpy(getelementref(r, ci + 1), data_d, sizeof(double));  // yValue
                data_d++;
                memcpy(getelementref(r, ci + 2), data_d, sizeof(double));  // zValue
                data_d++;

                data_i = reinterpret_cast<int*>(data_d);
            }
        }
        return;
        break;
    case mtx_lower:
        assert(_mem_rows == _mem_cols && "ReadMappedFileRegion(mtx_lower): matrix must be square");
        data_d = reinterpret_cast<double*>(data_U);

        if (_packed) {
            std::size_t ps = packed_size(_mem_rows);
            memcpy(_buffer, data_d, ps * sizeof(double));
            data_d += ps;
        } else {
            for (c = 0; c < _mem_cols; ++c) {
                memcpy(getbuffer(c, c), data_d, (_mem_rows - c) * sizeof(double));
                data_d += (_mem_rows - c);
            }
        }

        _symmetric = true;
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
        *data_U++ = 0;  // alignment padding
        break;
    }

    // In-place mmap: data is already in the region, only update footer
    if (!_owns_buffer) {
        // Calculate footer position by skipping past data
        double* data_d;
        switch (_matrixType) {
        case mtx_lower:
            if (_packed)
                data_d = reinterpret_cast<double*>(data_U) + packed_size(_mem_rows);
            else
                data_d = reinterpret_cast<double*>(data_U) + sumOfConsecutiveIntegers(_mem_rows);
            break;
        case mtx_full:
        default:
            data_d = reinterpret_cast<double*>(data_U) + static_cast<std::size_t>(_mem_rows) * _mem_cols;
            break;
        }
        PUINT32 footer = reinterpret_cast<UINT32*>(data_d);
        *footer++ = _maxvalRow;
        *footer = _maxvalCol;
        return;
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

                memcpy(data_d, getelementref(r, ci), sizeof(double));  // xValue
                data_d++;
                memcpy(data_d, getelementref(r, ci + 1), sizeof(double));  // yValue
                data_d++;
                memcpy(data_d, getelementref(r, ci + 2), sizeof(double));  // zValue
                data_d++;

                data_i = reinterpret_cast<int*>(data_d);
            }
        }
        return;
        break;
    case mtx_lower:
        data_d = reinterpret_cast<double*>(data_U);

        if (_packed) {
            std::size_t ps = packed_size(_mem_rows);
            memcpy(data_d, _buffer, ps * sizeof(double));
            data_d += ps;
        } else {
            for (c = 0; c < _mem_cols; ++c) {
                memcpy(data_d, getbuffer(c, c), (_mem_rows - c) * sizeof(double));
                data_d += (_mem_rows - c);
            }
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


void matrix_2d::AttachMappedFileRegion(void* addr) {
    // Read header metadata — same layout as ReadMappedFileRegion
    PUINT32 data_U = reinterpret_cast<PUINT32>(addr);
    _matrixType = *data_U++;
    _rows = *data_U++;
    _cols = *data_U++;

    switch (_matrixType) {
    case mtx_sparse:
        _mem_rows = _rows;
        _mem_cols = _cols;
        // Sparse cannot use in-place — fall back to copy path
        ReadMappedFileRegion(addr);
        return;
    case mtx_lower:
    case mtx_full:
    default:
        _mem_rows = *data_U++;
        _mem_cols = *data_U++;
        ++data_U;  // skip alignment padding
        break;
    }

    // Release any existing buffer
    deallocate();

    // Point _buffer directly at the data region in the mmap
    _buffer = reinterpret_cast<double*>(data_U);
    _owns_buffer = false;

    // Set packed/symmetric flags for lower-triangular matrices
    if (_matrixType == mtx_lower && _mem_rows == _mem_cols) {
        _packed = true;
        _symmetric = true;
    } else {
        _packed = false;
    }

    // Read footer (maxvalRow, maxvalCol) from after the data
    double* data_d;
    switch (_matrixType) {
    case mtx_lower:
        data_d = _buffer + packed_size(_mem_rows);
        break;
    case mtx_full:
    default:
        data_d = _buffer + static_cast<std::size_t>(_mem_rows) * _mem_cols;
        break;
    }

    PUINT32 footer = reinterpret_cast<UINT32*>(data_d);
    _maxvalRow = *footer++;
    _maxvalCol = *footer;
}


void matrix_2d::DetachMappedFileRegion(void* addr) {
    if (_owns_buffer || _buffer == nullptr)
        return;

    // Write updated footer (maxvalRow, maxvalCol) to the mmap region.
    // Header: 6 UINT32s (matrixType, rows, cols, mem_rows, mem_cols, pad) = 24 bytes.
    // Skip past data to find the footer position.
    double* data_d;
    switch (_matrixType) {
    case mtx_lower:
        data_d = _buffer + packed_size(_mem_rows);
        break;
    case mtx_full:
    default:
        data_d = _buffer + static_cast<std::size_t>(_mem_rows) * _mem_cols;
        break;
    }

    PUINT32 footer = reinterpret_cast<UINT32*>(data_d);
    *footer++ = _maxvalRow;
    *footer = _maxvalCol;

    // Detach without freeing the mmap memory
    _buffer = nullptr;
    _owns_buffer = true;
}


void matrix_2d::allocate() {
    if (_matrixType == mtx_lower && _mem_rows == _mem_cols && _mem_rows > 0) {
        deallocate();
        _packed = true;
        _symmetric = true;
        std::size_t ps = packed_size(_mem_rows);
        _buffer = static_cast<double*>(std::calloc(ps, sizeof(double)));
        if (!_buffer) {
            std::stringstream ss;
            ss << "Insufficient memory for a packed " << _mem_rows << " x " << _mem_rows << " matrix.";
            throw NetMemoryException(ss.str());
        }
        return;
    }
    allocate(_mem_rows, _mem_cols);
}

// creates memory for desired "memory size", not matrix dimensions
void matrix_2d::allocate(const UINT32& rows, const UINT32& columns) {
    deallocate();
    buy(rows, columns, &_buffer);
}

// creates memory for desired "memory size", not matrix dimensions
void matrix_2d::buy(const UINT32& rows, const UINT32& columns, double** mem_space) {
    //_method_ = "buy";

    // set globals for new_memory_handler function
    __row__ = rows;
    __col__ = columns;

    // calloc returns zero-initialised memory.  For large allocations glibc
    // satisfies calloc via mmap, whose anonymous pages are already zero-filled
    // by the kernel — so calloc skips the redundant memset that was previously
    // triggering page-fault storms (97 % of CPU time on the NSW benchmark).
    std::size_t total_size = static_cast<std::size_t>(rows) * static_cast<std::size_t>(columns);
    (*mem_space) = static_cast<double*>(std::calloc(total_size, sizeof(double)));

    if ((*mem_space) == nullptr) {
        std::stringstream ss;
        ss << "Insufficient memory for a " << rows << " x " << columns << " matrix.";
        throw NetMemoryException(ss.str());
    }
}

void matrix_2d::deallocate() {
    if (_buffer != nullptr) {
        if (_owns_buffer)
            std::free(_buffer);
        _buffer = nullptr;
    }
    _owns_buffer = true;
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
    if (_packed) {
        deallocate();
        _packed = false;
    }
    _symmetric = false;
    // if new matrix size is smaller than or equal to the previous
    // matrix size, then simply change dimensions and return
    if (_buffer != nullptr && rows <= _mem_rows && columns <= _mem_cols) {
        // Zero out the unused portions when reusing buffer
        // Zero partial columns (rows beyond new row count)
        for (UINT32 col = 0; col < columns && col < _mem_cols; ++col) {
            for (UINT32 row = rows; row < _mem_rows; ++row) {
                _buffer[col * _mem_rows + row] = 0.0;
            }
        }
        // Zero full columns beyond new column count
        if (columns < _mem_cols) {
            std::size_t start_idx = columns * _mem_rows;
            std::size_t count = (_mem_cols - columns) * _mem_rows;
            std::memset(_buffer + start_idx, 0, count * sizeof(double));
        }
        
        _rows = rows;
        _cols = columns;
        return;
    }

    //_method_ = "redim";

    
    // Save the current buffer and dimensions
    double* old_buffer = _buffer;
    UINT32 old_rows = _rows;
    UINT32 old_cols = _cols;
    
    // Allocate new buffer
    double* new_buffer;
    std::set_new_handler(out_of_memory_handler);
    buy(rows, columns, &new_buffer);
    
    // Copy old data to new buffer if there was any
    bool old_owns = _owns_buffer;
    if (old_buffer != nullptr && old_rows > 0 && old_cols > 0) {
        for (UINT32 col = 0; col < old_cols && col < columns; ++col) {
            for (UINT32 row = 0; row < old_rows && row < rows; ++row) {
                new_buffer[col * rows + row] = old_buffer[col * old_rows + row];
            }
        }
        if (old_owns)
            std::free(old_buffer);
    }
    
    _buffer = new_buffer;
    _owns_buffer = true;

    _rows = _mem_rows = rows;
    _cols = _mem_cols = columns;
}

void matrix_2d::redim_packed(const UINT32& n) {
    std::size_t ps = packed_size(n);

    if (_buffer != nullptr && _packed && n <= _mem_rows) {
        std::memset(_buffer, 0, ps * sizeof(double));
        _rows = _cols = _mem_rows = _mem_cols = n;
        return;
    }

    deallocate();
    _rows = _cols = _mem_rows = _mem_cols = n;
    _packed = true;
    _symmetric = true;
    _matrixType = mtx_lower;
    _buffer = static_cast<double*>(std::calloc(ps, sizeof(double)));
    if (!_buffer) {
        std::stringstream ss;
        ss << "Insufficient memory for a packed " << n << " x " << n << " matrix.";
        throw NetMemoryException(ss.str());
    }
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
    if (rows == _mem_rows && columns == _mem_cols &&
        oldmat.memRows() == _mem_rows && oldmat.memColumns() == _mem_cols) {
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

void matrix_2d::copyelements_generic(const UINT32& row_dest, const UINT32& column_dest, const matrix_2d& src,
                                     const UINT32& row_src, const UINT32& column_src, const UINT32& rows,
                                     const UINT32& columns) {
    assert(_buffer != nullptr && src._buffer != nullptr);
    assert(row_dest + rows <= _mem_rows && "copyelements_generic(): dest row overflow");
    assert(column_dest + columns <= _mem_cols && "copyelements_generic(): dest col overflow");
    assert(row_src + rows <= src._mem_rows && "copyelements_generic(): src row overflow");
    assert(column_src + columns <= src._mem_cols && "copyelements_generic(): src col overflow");
    // Fast path: packed source → dense dest (avoids per-element get/put overhead)
    if (src._packed && !_packed) {
        for (UINT32 c = 0; c < columns; ++c) {
            double* d = _buffer + static_cast<std::size_t>(column_dest + c) * _mem_rows + row_dest;
            UINT32 sc = column_src + c;
            for (UINT32 r = 0; r < rows; ++r) {
                UINT32 sr = row_src + r, sj = sc;
                if (sr < sj) std::swap(sr, sj);
                d[r] = src._buffer[packed_index(src._rows, sr, sj)];
            }
        }
        return;
    }
    // Fast path: dense source → packed dest
    if (_packed && !src._packed && !src._symmetric) {
        for (UINT32 c = 0; c < columns; ++c) {
            const double* s = src._buffer + static_cast<std::size_t>(column_src + c) * src._mem_rows + row_src;
            UINT32 dc = column_dest + c;
            for (UINT32 r = 0; r < rows; ++r) {
                UINT32 dr = row_dest + r;
                if (dr >= dc)
                    _buffer[packed_index(_rows, dr, dc)] = s[r];
            }
        }
        return;
    }
    // Fallback for remaining packed/symmetric combinations
    if (src._symmetric || src._packed || _packed) {
        for (UINT32 r = 0; r < rows; ++r)
            for (UINT32 c = 0; c < columns; ++c)
                put(row_dest + r, column_dest + c, src.get(row_src + r, column_src + c));
        return;
    }
    UINT32 cd(0), cs(0), colend_dest(column_dest + columns);
    for (cd = column_dest, cs = column_src; cd < colend_dest; ++cd, ++cs)
        memcpy(getelementref(row_dest, cd), src.getbuffer(row_src, cs),
               static_cast<std::size_t>(rows) * sizeof(double));
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
    for (k = k1; k < k2; k++)       //	for (k = k1; k <= k2; k++)
    {                               //	{
        if (fabs(get(k, k)) < eps)  //		if ( fabs( a(k, k) ) < eps)
        {
            for (it = 0; it < _rows; it++)  //			for (it = 1; it <= n; it++)
            {
                put(it, k, 0.);
                put(k, it, 0.);  //				a(it, k) = a(k, it) = 0.0;
            }
        } else {                                   //		else {
            d = 1.0 / get(k, k);                   //			d = 1.0 / a(k, k);
            put(k, k, d);                          //			a(k, k) = d;
            for (i = 0; i < _rows; i++)            //			for (i = 1; i <= n; i++)
                if (i != k)                        //				if (i != k)
                    *getelementref(i, k) *= -d;    //					a(i, k) *= (T) - d;
            for (j = 0; j < _rows; j++)            //			for (j = 1; j <= n; j++)
                if (j != k)                        //				if (j != k)
                    *getelementref(k, j) *= d;     //					a(k, j) *= (T) d;
            for (i = 0; i < _rows; i++) {          //			for (i = 1; i <= n; i++) {
                if (i != k) {                      //				if (i != k) {
                    for (j = 0; j < _rows; j++) {  //					for (j = 1; j <= n; j++) {
                        if (j != k)                //						if (j != k)
                            *getelementref(i, j) +=
                                get(i, k) * get(k, j) / d;  //							a(i, j) += a(i, k) *a(k, j) / d;
                    }  // end for j													//					} // end for j
                }  // end for i != k													//				} // end for i
                   // != k
            }  // end for i															//			} // end for i
        }  // end else																//		} // end else
    }  // end for k																	//	} // end for k
}

matrix_2d matrix_2d::sweepinverse() {
    if (_rows != _cols) throw std::runtime_error("sweepinverse: Matrix is not square.");

    sweep(0, _rows);
    return *this;
}

matrix_2d matrix_2d::cholesky_inverse(bool LOWER_IS_CLEARED /*=false*/, bool mark_symmetric /*=false*/) {
    if (_rows < 1) return *this;
    if (_rows != _cols) throw std::runtime_error("cholesky_inverse(): Matrix is not square.");

    assert(_buffer != nullptr && "cholesky_inverse(): null buffer");
    assert(_mem_rows >= _rows && "cholesky_inverse(): mem_rows < rows");
    assert(_mem_cols >= _cols && "cholesky_inverse(): mem_cols < cols");
    assert(!(mark_symmetric && LOWER_IS_CLEARED) &&
           "cholesky_inverse(): mark_symmetric with LOWER_IS_CLEARED not supported");

    if (_packed) {
        // Unpack to full format for fast blocked dpotrf/dpotri, then repack.
        // dpptrf/dpptri are element-by-element and orders of magnitude slower.
        lapack_int n = _rows;
        std::size_t full_size = static_cast<std::size_t>(n) * n;

        // Reuse thread-local workspace to avoid repeated allocation/deallocation.
        // dpotrf/dpotri with uplo='L' only read the lower triangle, so no
        // memset is needed — the unpack loop sets all lower-triangle elements.
        static thread_local std::vector<double> chol_workspace;
        if (chol_workspace.size() < full_size)
            chol_workspace.resize(full_size);
        double* full = chol_workspace.data();

        for (UINT32 j = 0; j < _rows; ++j)
            for (UINT32 i = j; i < _rows; ++i)
                full[static_cast<std::size_t>(j) * n + i] = _buffer[packed_index(_rows, i, j)];

        int blas_threads = blas_threads_for(n);
#if defined(USE_MKL) || defined(__MKL__)
        mkl_set_num_threads(blas_threads);
#elif !defined(__APPLE__)
        openblas_set_num_threads(blas_threads);
#endif

        char uplo = LOWER_TRIANGLE;
        lapack_int info, lda = n;
        LAPACK_FUNC(dpotrf)(&uplo, &n, full, &lda, &info);
        if (info != 0) throw MatrixInversionFailure("Matrix inversion failed, the matrix is singular.");
        LAPACK_FUNC(dpotri)(&uplo, &n, full, &lda, &info);
        if (info != 0) throw MatrixInversionFailure("Matrix inversion failed, the matrix is singular.");

        for (UINT32 j = 0; j < _rows; ++j)
            for (UINT32 i = j; i < _rows; ++i)
                _buffer[packed_index(_rows, i, j)] = full[static_cast<std::size_t>(j) * n + i];
        return *this;
    }

    char uplo(LOWER_TRIANGLE);
    if (LOWER_IS_CLEARED) uplo = UPPER_TRIANGLE;

    lapack_int info, n = _rows;
    lapack_int lda = _mem_rows;

    int blas_threads = blas_threads_for(n);
#if defined(USE_MKL) || defined(__MKL__)
    mkl_set_num_threads(blas_threads);
#elif !defined(__APPLE__)
    openblas_set_num_threads(blas_threads);
#endif

    // Perform Cholesky factorisation
    LAPACK_FUNC(dpotrf)(&uplo, &n, _buffer, &lda, &info);

    if (info != 0)
        throw MatrixInversionFailure("Matrix inversion failed, the matrix is singular.");

    // Perform Cholesky inverse
    LAPACK_FUNC(dpotri)(&uplo, &n, _buffer, &lda, &info);

    if (info != 0)
        throw MatrixInversionFailure("Matrix inversion failed, the matrix is singular.");

    if (mark_symmetric) {
        _symmetric = true;
    } else if (LOWER_IS_CLEARED) {
        filllower();
    } else {
        fillupper();
    }

    return *this;
}

matrix_2d matrix_2d::cholesky_factor(bool LOWER_IS_CLEARED /*=false*/) {
    if (_rows < 1) return *this;
    if (_rows != _cols) throw std::runtime_error("cholesky_factor(): Matrix is not square.");

    assert(_buffer != nullptr && "cholesky_factor(): null buffer");
    assert(_mem_rows >= _rows && "cholesky_factor(): mem_rows < rows");
    assert(_mem_cols >= _cols && "cholesky_factor(): mem_cols < cols");

    if (_packed) {
        // Unpack to full format for dpotrf (same approach as cholesky_inverse).
        // The factor L is stored back in packed format.
        lapack_int n = _rows;
        std::size_t full_size = static_cast<std::size_t>(n) * n;

        // Reuse thread-local workspace; no memset needed (dpotrf reads lower triangle only)
        static thread_local std::vector<double> chol_workspace;
        if (chol_workspace.size() < full_size)
            chol_workspace.resize(full_size);
        double* full = chol_workspace.data();

        for (UINT32 j = 0; j < _rows; ++j)
            for (UINT32 i = j; i < _rows; ++i)
                full[static_cast<std::size_t>(j) * n + i] = _buffer[packed_index(_rows, i, j)];

        int blas_threads = blas_threads_for(n);
#if defined(USE_MKL) || defined(__MKL__)
        mkl_set_num_threads(blas_threads);
#elif !defined(__APPLE__)
        openblas_set_num_threads(blas_threads);
#endif

        char uplo = LOWER_TRIANGLE;
        lapack_int info, lda = n;
        LAPACK_FUNC(dpotrf)(&uplo, &n, full, &lda, &info);
        if (info != 0) throw MatrixInversionFailure("Cholesky factorisation failed, the matrix is singular.");

        // Store the L factor back in packed format
        for (UINT32 j = 0; j < _rows; ++j)
            for (UINT32 i = j; i < _rows; ++i)
                _buffer[packed_index(_rows, i, j)] = full[static_cast<std::size_t>(j) * n + i];
        return *this;
    }

    char uplo(LOWER_TRIANGLE);
    if (LOWER_IS_CLEARED) uplo = UPPER_TRIANGLE;

    lapack_int info, n = _rows;
    lapack_int lda = _mem_rows;

    int blas_threads = blas_threads_for(n);
#if defined(USE_MKL) || defined(__MKL__)
    mkl_set_num_threads(blas_threads);
#elif !defined(__APPLE__)
    openblas_set_num_threads(blas_threads);
#endif

    LAPACK_FUNC(dpotrf)(&uplo, &n, _buffer, &lda, &info);

    if (info != 0)
        throw MatrixInversionFailure("Cholesky factorisation failed, the matrix is singular.");

    return *this;
}

void matrix_2d::cholesky_solve(matrix_2d& rhs, bool LOWER_IS_CLEARED /*=false*/) {
    assert(_rows == _cols && "cholesky_solve(): factor matrix is not square");
    assert(_rows == rhs._rows && "cholesky_solve(): dimension mismatch");
    assert(_buffer != nullptr && "cholesky_solve(): null factor buffer");
    assert(rhs._buffer != nullptr && "cholesky_solve(): null rhs buffer");

    if (_packed) {
        // Unpack the L factor to full format for dpotrs.
        lapack_int n = _rows;
        std::size_t full_size = static_cast<std::size_t>(n) * n;

        // Reuse thread-local workspace; no memset needed (dpotrs reads lower triangle only)
        static thread_local std::vector<double> chol_workspace;
        if (chol_workspace.size() < full_size)
            chol_workspace.resize(full_size);
        double* full = chol_workspace.data();

        for (UINT32 j = 0; j < _rows; ++j)
            for (UINT32 i = j; i < _rows; ++i)
                full[static_cast<std::size_t>(j) * n + i] = _buffer[packed_index(_rows, i, j)];

        char uplo = LOWER_TRIANGLE;
        lapack_int nrhs = rhs._cols;
        lapack_int lda = n;
        lapack_int ldb = rhs._mem_rows;
        lapack_int info;

        LAPACK_FUNC(dpotrs)(&uplo, &n, &nrhs, full, &lda, rhs._buffer, &ldb, &info);

        if (info != 0)
            throw MatrixInversionFailure("Cholesky solve failed.");
        return;
    }

    char uplo(LOWER_TRIANGLE);
    if (LOWER_IS_CLEARED) uplo = UPPER_TRIANGLE;

    lapack_int n = _rows;
    lapack_int nrhs = rhs._cols;
    lapack_int lda = _mem_rows;
    lapack_int ldb = rhs._mem_rows;
    lapack_int info;

    LAPACK_FUNC(dpotrs)(&uplo, &n, &nrhs, _buffer, &lda, rhs._buffer, &ldb, &info);

    if (info != 0)
        throw MatrixInversionFailure("Cholesky solve failed.");
}

double matrix_2d::dot(const matrix_2d& other) const {
    assert(_cols == 1 && other._cols == 1 && "dot(): both matrices must be column vectors");
    assert(_rows == other._rows && "dot(): dimension mismatch");
    assert(_buffer != nullptr && other._buffer != nullptr && "dot(): null buffer");

    double result = 0.0;
    for (UINT32 i = 0; i < _rows; ++i)
        result += get(i, 0) * other.get(i, 0);
    return result;
}

matrix_2d matrix_2d::scale(const double& scalar) {
    if (_packed) {
        std::size_t ps = packed_size(_rows);
        for (std::size_t k = 0; k < ps; ++k)
            _buffer[k] *= scalar;
        return *this;
    }
    UINT32 i, j;
    for (i = 0; i < _rows; ++i)
        for (j = 0; j < _cols; ++j) *getelementref(i, j) *= scalar;
    return *this;
}

void matrix_2d::scale_symmetric_diagonal(const double* diag) {
    assert(_symmetric && "scale_symmetric_diagonal(): matrix must be symmetric");
    if (_packed) {
        for (UINT32 j = 0; j < _rows; ++j)
            for (UINT32 i = j; i < _rows; ++i)
                _buffer[packed_index(_rows, i, j)] *= diag[i] * diag[j];
        return;
    }
    for (UINT32 j = 0; j < _rows; ++j)
        for (UINT32 i = j; i < _rows; ++i) {
            double s = diag[i] * diag[j];
            *getelementref(i, j) *= s;
        }
}

void matrix_2d::blockadd_generic(const UINT32& row_dest, const UINT32& col_dest, const matrix_2d& mat_src,
                                 const UINT32& row_src, const UINT32& col_src, const UINT32& rows, const UINT32& cols) {
    // Fast path: dense source → packed dest (avoids per-element get/elementadd overhead)
    if (_packed && !mat_src._packed && !mat_src._symmetric) {
        for (UINT32 c = 0; c < cols; ++c) {
            const double* s = mat_src._buffer + static_cast<std::size_t>(col_src + c) * mat_src._mem_rows + row_src;
            UINT32 dc = col_dest + c;
            for (UINT32 r = 0; r < rows; ++r) {
                UINT32 dr = row_dest + r;
                if (dr >= dc)
                    _buffer[packed_index(_rows, dr, dc)] += s[r];
            }
        }
        return;
    }
    // Fast path: both dense (non-packed, non-symmetric) — use direct buffer arithmetic
    if (!_packed && !mat_src._packed && !mat_src._symmetric) {
        for (UINT32 c = 0; c < cols; ++c) {
            double* d = _buffer + static_cast<std::size_t>(col_dest + c) * _mem_rows + row_dest;
            const double* s = mat_src._buffer + static_cast<std::size_t>(col_src + c) * mat_src._mem_rows + row_src;
            for (UINT32 r = 0; r < rows; ++r)
                d[r] += s[r];
        }
        return;
    }
    // Fallback: generic per-element path
    UINT32 i_dest, j_dest, i_src, j_src;
    UINT32 i_dest_end(row_dest + rows), j_dest_end(col_dest + cols);

    for (i_dest = row_dest, i_src = row_src; i_dest < i_dest_end; ++i_dest, ++i_src)
        for (j_dest = col_dest, j_src = col_src; j_dest < j_dest_end; ++j_dest, ++j_src)
            elementadd(i_dest, j_dest, mat_src.get(i_src, j_src));
}

void matrix_2d::blockTadd_generic(const UINT32& row_dest, const UINT32& col_dest, const matrix_2d& mat_src,
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

void matrix_2d::clearlower() {
    assert(_buffer != nullptr);
    assert(!_packed && "clearlower(): not valid for packed storage");
    // Sets lower triangle elements to zero
    UINT32 col, row;
    for (row = 1, col = 0; col < _mem_cols && row < _mem_rows; ++col, ++row)
        memset(getelementref(row, col), 0, (static_cast<std::size_t>(_mem_rows) - row) * sizeof(double));
}

void matrix_2d::clearupper() {
    assert(!_packed && "clearupper(): not valid for packed storage");
    // Sets upper triangle elements to zero
    // Column-major: for each column, zero the rows above the diagonal
    for (UINT32 col = 1; col < _cols; ++col)
        memset(_buffer + col * _mem_rows, 0, col * sizeof(double));
}

void matrix_2d::filllower() {
    assert(_buffer != nullptr);
    assert(!_packed && "filllower(): not valid for packed storage");
    assert(_rows == _cols && "filllower(): matrix must be square");
    // Copy upper triangle to lower: A(row,col) = A(col,row) for row > col
    // Column-major: A(r,c) = _buffer[c * _mem_rows + r]
    // Outer loop over destination columns so writes are sequential in memory.
    const std::size_t mr = _mem_rows;
    constexpr UINT32 BLK = 64;
    for (UINT32 cb = 0; cb < _cols; cb += BLK) {
        UINT32 ce = std::min(cb + BLK, _cols);
        for (UINT32 rb = ce; rb < _rows; rb += BLK) {
            UINT32 re = std::min(rb + BLK, _rows);
            for (UINT32 col = cb; col < ce; col++) {
                double* dst = _buffer + col * mr;
                for (UINT32 row = rb; row < re; row++)
                    dst[row] = _buffer[row * mr + col];
            }
        }
    }
    // Diagonal blocks: row > col within the same block
    for (UINT32 cb = 0; cb < _cols; cb += BLK) {
        UINT32 ce = std::min(cb + BLK, _cols);
        for (UINT32 col = cb; col < ce; col++) {
            double* dst = _buffer + col * mr;
            for (UINT32 row = col + 1; row < ce; row++)
                dst[row] = _buffer[row * mr + col];
        }
    }
}

void matrix_2d::fillupper() {
    assert(_buffer != nullptr);
    assert(!_packed && "fillupper(): not valid for packed storage");
    assert(_rows == _cols && "fillupper(): matrix must be square");
    // Copy lower triangle to upper: A(r,c) = A(c,r) for r < c
    // Column-major: A(r,c) = _buffer[c * _mem_rows + r]
    //
    // Outer loop over destination columns so writes are sequential in memory.
    // Reads are strided (one element per source column) but write-combining
    // and store buffers make sequential writes much faster than sequential reads
    // for large matrices.
    const std::size_t mr = _mem_rows;
    constexpr UINT32 BLK = 64;
    for (UINT32 cb = 0; cb < _cols; cb += BLK) {
        UINT32 ce = std::min(cb + BLK, _cols);
        for (UINT32 rb = 0; rb < ce; rb += BLK) {
            UINT32 re = std::min(rb + BLK, ce);
            for (UINT32 col = cb; col < ce; col++) {
                double* dst = _buffer + col * mr;
                UINT32 r0 = rb;
                UINT32 r1 = std::min(re, col);
                for (UINT32 row = r0; row < r1; row++)
                    dst[row] = _buffer[row * mr + col];
            }
        }
    }
}

void matrix_2d::zero() {
    assert(_buffer != nullptr && "zero(): null buffer");
    if (_packed)
        memset(_buffer, 0, packed_size(_mem_rows) * sizeof(double));
    else
        memset(_buffer, 0, buffersize());
}

// zero()
void matrix_2d::zero(const UINT32& row_begin, const UINT32& col_begin, const UINT32& rows, const UINT32& columns) {
    assert(_buffer != nullptr && "zero(sub): null buffer");
    assert(row_begin + rows <= _mem_rows && "zero(sub): row overflow");
    assert(col_begin + columns <= _mem_cols && "zero(sub): col overflow");
    UINT32 col(0), col_end(col_begin + columns);
    for (col = col_begin; col < col_end; ++col) memset(getelementref(row_begin, col), 0, rows * sizeof(double));
}


matrix_2d& matrix_2d::operator=(const matrix_2d& rhs) {
    if (this == &rhs) return *this;

    if (rhs._packed) {
        std::size_t ps = packed_size(rhs._rows);
        if (_packed && _mem_rows >= rhs._rows) {
            _rows = _cols = rhs._rows;
            memcpy(_buffer, rhs._buffer, ps * sizeof(double));
        } else {
            deallocate();
            _rows = _cols = _mem_rows = _mem_cols = rhs._rows;
            _buffer = static_cast<double*>(std::malloc(ps * sizeof(double)));
            if (!_buffer) throw NetMemoryException("Insufficient memory for packed matrix assignment.");
            memcpy(_buffer, rhs._buffer, ps * sizeof(double));
        }
        _packed = true;
        _symmetric = true;
        _matrixType = rhs._matrixType;
        _maxvalCol = rhs.maxvalueCol();
        _maxvalRow = rhs.maxvalueRow();
        return *this;
    }

    // rhs is not packed
    if (_packed) {
        deallocate();
        _packed = false;
        _mem_rows = _mem_cols = 0;
    }

    if (_mem_rows >= rhs.rows() && _mem_cols >= rhs.columns()) {
        _rows = rhs.rows();
        _cols = rhs.columns();
        copybuffer(_rows, _cols, rhs);

        _maxvalCol = rhs.maxvalueCol();
        _maxvalRow = rhs.maxvalueRow();
        _symmetric = rhs._symmetric;

        return *this;
    }

    deallocate();
    _mem_rows = rhs.memRows();
    _mem_cols = rhs.memColumns();
    _rows = rhs.rows();
    _cols = rhs.columns();
    allocate(_mem_rows, _mem_cols);
    copybuffer(_rows, _cols, rhs);

    _maxvalCol = rhs.maxvalueCol();
    _maxvalRow = rhs.maxvalueRow();
    _symmetric = rhs._symmetric;

    return *this;
}

// Multiplication operator
matrix_2d matrix_2d::operator*(const double& rhs) const {
    // Answer
    matrix_2d m(_rows, _cols);

    UINT32 row, column;
    for (row = 0; row < _rows; row++)
        for (column = 0; column < _cols; ++column) m.put(row, column, get(row, column) * rhs);
    return m;
}

matrix_2d matrix_2d::add(const matrix_2d& rhs) {
    if (_rows != rhs.rows() || _cols != rhs.columns())
        throw std::runtime_error("add: Result matrix dimensions are incompatible.");

    UINT32 row, column;
    for (row = 0; row < _rows; row++) {
        for (column = 0; column < _cols; ++column) { *getelementref(row, column) += rhs.get(row, column); }
    }
    return *this;
}

// multiplies this matrix by rhs and stores the result in a new matrix
// Uses Intel MKL dgemm
matrix_2d matrix_2d::multiply(const char* lhs_trans, const matrix_2d& rhs, const char* rhs_trans) {
    assert(!_symmetric && "multiply(dgemm): LHS is symmetric — use multiply_sym instead");
    assert(!rhs._symmetric && "multiply(dgemm): RHS is symmetric — upper triangle is unpopulated");
    assert(_buffer != nullptr && rhs._buffer != nullptr);

    matrix_2d m(_rows, rhs.columns());

    const double one = 1.0;
    const double zero = 0.0;

    lapack_int lhs_rows(rows()), rhs_cols(rhs.columns());
    lapack_int lhs_cols(columns()), rhs_rows(rhs.rows());
    lapack_int new_mem_rows(memRows());
    lapack_int lhs_mem_rows(memRows());
    lapack_int rhs_mem_rows(memRows());

    if (strcmp(lhs_trans, "T") == 0) {
        lhs_rows = columns();  // transpose
        lhs_cols = rows();     // transpose
    }

    if (strcmp(rhs_trans, "T") == 0) {
        rhs_rows = rhs.columns();  // transpose
        rhs_cols = rhs.rows();     // transpose
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
    assert(!lhs._symmetric && "multiply(dgemm): LHS is symmetric — use multiply_sym instead");
    assert(!rhs._symmetric && "multiply(dgemm): RHS is symmetric — upper triangle is unpopulated");
    assert(lhs._buffer != nullptr && rhs._buffer != nullptr && _buffer != nullptr);

    const double one = 1.0;
    const double zero = 0.0;

    lapack_int lhs_rows(lhs.rows()), rhs_cols(rhs.columns());
    lapack_int lhs_cols(lhs.columns()), rhs_rows(rhs.rows());
    lapack_int new_mem_rows(memRows());
    lapack_int lhs_mem_rows(lhs.memRows());
    lapack_int rhs_mem_rows(rhs.memRows());

    if (strncmp(lhs_trans, "T", 1) == 0) {
        lhs_rows = lhs.columns();  // transpose
        lhs_cols = lhs.rows();     // transpose
    }

    if (strncmp(rhs_trans, "T", 1) == 0) {
        rhs_rows = rhs.columns();  // transpose
        rhs_cols = rhs.rows();     // transpose
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
}  // Multiply()

// C (this) = sym_lhs * rhs, using dsymm or dspmv (packed)
matrix_2d matrix_2d::multiply_sym(const matrix_2d& sym_lhs, const matrix_2d& rhs) {
    lapack_int m = sym_lhs.rows();
    lapack_int n = rhs.columns();

    assert(sym_lhs._symmetric && "multiply_sym(): LHS must be marked symmetric");
    assert(sym_lhs._buffer != nullptr && "multiply_sym(): LHS null buffer");
    assert(rhs._buffer != nullptr && "multiply_sym(): RHS null buffer");
    assert(_buffer != nullptr && "multiply_sym(): result null buffer");
    assert(_buffer != sym_lhs._buffer && _buffer != rhs._buffer &&
           "multiply_sym(): result must not alias inputs");

    if (sym_lhs.columns() != sym_lhs.rows())
        throw std::runtime_error("multiply_sym(): LHS matrix is not square.");
    if (static_cast<lapack_int>(rhs.rows()) != m)
        throw std::runtime_error("multiply_sym(): Matrix dimensions are incompatible.");
    if (_rows != static_cast<UINT32>(m) || _cols != static_cast<UINT32>(n))
        throw std::runtime_error("multiply_sym(): Result matrix dimensions are incompatible.");

    if (sym_lhs._packed) {
        for (lapack_int j = 0; j < n; ++j) {
            BLAS_FUNC(dspmv)(CblasColMajor, CblasLower,
                             m, 1.0, sym_lhs._buffer,
                             rhs._buffer + j * rhs._mem_rows, 1,
                             0.0, _buffer + j * _mem_rows, 1);
        }
        return *this;
    }

    assert(sym_lhs.memRows() >= sym_lhs.rows() && "multiply_sym(): LHS LDA < M");
    assert(rhs.memRows() >= rhs.rows() && "multiply_sym(): RHS LDA < M");
    assert(_mem_rows >= _rows && "multiply_sym(): result LDA < M");

    BLAS_FUNC(dsymm)(CblasColMajor, CblasLeft, CblasLower,
                     m, n, 1.0,
                     sym_lhs.getbuffer(), sym_lhs.memRows(),
                     rhs.getbuffer(), rhs.memRows(),
                     0.0, _buffer, _mem_rows);

    return *this;
}

// Transpose()
matrix_2d matrix_2d::transpose(const matrix_2d& matA) {
    if ((matA.columns() != _rows) || (matA.rows() != _cols))
        throw std::runtime_error("transpose: Matrix dimensions are incompatible.");

    UINT32 column, row;
    for (row = 0; row < _rows; row++)
        for (column = 0; column < _cols; column++) *getelementref(row, column) = matA.get(column, row);
    return *this;
}  // Transpose()

// Transpose()
matrix_2d matrix_2d::transpose() {
    matrix_2d m(_cols, _rows);
    UINT32 column, row;
    for (row = 0; row < _rows; row++)
        for (column = 0; column < _cols; column++) m.put(column, row, get(row, column));
    return m;
}  // Transpose()

double matrix_2d::compute_maximum_value() {
    _maxvalCol = _maxvalRow = 0;
    if (_packed) {
        for (UINT32 j = 0; j < _rows; ++j) {
            for (UINT32 i = j; i < _rows; ++i) {
                if (fabs(get(i, j)) > fabs(get(_maxvalRow, _maxvalCol))) {
                    _maxvalRow = i;
                    _maxvalCol = j;
                }
            }
        }
        return get(_maxvalRow, _maxvalCol);
    }
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
        for (j = 0; j < _cols; ++j) TRACE(format.c_str(), get(i, j));
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
        for (j = col_begin; j < col_end; ++j) TRACE(format.c_str(), get(i, j));
        TRACE("\n");
    }
    TRACE("\n");
}
#endif

}  // namespace math
}  // namespace dynadjust
