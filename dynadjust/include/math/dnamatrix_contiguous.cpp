#include <include/math/dnamatrix_contiguous.hpp>

namespace dynadjust {
namespace math {

std::ostream& operator<<(std::ostream& os, const matrix_2d& rhs) {
  if (os.iword(0) == binary) {
    os.write(reinterpret_cast<const char*>(&rhs._matrixType), sizeof(index_t));
    os.write(reinterpret_cast<const char*>(&rhs._rows), sizeof(index_t));
    os.write(reinterpret_cast<const char*>(&rhs._cols), sizeof(index_t));
    os.write(reinterpret_cast<const char*>(&rhs._mem_rows), sizeof(index_t));
    os.write(reinterpret_cast<const char*>(&rhs._mem_cols), sizeof(index_t));
    index_t c, r;
    switch (rhs._matrixType) {
    case mtx_lower:
      if (rhs._mem_rows != rhs._mem_cols) {
        throw boost::enable_current_exception(
            std::runtime_error("matrix_2d operator<< (): Matrix not square."));
      }
      for (c = 0; c < rhs._mem_cols; ++c) {
        os.write(reinterpret_cast<const char*>(rhs.getelementref(c, c)),
                 (rhs._mem_rows - c) * sizeof(double));
      }
      break;
    case mtx_sparse: break;
    case mtx_full:
    default:
      for (r = 0; r < rhs._mem_rows; ++r) {
        for (c = 0; c < rhs._mem_cols; ++c) {
          os.write(reinterpret_cast<const char*>(rhs.getelementref(r, c)),
                   sizeof(double));
        }
      }
      break;
    }
    os.write(reinterpret_cast<const char*>(&rhs._maxvalRow), sizeof(index_t));
    os.write(reinterpret_cast<const char*>(&rhs._maxvalCol), sizeof(index_t));
  } else {
    os << rhs._matrixType << " " << rhs._rows << " " << rhs._cols << " "
       << rhs._mem_rows << " " << rhs._mem_cols << std::endl;
    for (index_t r = 0; r < rhs._mem_rows; ++r) {
      for (index_t c = 0; c < rhs._mem_cols; ++c) {
        os << std::scientific << std::setprecision(16) << rhs.get(r, c) << " ";
      }
      os << std::endl;
    }
    os << rhs._maxvalRow << " " << rhs._maxvalCol << std::endl;
    os << std::endl;
  }
  return os;
}

index_t __row__;
index_t __col__;

void out_of_memory_handler() {
  std::stringstream ss("");
  index_t mem(std::max(__row__, __col__));
  mem *= std::min(__row__, __col__) * sizeof(index_t);
  ss << "Insufficient memory available to create a " << __row__ << " x "
     << __col__ << " matrix (" << std::fixed << std::setprecision(2);
  if (mem < MEGABYTE_SIZE)
    ss << (mem / KILOBYTE_SIZE) << " KB).";
  else if (mem < GIGABYTE_SIZE)
    ss << (mem / MEGABYTE_SIZE) << " MB).";
  else
    ss << (mem / GIGABYTE_SIZE) << " GB).";
  throw boost::enable_current_exception(NetMemoryException(ss.str()));
}

matrix_2d::matrix_2d()
    : _mem_cols(0), _mem_rows(0), _cols(0), _rows(0), _buffer(0), _maxvalCol(0),
      _maxvalRow(0), _matrixType(mtx_full) {
  std::set_new_handler(out_of_memory_handler);
}

matrix_2d::matrix_2d(const index_t& rows, const index_t& columns)
    : _mem_cols(columns), _mem_rows(rows), _cols(columns), _rows(rows),
      _buffer(0), _maxvalCol(0), _maxvalRow(0), _matrixType(mtx_full) {
  std::set_new_handler(out_of_memory_handler);
  allocate(_rows, _cols);
}

matrix_2d::matrix_2d(const index_t& rows, const index_t& columns,
                     const double data[], const index_t& data_size,
                     const index_t& matrix_type)
    : _mem_cols(columns), _mem_rows(rows), _cols(columns), _rows(rows),
      _buffer(0), _maxvalCol(0), _maxvalRow(0), _matrixType(matrix_type) {
  std::set_new_handler(out_of_memory_handler);
  std::stringstream ss;
  index_t upperMatrixElements(sumOfConsecutiveIntegers(rows));
  index_t j;
  const double* dataptr = &data[0];
  switch (matrix_type) {
  case mtx_lower:
    if (upperMatrixElements != data_size) {
      ss << "Data size must match upper matrix element count for " << rows
         << " x " << columns << ".";
      throw boost::enable_current_exception(std::runtime_error(ss.str()));
    }
    allocate(_rows, _cols);
    for (j = 0; j < columns; ++j) {
      memcpy(getelementref(j, j), dataptr, (rows - j) * sizeof(double));
      dataptr += (rows - j);
    }
    fillupper();
    break;

  case mtx_sparse:
    ss << "matrix_2d(): Cannot initialise sparse matrix from double array.";
    throw boost::enable_current_exception(std::runtime_error(ss.str()));
    break;

  case mtx_full:
  default:
    if (data_size != rows * columns) {
      ss << "Data size must match matrix dimensions (" << rows << " x "
         << columns << ").";
      throw boost::enable_current_exception(std::runtime_error(ss.str()));
    }
    allocate(_rows, _cols);
    memcpy(_buffer, data, data_size * sizeof(double));
    break;
  }
}

matrix_2d::matrix_2d(const matrix_2d& newmat)
    : _mem_cols(newmat.memColumns()), _mem_rows(newmat.memRows()),
      _cols(newmat.columns()), _rows(newmat.rows()), _buffer(0),
      _maxvalCol(newmat.maxvalueCol()), _maxvalRow(newmat.maxvalueRow()),
      _matrixType(newmat.matrixType()) {
  std::set_new_handler(out_of_memory_handler);
  allocate(_mem_rows, _mem_cols);
  const double* ptr = newmat.getbuffer();
  memcpy(_buffer, ptr, newmat.buffersize());
}

matrix_2d::~matrix_2d() { deallocate(); }

std::size_t matrix_2d::get_size() {
  size_t size = (7 * sizeof(index_t));
  switch (_matrixType) {
  case mtx_lower:
    size += sumOfConsecutiveIntegers(_mem_rows) * sizeof(double);
    break;
  case mtx_sparse: break;
  case mtx_full:
  default: size += buffersize(); break;
  }
  return size;
}

void matrix_2d::ReadMappedFileRegion(void* addr) {
  index_t* data_U = reinterpret_cast<index_t*>(addr);
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
    break;
  }

  allocate(_mem_rows, _mem_cols);

  double* data_d;
  int* data_i;
  index_t c, r, i;
  int ci;

  switch (_matrixType) {
  case mtx_sparse:
    data_i = reinterpret_cast<int*>(data_U);
    for (r = 0; r < _rows; ++r) {
      for (i = 0; i < 3; ++i) {
        ci = *data_i++;
        data_d = reinterpret_cast<double*>(data_i);
        if (ci < 0) {
          data_d += 3;
          data_i = reinterpret_cast<int*>(data_d);
          continue;
        }
        memcpy(getelementref(r, ci), data_d, sizeof(double));
        data_d++;
        memcpy(getelementref(r, ci + 1), data_d, sizeof(double));
        data_d++;
        memcpy(getelementref(r, ci + 2), data_d, sizeof(double));
        data_d++;
        data_i = reinterpret_cast<int*>(data_d);
      }
    }
    return;
    break;

  case mtx_lower:
    data_d = reinterpret_cast<double*>(data_U);
    for (c = 0; c < _mem_cols; ++c) {
      memcpy(getelementref(c, c), data_d, (_mem_rows - c) * sizeof(double));
      data_d += (_mem_rows - c);
    }
    fillupper();
    break;

  case mtx_full:
  default:
    data_d = reinterpret_cast<double*>(data_U);
    memcpy(_buffer, data_d, buffersize());
    data_d += _mem_rows * _mem_cols;
    break;
  }
  data_U = reinterpret_cast<index_t*>(data_d);
  _maxvalRow = *data_U++;
  _maxvalCol = *data_U;
}

void matrix_2d::WriteMappedFileRegion(void* addr) {
  index_t* data_U = reinterpret_cast<index_t*>(addr);
  *data_U++ = _matrixType;
  *data_U++ = _rows;
  *data_U++ = _cols;

  switch (_matrixType) {
  case mtx_sparse: break;
  case mtx_lower:
  case mtx_full:
  default:
    *data_U++ = _mem_rows;
    *data_U++ = _mem_cols;
    break;
  }

  double* data_d;
  int* data_i;
  index_t c, r, i;
  int ci;

  switch (_matrixType) {
  case mtx_sparse:
    data_i = reinterpret_cast<int*>(data_U);
    for (r = 0; r < _rows; ++r) {
      for (i = 0; i < 3; ++i) {
        ci = *data_i++;
        data_d = reinterpret_cast<double*>(data_i);
        if (ci < 0) {
          data_d += 3;
          data_i = reinterpret_cast<int*>(data_d);
          continue;
        }
        memcpy(data_d, getelementref(r, ci), sizeof(double));
        data_d++;
        memcpy(data_d, getelementref(r, ci + 1), sizeof(double));
        data_d++;
        memcpy(data_d, getelementref(r, ci + 2), sizeof(double));
        data_d++;
        data_i = reinterpret_cast<int*>(data_d);
      }
    }
    return;
    break;

  case mtx_lower:
    data_d = reinterpret_cast<double*>(data_U);
    for (c = 0; c < _mem_cols; ++c) {
      memcpy(data_d, getbuffer(c, c), (_mem_rows - c) * sizeof(double));
      data_d += (_mem_rows - c);
    }
    break;

  case mtx_full:
  default:
    data_d = reinterpret_cast<double*>(data_U);
    memcpy(data_d, _buffer, _mem_rows * _mem_cols * sizeof(double));
    data_d += _mem_rows * _mem_cols;
    break;
  }

  data_U = reinterpret_cast<index_t*>(data_d);
  *data_U++ = _maxvalRow;
  *data_U = _maxvalCol;
}

void matrix_2d::allocate() { allocate(_mem_rows, _mem_cols); }

void matrix_2d::allocate(const index_t& rows, const index_t& columns) {
  deallocate();
  buy(rows, columns, &_buffer);
}

void matrix_2d::buy(const index_t& rows, const index_t& columns,
                    double** mem_space) {
  __row__ = rows;
  __col__ = columns;
  (*mem_space) = new double[rows * columns];
  if ((*mem_space) == NULL) {
    std::stringstream ss;
    ss << "Insufficient memory for a " << rows << " x " << columns
       << " matrix.";
    throw boost::enable_current_exception(NetMemoryException(ss.str()));
  }
  memset(*mem_space, 0, byteSize<double>(rows * columns));
}

void matrix_2d::deallocate() {
  if (_buffer != NULL) {
    delete[] _buffer;
    _buffer = 0;
  }
}

matrix_2d
matrix_2d::submatrix(const index_t& row_begin, const index_t& col_begin,
                     const index_t& rows, const index_t& columns) const {
  matrix_2d b(rows, columns);
  submatrix(row_begin, col_begin, &b, rows, columns);
  return b;
}

void matrix_2d::submatrix(const index_t& row_begin, const index_t& col_begin,
                          matrix_2d* dest, const index_t& subrows,
                          const index_t& subcolumns) const {
  if (row_begin >= _rows) {
    std::stringstream ss;
    ss << row_begin << ", " << col_begin << " lies outside matrix range ("
       << _rows << ", " << _cols << ").";
    throw boost::enable_current_exception(std::runtime_error(ss.str()));
  }
  if (col_begin >= _cols) {
    std::stringstream ss;
    ss << row_begin << ", " << col_begin << " lies outside matrix range ("
       << _rows << ", " << _cols << ").";
    throw boost::enable_current_exception(std::runtime_error(ss.str()));
  }
  if (subrows > dest->rows()) {
    std::stringstream ss;
    ss << subrows << ", " << subcolumns << " exceeds destination matrix ("
       << dest->rows() << ", " << dest->columns() << ").";
    throw boost::enable_current_exception(std::runtime_error(ss.str()));
  }
  if (subcolumns > dest->columns()) {
    std::stringstream ss;
    ss << subrows << ", " << subcolumns << " exceeds destination matrix ("
       << dest->rows() << ", " << dest->columns() << ").";
    throw boost::enable_current_exception(std::runtime_error(ss.str()));
  }
  if (row_begin + subrows > _rows) {
    std::stringstream ss;
    ss << row_begin + subrows << ", " << col_begin + subcolumns
       << " lies outside matrix range (" << _rows << ", " << _cols << ").";
    throw boost::enable_current_exception(std::runtime_error(ss.str()));
  }
  if (col_begin + subcolumns > _cols) {
    std::stringstream ss;
    ss << row_begin + subrows << ", " << col_begin + subcolumns
       << " lies outside matrix range (" << _rows << ", " << _cols << ").";
    throw boost::enable_current_exception(std::runtime_error(ss.str()));
  }

  index_t i, j, m(0), n(0), row_end(row_begin + subrows),
      col_end(col_begin + subcolumns);
  for (i = row_begin; i < row_end; ++i) {
    for (j = col_begin; j < col_end; ++j) {
      dest->put(m, n, get(i, j));
      n++;
    }
    n = 0;
    m++;
  }
}

void matrix_2d::redim(const index_t& rows, const index_t& columns) {
  if (rows <= _mem_rows && columns <= _mem_cols) {
    _rows = rows;
    _cols = columns;
    return;
  }
  std::set_new_handler(out_of_memory_handler);
  deallocate();
  allocate(rows, columns);
  _rows = _mem_rows = rows;
  _cols = _mem_cols = columns;
}

void matrix_2d::shrink(const index_t& rows, const index_t& columns) {
  if (rows > _rows || columns > _cols) {
    std::stringstream ss;
    ss << " " << std::endl;
    if (rows >= _rows) {
      ss << "    Cannot shrink by " << rows << " rows on matrix of " << _rows
         << " rows. " << std::endl;
    }
    if (columns >= _cols) {
      ss << "    Cannot shrink by " << columns << " columns on matrix of "
         << _cols << " columns.";
    }
    throw boost::enable_current_exception(std::runtime_error(ss.str()));
  }
  _rows -= rows;
  _cols -= columns;
}

void matrix_2d::grow(const index_t& rows, const index_t& columns) {
  if ((rows + _rows) > _mem_rows || (columns + _cols) > _mem_cols) {
    std::stringstream ss;
    ss << " " << std::endl;
    if (rows >= _rows) {
      ss << "    Cannot grow matrix by " << rows
         << " rows: exceeds row memory limit (" << _mem_rows << ").";
    }
    if (columns >= _cols) {
      ss << "    Cannot grow matrix by " << columns
         << " columns: exceeds column memory limit (" << _mem_cols << ").";
    }
    throw boost::enable_current_exception(std::runtime_error(ss.str()));
  }
  _rows += rows;
  _cols += columns;
}

void matrix_2d::setsize(const index_t& rows, const index_t& columns) {
  deallocate();
  _rows = _mem_rows = rows;
  _cols = _mem_cols = columns;
}

void matrix_2d::replace(const index_t& rowstart, const index_t& columnstart,
                        const matrix_2d& newmat) {
  copybuffer(rowstart, columnstart, newmat.rows(), newmat.columns(), newmat);
}

void matrix_2d::replace(const index_t& rowstart, const index_t& columnstart,
                        const index_t& rows, const index_t& columns,
                        const matrix_2d& newmat) {
  copybuffer(rowstart, columnstart, rows, columns, newmat);
}

void matrix_2d::copybuffer(const index_t& rows, const index_t& columns,
                           const matrix_2d& oldmat) {
  if (rows == _mem_rows && columns == _mem_cols) {
    memcpy(_buffer, oldmat.getbuffer(), buffersize());
    return;
  }
#if defined(DNAMATRIX_ROW_WISE)
  for (index_t row = 0; row < rows; row++) {
    memcpy(getelementref(row, 0), oldmat.getbuffer(row, 0),
           columns * sizeof(double));
  }
#else
  for (index_t column = 0; column < columns; column++) {
    memcpy(getelementref(0, column), oldmat.getbuffer(0, column),
           rows * sizeof(double));
  }
#endif
}

void matrix_2d::copybuffer(const index_t& rowstart, const index_t& columnstart,
                           const index_t& rows, const index_t& columns,
                           const matrix_2d& mat) {
  index_t rowend(rowstart + rows), columnend(columnstart + columns);
  if (rowend > _rows || columnend > _cols) {
    std::stringstream ss;
    ss << " " << std::endl;
    if (rowend >= _rows) {
      ss << "    Row index " << rowend << " exceeds matrix row count (" << _rows
         << "). " << std::endl;
    }
    if (columnend >= _cols) {
      ss << "    Column index " << columnend << " exceeds matrix column count ("
         << _cols << ").";
    }
    throw boost::enable_current_exception(std::runtime_error(ss.str()));
  }
#if defined(DNAMATRIX_ROW_WISE)
  index_t row(0), r(0);
  for (row = rowstart; row < rowend; ++row, ++r) {
    memcpy(getelementref(row, columnstart), mat.getbuffer(r, 0),
           columns * sizeof(double));
  }
#else
  index_t column(0), c(0);
  for (column = columnstart; column < columnend; ++column, ++c) {
    memcpy(getelementref(rowstart, column), mat.getbuffer(0, c),
           rows * sizeof(double));
  }
#endif
}

void matrix_2d::copyelements(const index_t& row_dest, const index_t& column_dest,
                             const matrix_2d& src, const index_t& row_src,
                             const index_t& column_src, const index_t& rows,
                             const index_t& columns) {
#if defined(DNAMATRIX_ROW_WISE)
  index_t rd(0), rs(0), rowend_dest(row_dest + rows);
  for (rd = row_dest, rs = row_src; rd < rowend_dest; ++rd, ++rs) {
    memcpy(getelementref(rd, column_dest), src.getbuffer(rs, column_src),
           columns * sizeof(double));
  }
#else
  index_t cd(0), cs(0), colend_dest(column_dest + columns);
  for (cd = column_dest, cs = column_src; cd < colend_dest; ++cd, ++cs) {
    memcpy(getelementref(row_dest, cd), src.getbuffer(row_src, cs),
           rows * sizeof(double));
  }
#endif
}

void matrix_2d::copyelements(const index_t& row_dest, const index_t& column_dest,
                             const matrix_2d* src, const index_t& row_src,
                             const index_t& column_src, const index_t& rows,
                             const index_t& columns) {
  copyelements(row_dest, column_dest, *src, row_src, column_src, rows, columns);
}

void matrix_2d::sweep(index_t k1, index_t k2) {
  double eps(1.0e-8), d;
  if (k2 < k1) {
    index_t k = k1;
    k1 = k2;
    k2 = k;
  }
  for (index_t k = k1; k < k2; k++) {
    if (fabs(get(k, k)) < eps) {
      for (index_t it = 0; it < _rows; it++) {
        put(it, k, 0.);
        put(k, it, 0.);
      }
    } else {
      d = 1.0 / get(k, k);
      put(k, k, d);
      for (index_t i = 0; i < _rows; i++) {
        if (i != k) {
          *getelementref(i, k) *= -d;
        }
      }
      for (index_t j = 0; j < _rows; j++) {
        if (j != k) {
          *getelementref(k, j) *= d;
        }
      }
      for (index_t i = 0; i < _rows; i++) {
        if (i != k) {
          for (index_t j = 0; j < _rows; j++) {
            if (j != k) {
              *getelementref(i, j) += get(i, k) * get(k, j) / d;
            }
          }
        }
      }
    }
  }
}

matrix_2d matrix_2d::sweepinverse() {
  if (_rows != _cols) {
    throw boost::enable_current_exception(
        std::runtime_error("sweepinverse(): Matrix is not square."));
  }
  sweep(0, _rows);
  return *this;
}

matrix_2d matrix_2d::cholesky_inverse(bool LOWER_IS_CLEARED) {
  if (_rows < 1)
    return *this;
  if (_rows != _cols) {
    throw boost::enable_current_exception(
        std::runtime_error("cholesky_inverse(): Matrix is not square."));
  }
  char uplo = LOWER_IS_CLEARED ? 'U' : 'L';

  lapack_int info(0);
  lapack_int n = static_cast<lapack_int>(_rows);

  info = LAPACKE_dpotrf(
#if defined(DNAMATRIX_ROW_WISE)
      LAPACK_ROW_MAJOR,
#else
      LAPACK_COL_MAJOR,
#endif
      uplo, n, _buffer, n);

  if (info != 0) {
    throw boost::enable_current_exception(std::runtime_error(
        "cholesky_inverse(): Cholesky factorisation failed."));
  }

  info = LAPACKE_dpotri(
#if defined(DNAMATRIX_ROW_WISE)
      LAPACK_ROW_MAJOR,
#else
      LAPACK_COL_MAJOR,
#endif
      uplo, n, _buffer, n);

  if (info != 0) {
    throw boost::enable_current_exception(std::runtime_error(
        "cholesky_inverse(): Cholesky inversion failed."));
  }

  if (LOWER_IS_CLEARED)
    filllower();
  else
    fillupper();

  return *this;
}

matrix_2d matrix_2d::scale(const double& scalar) {
  for (index_t i = 0; i < _rows; ++i) {
    for (index_t j = 0; j < _cols; ++j) {
      *getelementref(i, j) *= scalar;
    }
  }
  return *this;
}

void matrix_2d::blockadd(const index_t& row_dest, const index_t& col_dest,
                         const matrix_2d& mat_src, const index_t& row_src,
                         const index_t& col_src, const index_t& rows,
                         const index_t& cols) {
  for (index_t i_dest = row_dest, i_src = row_src; i_dest < row_dest + rows;
       ++i_dest, ++i_src) {
    for (index_t j_dest = col_dest, j_src = col_src; j_dest < col_dest + cols;
         ++j_dest, ++j_src) {
      elementadd(i_dest, j_dest, mat_src.get(i_src, j_src));
    }
  }
}

void matrix_2d::blockTadd(const index_t& row_dest, const index_t& col_dest,
                          const matrix_2d& mat_src, const index_t& row_src,
                          const index_t& col_src, const index_t& rows,
                          const index_t& cols) {
  for (index_t i_dest = row_dest, i_src = row_src; i_dest < row_dest + rows;
       ++i_dest, ++i_src) {
    for (index_t j_dest = col_dest, j_src = col_src; j_dest < col_dest + cols;
         ++j_dest, ++j_src) {
      elementadd(i_dest, j_dest, mat_src.get(j_src, i_src));
    }
  }
}

void matrix_2d::blocksubtract(const index_t& row_dest, const index_t& col_dest,
                              const matrix_2d& mat_src, const index_t& row_src,
                              const index_t& col_src, const index_t& rows,
                              const index_t& cols) {
  for (index_t i_dest = row_dest, i_src = row_src; i_dest < row_dest + rows;
       ++i_dest, ++i_src) {
    for (index_t j_dest = col_dest, j_src = col_src; j_dest < col_dest + cols;
         ++j_dest, ++j_src) {
      elementsubtract(i_dest, j_dest, mat_src.get(i_src, j_src));
    }
  }
}

void matrix_2d::clearlower() {
  for (index_t col = 0, row = 1; col < _mem_cols; ++col, ++row) {
    memset(getelementref(row, col), 0, (_mem_rows - row) * sizeof(double));
  }
}

void matrix_2d::filllower() {
  for (index_t row = 1; row < _rows; row++) {
    for (index_t column = 0; column < row; column++) {
      put(row, column, get(column, row));
    }
  }
}

void matrix_2d::fillupper() {
  for (index_t row = 1; row < _rows; row++) {
    for (index_t column = 0; column < row; column++) {
      put(column, row, get(row, column));
    }
  }
}

void matrix_2d::zero() { memset(_buffer, 0, buffersize()); }

void matrix_2d::zero(const index_t& row_begin, const index_t& col_begin,
                     const index_t& rows, const index_t& columns) {
  for (index_t col = col_begin; col < col_begin + columns; ++col) {
    memset(getelementref(row_begin, col), 0, rows * sizeof(double));
  }
}

matrix_2d matrix_2d::operator=(const matrix_2d& rhs) {
  if (this == &rhs) {
    return *this;
  }
  if (_mem_rows >= rhs.rows() && _mem_cols >= rhs.columns()) {
    _rows = rhs.rows();
    _cols = rhs.columns();
    copybuffer(_rows, _cols, rhs);
    _maxvalCol = rhs.maxvalueCol();
    _maxvalRow = rhs.maxvalueRow();
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
  return *this;
}

matrix_2d matrix_2d::operator*(const double& rhs) const {
  matrix_2d m(_rows, _cols);
  for (index_t row = 0; row < _rows; row++) {
    for (index_t column = 0; column < _cols; ++column) {
      m.put(row, column, get(row, column) * rhs);
    }
  }
  return m;
}

matrix_2d matrix_2d::add(const matrix_2d& rhs) {
  if (_rows != rhs.rows() || _cols != rhs.columns()) {
    throw boost::enable_current_exception(
        std::runtime_error("add(): Dimensions are incompatible."));
  }
  for (index_t row = 0; row < _rows; row++) {
    for (index_t column = 0; column < _cols; ++column) {
      *getelementref(row, column) += rhs.get(row, column);
    }
  }
  return *this;
}

matrix_2d matrix_2d::add(const matrix_2d& lhs, const matrix_2d& rhs) {
  if (lhs.rows() != rhs.rows() || lhs.columns() != rhs.columns()) {
    throw boost::enable_current_exception(
        std::runtime_error("add(): Dimensions are incompatible."));
  }
  *this = lhs;
  return add(rhs);
}

matrix_2d matrix_2d::multiply_mkl(const char* lhs_trans, const matrix_2d& rhs,
                                  const char* rhs_trans) {
  matrix_2d m(_rows, rhs.columns());
  lapack_int lhs_rows = static_cast<lapack_int>(rows());
  lapack_int rhs_cols = static_cast<lapack_int>(rhs.columns());
  lapack_int lhs_cols = static_cast<lapack_int>(columns());
  lapack_int rhs_rows = static_cast<lapack_int>(rhs.rows());

  CBLAS_TRANSPOSE tA =
      (strcmp(lhs_trans, "T") == 0) ? CblasTrans : CblasNoTrans;
  CBLAS_TRANSPOSE tB =
      (strcmp(rhs_trans, "T") == 0) ? CblasTrans : CblasNoTrans;

  if (tA == CblasTrans) {
    lhs_rows = columns();
    lhs_cols = rows();
  }
  if (tB == CblasTrans) {
    rhs_rows = rhs.columns();
    rhs_cols = rhs.rows();
  }
  if (lhs_cols != rhs_rows) {
    throw boost::enable_current_exception(
        std::runtime_error("multiply_mkl(): Incompatible matrix dimensions."));
  }
  if (m.rows() != (index_t)lhs_rows || m.columns() != (index_t)rhs_cols) {
    throw boost::enable_current_exception(std::runtime_error(
        "multiply_mkl(): Result matrix dimensions mismatch."));
  }

#if defined(DNAMATRIX_ROW_WISE)
  cblas_dgemm(CblasRowMajor, tA, tB, lhs_rows, rhs_cols, lhs_cols, 1.0, _buffer,
              _mem_cols, rhs.getbuffer(), rhs.memColumns(), 0.0, m.getbuffer(),
              m.memColumns());
#else
  cblas_dgemm(CblasColMajor, tA, tB, lhs_rows, rhs_cols, lhs_cols, 1.0, _buffer,
              _mem_rows, rhs.getbuffer(), rhs.memRows(), 0.0, m.getbuffer(),
              m.memRows());
#endif

  // return (*this = m);
  return m;
}

matrix_2d matrix_2d::multiply_mkl(const matrix_2d& lhs, const char* lhs_trans,
                                  const matrix_2d& rhs, const char* rhs_trans) {
  lapack_int lhs_rows = static_cast<lapack_int>(lhs.rows());
  lapack_int rhs_cols = static_cast<lapack_int>(rhs.columns());
  lapack_int lhs_cols = static_cast<lapack_int>(lhs.columns());
  lapack_int rhs_rows = static_cast<lapack_int>(rhs.rows());

  CBLAS_TRANSPOSE tA =
      (strncmp(lhs_trans, "T", 1) == 0) ? CblasTrans : CblasNoTrans;
  CBLAS_TRANSPOSE tB =
      (strncmp(rhs_trans, "T", 1) == 0) ? CblasTrans : CblasNoTrans;

  if (tA == CblasTrans) {
    lhs_rows = lhs.columns();
    lhs_cols = lhs.rows();
  }
  if (tB == CblasTrans) {
    rhs_rows = rhs.columns();
    rhs_cols = rhs.rows();
  }
  if (lhs_cols != rhs_rows) {
    throw boost::enable_current_exception(
        std::runtime_error("multiply_mkl(): Incompatible matrix dimensions."));
  }
  if (_rows != (index_t)lhs_rows || _cols != (index_t)rhs_cols) {
    throw boost::enable_current_exception(std::runtime_error(
        "multiply_mkl(): Result matrix dimensions mismatch."));
  }

#if defined(DNAMATRIX_ROW_WISE)
  cblas_dgemm(CblasRowMajor, tA, tB, lhs_rows, rhs_cols, lhs_cols, 1.0,
              lhs.getbuffer(), lhs.memColumns(), rhs.getbuffer(),
              rhs.memColumns(), 0.0, _buffer, _mem_cols);
#else
  cblas_dgemm(CblasColMajor, tA, tB, lhs_rows, rhs_cols, lhs_cols, 1.0,
              lhs.getbuffer(), lhs.memRows(), rhs.getbuffer(), rhs.memRows(),
              0.0, _buffer, _mem_rows);
#endif

  return *this;
}

matrix_2d matrix_2d::transpose(const matrix_2d& matA) {
  if ((matA.columns() != _rows) || (matA.rows() != _cols)) {
    throw boost::enable_current_exception(
        std::runtime_error("transpose(): Dimensions incompatible."));
  }
  for (index_t row = 0; row < _rows; row++) {
    for (index_t column = 0; column < _cols; column++) {
      *getelementref(row, column) = matA.get(column, row);
    }
  }
  return *this;
}

matrix_2d matrix_2d::transpose() {
  matrix_2d m(_cols, _rows);
  for (index_t row = 0; row < _rows; row++) {
    for (index_t column = 0; column < _cols; column++) {
      m.put(column, row, get(row, column));
    }
  }
  return m;
}

double matrix_2d::compute_maximum_value() {
  _maxvalCol = _maxvalRow = 0;
  for (index_t row = 0; row < _rows; ++row) {
    for (index_t col = 0; col < _cols; col++) {
      if (std::fabs(get(row, col)) > std::fabs(get(_maxvalRow, _maxvalCol))) {
        _maxvalCol = col;
        _maxvalRow = row;
      }
    }
  }
  return get(_maxvalRow, _maxvalCol);
}

#ifdef _MSDEBUG
void matrix_2d::trace(const std::string& comment,
                      const std::string& format) const {
  if (comment.empty()) {
    TRACE("%d %d\n", _rows, _cols);
  } else {
    TRACE("%s (%d, %d):\n", comment.c_str(), _rows, _cols);
  }
  for (index_t i = 0; i < _rows; ++i) {
    for (index_t j = 0; j < _cols; ++j) {
      TRACE(format.c_str(), get(i, j));
    }
    TRACE("\n");
  }
  TRACE("\n");
}
#endif

#ifdef _MSDEBUG
void matrix_2d::trace(const std::string& comment,
                      const std::string& submat_comment,
                      const std::string& format, const index_t& row_begin,
                      const index_t& col_begin, const index_t& rows,
                      const index_t& columns) const {
  if (row_begin >= _rows) {
    TRACE("%d %d outside matrix (%d %d)\n", row_begin, col_begin, _rows, _cols);
    return;
  }
  if (col_begin >= _cols) {
    TRACE("%d %d outside matrix (%d %d)\n", row_begin, col_begin, _rows, _cols);
    return;
  }
  if (row_begin + rows > _rows) {
    TRACE("%d %d outside matrix (%d %d)\n", row_begin, col_begin, _rows, _cols);
    return;
  }
  if (col_begin + columns > _cols) {
    TRACE("%d %d outside matrix (%d %d)\n", row_begin, col_begin, _rows, _cols);
    return;
  }
  if (comment.empty()) {
    TRACE("%d %d, %s submatrix (%d, %d, %d*%d)\n", _rows, _cols,
          submat_comment.c_str(), row_begin, col_begin, rows, columns);
  } else {
    TRACE("%s (%d, %d), %s submatrix (%d, %d, %d*%d)\n", comment.c_str(), _rows,
          _cols, submat_comment.c_str(), row_begin, col_begin, rows, columns);
  }

  index_t i, j;
  index_t row_end = row_begin + rows, col_end = col_begin + columns;
  for (i = row_begin; i < row_end; ++i) {
    for (j = col_begin; j < col_end; ++j) {
      TRACE(format.c_str(), get(i, j));
    }
    TRACE("\n");
  }
  TRACE("\n");
}
#endif

} // namespace math
} // namespace dynadjust
