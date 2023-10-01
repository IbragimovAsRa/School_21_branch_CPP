#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : S21Matrix(SIZE_DEFAULT_MATRIX, SIZE_DEFAULT_MATRIX) {}

S21Matrix::~S21Matrix() { FreeMatrix(); }

S21Matrix::S21Matrix(int rows, int cols)
    : rows_(rows), cols_(cols), matrix_(nullptr) {
  try {
    if (rows <= 0 || cols <= 0) {
      throw std::runtime_error(
          "S21Matrix(int rows, int cols)- Некоректный размер матрицы");
    }
    CreateMatrix();
  } catch (const std::exception &e) {
    this->rows_ = 1;
    this->cols_ = 1;
    CreateMatrix();
    std::cerr << "Exception: " << e.what() << std::endl;
  }
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.rows_), cols_(other.cols_) {
  if (this != &other) {
    CreateMatrix();
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++) {
        this->matrix_[i][j] = other.matrix_[i][j];
      }
    }
  }
}

S21Matrix::S21Matrix(S21Matrix &&other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

int S21Matrix::getRows_() const { return rows_; }

int S21Matrix::getCols_() const { return cols_; }

void S21Matrix::setRows_(int rows) {
  try {
    if (rows <= 0 || this->cols_ <= 0 || this->rows_ <= 0) {
      throw std::runtime_error(
          "setRows_(int rows)- Некоректное количество строк");
    }
    if (this->rows_ < rows) { // увеличить кол-во строк
      double **tmpRows_ = new double *[rows];
      for (int i = 0; i < this->rows_; i++) {
        tmpRows_[i] = matrix_[i];
      }
      for (int i = this->rows_; i < rows; i++) {
        tmpRows_[i] = new double[this->cols_];
        for (int j = 0; j < this->cols_; j++) {
          tmpRows_[i][j] = 0;
        }
      }
      delete[] matrix_;
      this->matrix_ = tmpRows_;
    } else if (this->rows_ > rows) { // уменьшить кол-во строк
      for (int i = rows; i < this->rows_; i++) {
        delete[] this->matrix_[i];
      }
    }
    this->rows_ = rows;
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
  }
}

void S21Matrix::setCols_(int cols) {
  try {
    if (cols <= 0 || this->cols_ <= 0 || this->rows_ <= 0) {
      throw std::runtime_error(
          "setCols_(int rows)- Некоректное количество столбцов");
    }
    if (this->cols_ < cols) { // увеличить кол-во столбцов
      double *tmpCols = NULL;
      for (int i = 0; i < this->rows_; i++) {
        tmpCols = new double[cols];
        for (int j = 0; j < this->cols_; j++) {
          tmpCols[j] = this->matrix_[i][j];
        }
        for (int j = this->cols_; j < cols; j++) {
          tmpCols[j] = 0;
        }
        delete[] this->matrix_[i];
        this->matrix_[i] = tmpCols;
      }
    }
    this->cols_ = cols;
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
  }
}

double S21Matrix::Determinant() {
  double dtrm = 0;
  try {
    if (this->rows_ != this->cols_ || this->cols_ <= 0 || this->rows_ <= 0) {
      throw std::runtime_error("Determinant()- Матрица не является квадратной");
    }
    double tmpDtrm = 0;
    int sign = 1;

    S21Matrix tmpMatrix;

    if (this->rows_ == 1) {
      dtrm = this->matrix_[0][0];
    } else {
      for (int j = 0; j < this->rows_; j++) {
        tmpMatrix = *this;
        tmpMatrix.CutLines(0, j);
        tmpDtrm = tmpMatrix.Determinant();
        dtrm += sign * this->matrix_[0][j] * tmpDtrm;
        sign = -sign;
      }
    }
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
  }
  return dtrm;
}

S21Matrix S21Matrix::InverseMatrix() {
  try {
    if (abs(this->Determinant()) < EPS) {
      throw std::runtime_error("InverseMatrix()- Определитель равен нулю");
    }
    S21Matrix resMatrix(this->CalcComplements().Transpose());
    resMatrix.MulNumber(1.0 / this->Determinant());
    return resMatrix;
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
    return *this;
  }
}

S21Matrix S21Matrix::CalcComplements() {
  try {
    if (this->rows_ != this->cols_) {
      throw std::runtime_error(
          "CalcComplements()- Матрица не является квадратной");
    }
    int sign(1);
    S21Matrix tmpMatrix;
    S21Matrix resMatrix(this->rows_, this->cols_);

    for (int i = 0; i < resMatrix.rows_; i++) {
      for (int j = 0; j < resMatrix.cols_; j++) {
        tmpMatrix = *this;
        tmpMatrix.CutLines(i, j);
        resMatrix.matrix_[i][j] = sign * tmpMatrix.Determinant();
        sign = -sign;
      }
      sign = -sign;
    }
    return resMatrix;
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
    return *this;
  }
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix tr(this->rows_, this->cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      tr.matrix_[i][j] = this->matrix_[j][i];
    }
  }
  return tr;
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  try {
    if (this->cols_ != other.rows_) {
      throw std::runtime_error(
          "MulMatrix(const S21Matrix& other)- число столбцов первой "
          "матрицы не равно числу строк второй матрицы");
    }
    int rows_res_mat = this->rows_;
    int cols_res_mat = other.cols_;

    S21Matrix tmp(*this);
    this->setRows_(rows_res_mat);
    this->setCols_(cols_res_mat);

    this->ToZeros();

    int m = tmp.cols_;
    for (int i = 0; i < rows_res_mat; i++) {
      for (int j = 0; j < cols_res_mat; j++) {
        for (int r = 0; r < m; r++) {
          this->matrix_[i][j] += tmp.matrix_[i][r] * other.matrix_[r][j];
        }
      }
    }
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      this->matrix_[i][j] = this->matrix_[i][j] * num;
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  try {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
      throw std::runtime_error(
          "SubMatrix(const S21Matrix& other)- Несовпадают размерности "
          "матриц");
    }
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        this->matrix_[i][j] = this->matrix_[i][j] - other.matrix_[i][j];
      }
    }
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
  }
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  try {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
      throw std::runtime_error(
          "SumMatrix(const S21Matrix& other)- Несовпадают размерности "
          "матриц");
    }
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        this->matrix_[i][j] = this->matrix_[i][j] + other.matrix_[i][j];
      }
    }
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
  }
}

bool S21Matrix::EqMatrix(const S21Matrix &other) {
  bool result = true;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (std::abs(this->matrix_[i][j] - other.matrix_[i][j]) > EPS) {
        result = false;
        break;
      }
    }
    if (!result)
      break;
  }
  return result;
}

// Перегрузка операторов
S21Matrix S21Matrix::operator+(const S21Matrix &m) {
  S21Matrix resMatrix(*this);
  resMatrix.SumMatrix(m);
  return resMatrix;
}

S21Matrix S21Matrix::operator-(const S21Matrix &m) {
  S21Matrix resMatrix(*this);
  resMatrix.SubMatrix(m);
  return resMatrix;
}

S21Matrix S21Matrix::operator*(const S21Matrix &m) {
  S21Matrix resMatrix(*this);
  resMatrix.MulMatrix(m);
  return resMatrix;
}

S21Matrix S21Matrix::operator*(const double &n) {
  S21Matrix resMatrix(*this);
  resMatrix.MulNumber(n);
  return resMatrix;
}

bool S21Matrix::operator==(const S21Matrix &m) { return this->EqMatrix(m); }

S21Matrix &S21Matrix::operator=(const S21Matrix &m) {
  this->setRows_(m.rows_);
  this->setCols_(m.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      this->matrix_[i][j] = m.matrix_[i][j];
    }
  }
  return *this;
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &m) {
  this->SumMatrix(m);
  return *this;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &m) {
  this->SubMatrix(m);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &m) {
  this->MulMatrix(m);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const double n) {
  this->MulNumber(n);
  return *this;
}

double &S21Matrix::operator()(int i, int j) {
  try {
    if (i >= rows_ || j >= cols_ || i < 0 || j < 0) {
      throw std::runtime_error(
          "operator()(int i, int j)- Индекс за пределами матрицы");
    }
    return this->matrix_[i][j];
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
    return this->matrix_[0][0];
  }
}

// Вспомогательные функции
void S21Matrix::ToZeros() {
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] = 0;
    }
  }
}

void S21Matrix::CreateMatrix() {
  if (rows_ > 0 && cols_ > 0) {
    matrix_ = new double *[rows_];
    for (int i = 0; i < rows_; i++) {
      matrix_[i] = new double[cols_];
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] = 0.0;
      }
    }
  }
}

void S21Matrix::FreeMatrix() {
  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
  }
  this->rows_ = 0;
  this->cols_ = 0;
}

S21Matrix::S21Matrix(int rows, int cols, double *m) : rows_(rows), cols_(cols) {
  CreateMatrix();
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] = m[i * this->cols_ + j];
    }
  }
}

void S21Matrix::CutLines(int cut_row, int cut_cols) {
  int rows_new = this->rows_;
  int cols_new = this->cols_;

  // Если поступает отрицательное число, значит вырезать не надо
  if (cut_row >= 0) {
    rows_new--;
  }
  if (cut_cols >= 0) {
    cols_new--;
  }

  S21Matrix m(rows_new, cols_new);
  int i_m = 0;
  int j_m = 0;

  for (int i = 0; i < this->rows_; i++) {
    if (i != cut_row) {
      for (int j = 0; j < this->cols_; j++) {
        if (j != cut_cols) {
          m.matrix_[i_m][j_m] = this->matrix_[i][j];
          j_m++;
        }
      }
      i_m++;
    }
    j_m = 0;
  }

  *this = m;
}
