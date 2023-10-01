#ifndef S21_MATRIX_OOP_H
#define S21_MATRIX_OOP_H

#include <exception>
#include <iostream>

#define SIZE_DEFAULT_MATRIX 3
#define EPS 10e-6

// Добавить обработку ошибок
// Написать тесты

class S21Matrix {
private:
  int rows_;
  int cols_;
  double **matrix_;

public:
  // Конструкторы
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix &other);
  S21Matrix(S21Matrix &&other);

  // Деструктор
  ~S21Matrix();

  // Аксессоры
  int getRows_() const;
  int getCols_() const;

  // Мутаторы
  void setRows_(int rows);
  void setCols_(int cols);

  // Операции
  bool EqMatrix(const S21Matrix &other);
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix &other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  // Перегрузки операторов
  S21Matrix &operator=(const S21Matrix &m);
  S21Matrix operator+(const S21Matrix &m);
  S21Matrix operator-(const S21Matrix &m);
  S21Matrix operator*(const S21Matrix &m);
  S21Matrix operator*(const double &n);

  bool operator==(const S21Matrix &m);
  S21Matrix &operator+=(const S21Matrix &m);
  S21Matrix &operator-=(const S21Matrix &m);
  S21Matrix &operator*=(const S21Matrix &m);
  S21Matrix &operator*=(const double n);
  double &operator()(int i, int j);

  // Вспомогательные методы
  S21Matrix(int rows, int cols, double *m);
  void CreateMatrix();
  void FreeMatrix();
  void ToZeros();
  void CutLines(int cut_row, int cut_cols);
};

#endif // S21_MATRIX_OOP_H
