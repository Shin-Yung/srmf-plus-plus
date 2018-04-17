/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* @Date:   2018-04-08 17:11:15
* @Last Modified by:   Amal Medhi, amedhi@macbook
* @Last Modified time: 2018-04-09 12:54:12
*----------------------------------------------------------------------------*/
#include "./strmatrix.h"
#include <fstream>
#include <sstream> 
#include <stdexcept>

strMatrix::strMatrix(const int& rows, const int& cols)
{
  matrix_.resize(rows);
  for (auto& rvec : matrix_) rvec.resize(cols);
  rows_ = rows;
  cols_ = cols;
}

strMatrix::strMatrix(const data_t& expr)
{
  rows_ = 1;
  cols_ = 1;
  matrix_.resize(1);
  matrix_[0].resize(1);
  matrix_[0][0] = expr;
}

strMatrix::strMatrix(const strMatrix::row_t& vec)
{
  rows_ = 1;
  cols_ = vec.size();
  matrix_.resize(rows_);
  matrix_[0] = vec;
}

void strMatrix::clear(void) 
{
  matrix_.clear();
  rows_ = 0;
  cols_ = 0;
}

void strMatrix::resize(const int& rows, const int& cols)
{
  if (rows_==rows && cols_==cols) return;
  matrix_.resize(rows);
  for (auto& rvec : matrix_) rvec.resize(cols);
  rows_ = rows;
  cols_ = cols;
}

void strMatrix::getfromtxt(const std::string& file)
{
  std::ifstream fs(file);
  if (!fs.is_open()) {
    throw std::runtime_error("strMatrix::getfromtxt: file open failed");
  }
  std::string word;
  std::string line;
  matrix_.clear();
  while (std::getline(fs, line)) {
    // # is a comment character
    line = line.substr(0, line.find('#'));
    // skip blank line
    if (line.size() == 0) continue;
    std::size_t first = line.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) continue;

    std::istringstream ss(line);
    strMatrix::row_t row_vector;
    while (ss >> word) row_vector.push_back(word);
    if (matrix_.size() > 0) {
      if (row_vector.size() != matrix_.back().size()) {
        throw std::runtime_error("strMatrix::getfromtxt: unexpected data format in file");
      }
    }
    matrix_.push_back(row_vector);
  }
  rows_ = matrix_.size();
  cols_ = matrix_[0].size();
  fs.close();
}

strMatrix& strMatrix::operator=(const data_t& expr)
{
  rows_ = 1;
  cols_ = 1;
  matrix_.resize(1);
  matrix_[0].resize(1);
  matrix_[0][0] = expr;
  return *this;
}

strMatrix& strMatrix::operator=(const strMatrix& expr_mat)
{
  rows_ = expr_mat.rows();
  cols_ = expr_mat.cols();
  matrix_.resize(rows_);
  for (unsigned i=0; i<rows_; ++i)
    matrix_[i] = expr_mat[i];
  return *this;
}

strMatrix::data_t& strMatrix::operator()(const int& row, const int& col)
{
  return matrix_[row][col];
}

const strMatrix::data_t& strMatrix::operator()(const int& row, const int& col) const
{
  return matrix_[row][col];
}

strMatrix::row_t& strMatrix::operator[](const int& row)
{
  return matrix_[row];
}

const strMatrix::row_t& strMatrix::operator[](const int& row) const
{
  return matrix_[row];
}

std::ostream& operator<<(std::ostream& os, const strMatrix& mat)
{
  int maxlen = 0;
  for (int i=0; i<mat.rows(); ++i) {
    for (int j=0; j<mat.cols(); ++j) {
      int len = mat(i,j).size();
      if (maxlen <= len) maxlen = len;
    }
  }
  os << std::left;
  for (int i=0; i<mat.rows(); ++i) {
    for (int j=0; j<mat.cols(); ++j) {
      os << std::setw(maxlen) << mat(i,j) << "   ";
    }
    os << "\n";
  }
  os << std::right;
  return os;
}
