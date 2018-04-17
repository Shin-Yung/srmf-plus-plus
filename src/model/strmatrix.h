/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef STRMATRIX_H
#define STRMATRIX_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

class strMatrix 
{
public:
  using data_t = std::string;
  using row_t = std::vector<data_t>;
  strMatrix() {}
  strMatrix(const int& rows, const int& cols); 
  strMatrix(const data_t& expr); 
  strMatrix(const row_t& vec); 
  ~strMatrix() {}
  void clear(void);
  void resize(const int& rows, const int& cols);
  void getfromtxt(const std::string& file);
  const int& rows(void) const { return rows_; }
  const int& cols(void) const { return cols_; }
  int size(void) const { return rows_ * cols_; }
  strMatrix& operator=(const data_t& expr); 
  strMatrix& operator=(const strMatrix& expr_mat); 
  data_t& operator()(const int& row, const int& col);
  const data_t& operator()(const int& row, const int& col) const;
  row_t& operator[](const int& row);
  const row_t& operator[](const int& row) const;
  friend std::ostream& operator<<(std::ostream& os, const strMatrix& mat);
private:
  std::vector<std::vector<data_t> > matrix_;
  int rows_{0};
  int cols_{0};
}; 

#endif
