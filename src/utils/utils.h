/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-05 11:48:05
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-09 17:48:17
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <exception>
#include <deque>
#include <Eigen/Core>

namespace util {

// Mann-Kendall statisic of a time series
class MK_Statistic
{
public:
  using data_t = Eigen::VectorXd;
  MK_Statistic();
  MK_Statistic(const unsigned& size); 
  MK_Statistic(const unsigned& size, const unsigned& max_len); 
  ~MK_Statistic() {}
  void reset(void);
  void resize(const unsigned& size);
  void set_maxlen(const unsigned& maxlen);
  const MK_Statistic& add_data(const double& x);
  const MK_Statistic& add_data(const data_t& x);
  const MK_Statistic& operator<<(const double& x) { return add_data(x); }
  const MK_Statistic& operator<<(const data_t& x) { return add_data(x); }
  void get_series_avg(data_t& mean) const;
  const double& elem_max_trend(void) const { return mk_trend_max_; }
  const bool& is_full(void) const { return is_full_; }
private:
  unsigned data_size_{1};
  unsigned series_maxlen_{0};
  bool is_full_{false};
  std::deque<data_t> time_series_;
  Eigen::VectorXi mk_statistic_;
  Eigen::VectorXd mk_trend_;
  double mk_trend_max_{0};
};


} // end namespace util

#endif