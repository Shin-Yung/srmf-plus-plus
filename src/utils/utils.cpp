/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-05 11:49:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-06 17:27:18
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./utils.h"
#include <boost/math/special_functions/sign.hpp>
#include <cassert>

namespace util {

//------------------MK_Statistic class----------------------
MK_Statistic::MK_Statistic(void)
{
  resize(1);
}

MK_Statistic::MK_Statistic(const unsigned& size)
{
  resize(size);
}

MK_Statistic::MK_Statistic(const unsigned& size, const unsigned& maxlen)
  : series_maxlen_{maxlen}
{
  resize(size);
}

void MK_Statistic::resize(const unsigned& size)
{
  data_size_ = size;
  time_series_.clear();
  mk_statistic_.resize(data_size_);
  mk_statistic_.setZero();
  mk_trend_.resize(data_size_);
  if (series_maxlen_==0) is_full_ = true;
  else is_full_ = false;
}

void MK_Statistic::set_maxlen(const unsigned& maxlen)
{
  series_maxlen_ = maxlen;
  reset();
}

void MK_Statistic::reset(void)
{
  time_series_.clear(); 
  mk_statistic_.setZero();
  if (series_maxlen_==0) is_full_ = true;
  else is_full_ = false;
}


const MK_Statistic& MK_Statistic::add_data(const double& x)
{
  assert(data_size_==1);
  data_t data(1);
  data[0] = x;
  add_data(data_t(data));
  return *this;
}

const MK_Statistic& MK_Statistic::add_data(const data_t& y)
{
  assert(data_size_==y.size());
  // if time series length has exceeded limit 
  // (series_maxlen_=0 implies infinite lenth)
  if (series_maxlen_>0 && time_series_.size()>=series_maxlen_) {
    is_full_ = true;
    // remove the data front, update statistic  
    data_t x0(data_size_);
    x0 = time_series_.front();
    for (const auto& x : time_series_) {
      for (unsigned j=0; j<data_size_; ++j)
        mk_statistic_[j] -= boost::math::sign(x[j] - x0[j]);
    }
    time_series_.pop_front();
  }

  // add new data at back, update statistic
  for (const auto& x : time_series_) {
    for (unsigned i=0; i<data_size_; ++i)
      mk_statistic_[i] += boost::math::sign(y[i]-x[i]);
  } 
  time_series_.push_back(data_t(y));
  // Mann-Kendall trend
  int N = time_series_.size();
  /*std::cout << "time series len = " << time_series_.size() << "\n";
  for (const auto& x : time_series_)
  std::cout << "x = " << x << "\n";
  std::cout << "-----------" << "\n";
  */
  if (N>1) {
    int M = N*(N-1);
    for (unsigned i=0; i<data_size_; ++i) {
      mk_trend_[i] = 2.0 * static_cast<double>(std::abs(mk_statistic_[i]))/M;
    } 
  }
  else mk_trend_.setZero();
  mk_trend_max_ = mk_trend_.maxCoeff();
  return *this;
}

void MK_Statistic::get_series_avg(data_t& mean) const
{
  mean.resize(data_size_);
  mean.setZero();
  for (const auto& x : time_series_) mean += x;
  if (time_series_.size()>0) mean /= time_series_.size();
}
//------------------------------------------------------



} // end namespace util

