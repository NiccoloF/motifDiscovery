#include "FunBialign.hpp"

FunBialign::FunBialign(const KMA::matrix& data,
                       unsigned int portion_len):_outcols(portion_len),
                                                 _numerosity(data.n_rows,data.n_cols - portion_len + 1)
{
  const unsigned int totdim = data.n_cols;
  const unsigned int totobs = data.n_rows;
  const unsigned int totrows = (totdim - portion_len + 1) * totobs;
  const unsigned int totportion = totdim - portion_len + 1;
  
  // set the size for data structure
  _windowData.set_size(totrows,portion_len);
  if(totobs <= 1) // we have a single curve
  {
    _isSingle = true;
    // fill in my data 
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(arma::uword i = 0; i < totrows; ++i) 
      _windowData.row(i) = data.cols(i,i + portion_len - 1); //Each peace of curve is stored in a row
  }
  else // we have multiple curves
  {
    _isSingle = false;
    // fill in my data
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
      for(arma::uword k = 0; k < totobs; ++k)
        for (arma::uword i = 0; i < totportion; ++i)
          _windowData.row(totportion * k + i) = data(k,arma::span(i,i + portion_len - 1));
  }
  
  
  // Find indexes of rows containing at least one NaN
  const unsigned int& total_rows = _windowData.n_rows;
  arma::uvec NA_rows = arma::unique(arma::conv_to<arma::uvec>::from
                       (arma::find_nonfinite(_windowData)).transform
                       ([&total_rows](const arma::uword row)
                       {return row%total_rows;}));
  
  // Remove rows with NaN 
  _windowData.shed_rows(NA_rows);
  
  // compute the correct row 
  NA_rows.transform([&totportion](unsigned int row){return std::floor(row / totportion);});
  
  //counts number of occurrences for each curve
  for(const unsigned int row: NA_rows)
    _numerosity[row]--;
  
  // save number of rows
  _outrows = _windowData.n_rows;
  
  _scoreData.set_size(_outrows * (_outrows - 1)/2);
}




KMA::vector FunBialign::Hscore() const
{
  double outcols_inv = 1.0 / static_cast<double>(_outcols);
  const arma::colvec& vsum = arma::sum(_windowData,1);
  
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
  for (arma::uword i = 1; i < _outrows; ++i) { // for any functional observation
    const arma::rowvec& x_i = _windowData.row(i); // curve i 
    for (arma::uword j = 0; j < i; ++j) {
      const arma::rowvec& x_j = _windowData.row(j);
      const arma::rowvec& crossSum = x_i + x_j;
      
      // Precompute common terms
      const double& crossSumAccu = arma::accu(crossSum);
      const double& commonTerm = crossSumAccu * (outcols_inv * 0.5);
      const double& x_i_sum_term = vsum[i] * outcols_inv;
    
       // map from the lower triangular part to vector
       arma::uword index = j*(2*_outrows - j - 1)/2 + i - j - 1;
      _scoreData(index) = arma::accu(arma::square(
                                        x_i - x_i_sum_term - crossSum*0.5 
                                        + commonTerm))*outcols_inv;
    }
  }
  return _scoreData;
}


double FunBialign::fMSR_adj() const
{  
  if(!_outrows) // degenerate case with zero curves
    return 0.0;
  
  double score = arma::accu(arma::square(
                 (_windowData.each_col() - arma::mean(_windowData, 1)).each_row() - arma::mean(_windowData,0)
                 + arma::mean(arma::vectorise(_windowData))))/static_cast<double>(_outrows * _outcols);
  if (_outrows > 2) 
  {
    score /= arma::as_scalar(arma::prod(arma::regspace<KMA::vector>(2, 1, _outrows).transform(
      [](double k){return k * k / (k * k - 1);})));
  }
  return score;
}


KMA::vector FunBialign::createDistance() const
{
  this -> Hscore();

  // Modify accolites 
  double M = static_cast<double>(*std::max_element(_scoreData.begin(),_scoreData.end())) + 1000.0; // very large distance for accolites
  int overlap = static_cast<int>(std::floor(_outcols / 2.0)); // number of right/left accolites
  
  // Single curve
  if (_isSingle) {
    for (int i = 0; i < _outrows-1; ++i) {
      // check for right accolites
      arma::uword start = i+1;
      arma::uword end = std::min<int>(i+overlap,_outrows-1);
      start += i*(2*_outrows - i -1)/2 - i - 1; 
      end += i*(2*_outrows - i -1)/2 - i - 1; 
      _scoreData(arma::span(start,end)) += M;
    }
    
    /* TODO
    auto range = std::views::iota(0,static_cast<int>(_outrows) - 1);
    std::for_each(std::execution::par,range.begin(),range.end(),[&overlap,&M,this](arma::uword i)
      {arma::uword start = i+1;
       arma::uword end = std::min<int>(i+overlap,_outrows-1);
       start += i*(2*_outrows - i -1)/2 - i - 1; 
       end += i*(2*_outrows - i -1)/2 - i - 1; 
      _scoreData(arma::span(start,end)) += M;});
   */
    
  // Multiple curves
  } else {
    int until_here = 0;
    for(int j = 0; j <_numerosity.size();++j)
    {
      int num = _numerosity[j] - 1;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
      for(int i = until_here; i < until_here + _numerosity[j]-1;++i)
      {
        arma::uword start = i * (2*_outrows - i -1)/2 + (i+1) - i -1;
        arma::uword k = std::min<int>(std::min<int>(i+overlap ,i + (num--)),_outrows);
        arma::uword end = i * (2*_outrows - i -1)/2 + k - i -1;
        _scoreData(arma::span(start,end)) += M;
      }
      until_here += _numerosity[j];
    }
  }
  return _scoreData;
}




RCPP_EXPOSED_CLASS(FunBialign);

RCPP_MODULE(FunBialignModule) {
  Rcpp::class_<FunBialign>("FunBialign")
  .constructor<KMA::matrix,unsigned int>()
  .method("fMSR",&FunBialign::fMSR_adj)
  .method("createDistance",&FunBialign::createDistance);
}

