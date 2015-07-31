#include "DataSet_MatrixFlt.h"
void DataSet_MatrixFlt::WriteBuffer(CpptrajFile& outfile, SizeArray const& pIn) const {
  size_t x = (size_t)pIn[0];
  size_t y = (size_t)pIn[1];
  if ( x >= mat_.Ncols() || y >= mat_.Nrows() )
    outfile.Printf(data_format_, 0.0);
  else 
    outfile.Printf(data_format_, mat_.element(x,y));
}

double* DataSet_MatrixFlt::MatrixArray() const {
  double* matOut = new double[ mat_.size() ];
  for (size_t i = 0; i < mat_.size(); ++i)
    matOut[i] = (double)mat_[i];
  return matOut;
}
