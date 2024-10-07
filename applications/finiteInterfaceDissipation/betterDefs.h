using scalarValue = dealii::VectorizedArray<double>;
using scalarGrad = dealii::Tensor<1, dim, dealii::VectorizedArray<double>>;
#define constV(a) dealii::make_vectorized_array(a)