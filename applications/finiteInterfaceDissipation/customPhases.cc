// This file will be created or modified by AMMBER
#include <deal.II/base/vectorization.h>

#include "SystemContainer.h"

dealii::VectorizedArray<double> T = dealii::make_vectorized_array(320.);
dealii::VectorizedArray<double> R = dealii::make_vectorized_array(8.314);

// Define a function object to wrap the template
struct intpow
{
  template <int N, typename T>
  T
  operator()(const T t) const
  {
    return dealii::Utilities::fixed_power<N>(t);
  }
};

// Individual phases are derived classes of PhaseFieldContainer
template <int dim, int degree>
class LIQUID : public PhaseFieldContainer<dim, degree>
{
public:
  LIQUID(const IsothermalSystem                                          &isoSys,
         const std::string                                               &phase_name,
         const std::map<std::string, PhaseFieldContainer<dim, degree> *> &phase_fields)
    : PhaseFieldContainer<dim, degree>(isoSys, phase_name, phase_fields)
  {}

  dealii::VectorizedArray<double> GaIn0a = constV(5.14148219e+03);
  dealii::VectorizedArray<double> GaIn0b = constV(-7.89869349e-01);
  dealii::VectorizedArray<double> GaIn1a = constV(3.14491130e+03);
  dealii::VectorizedArray<double> GaIn1b = constV(-9.80783643e+00);

  dealii::VectorizedArray<double> InPb0a = constV(3.775e+03);
  dealii::VectorizedArray<double> InPb0b = constV(-1.285e+00);
  dealii::VectorizedArray<double> InPb1a = constV(1.830e+02);
  dealii::VectorizedArray<double> InPb1b = constV(3.810e-01);

  dealii::VectorizedArray<double> PbGa0a = constV(2.28142140e+04);
  dealii::VectorizedArray<double> PbGa0b = constV(-8.99292499e+00);
  dealii::VectorizedArray<double> PbGa1a = constV(-4.90642270e+03);
  dealii::VectorizedArray<double> PbGa1b = constV(5.96724182e+00);

  inline void
  calculate_free_energy() override
  {
    const FieldContainer<dim> &x_A = this->comp_data["A"].x_data;
    const FieldContainer<dim> &x_B = this->comp_data["B"].x_data;

    FieldContainer<dim> &dfdx_A = this->comp_data["A"].dfdx;
    FieldContainer<dim> &dfdx_B = this->comp_data["B"].dfdx;

    this->phase_free_energy =
      R * T *
        (x_A.val * log(x_A.val) + x_B.val * log(x_B.val) +
         (-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0)) *
           log(-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) +
      4.0 * x_A.val * x_B.val *
        (GaIn0a + GaIn0b * T + (GaIn1a + GaIn1b * T) * (x_A.val - 1.0 * x_B.val)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) -
          1.0 * dealii::Utilities::fixed_power<2>(x_A.val - 1.0 * x_B.val)) +
      x_A.val *
        (PbGa0a + PbGa0b * T +
         (PbGa1a + PbGa1b * T) * (-2.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) *
        (-4.0 * x_A.val - 4.0 * x_B.val + constV(4.0)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(
                                -1.0 * x_A.val - 0.5 * x_B.val + constV(0.5))) +
      4.0 * x_B.val *
        (InPb0a + InPb0b * T +
         (InPb1a + InPb1b * T) * (x_A.val + 2.0 * x_B.val + constV(-1.0))) *
        (-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(0.5 * x_A.val + x_B.val +
                                                                constV(-0.5)));

    dfdx_A.val =
      R * T * (log(x_A.val) - 1.0 * log(-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) +
      4.0 * x_A.val * x_B.val * (GaIn1a + GaIn1b * T) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) -
          1.0 * dealii::Utilities::fixed_power<2>(x_A.val - 1.0 * x_B.val)) +
      4.0 * x_A.val * x_B.val * (2.0 * x_A.val - 2.0 * x_B.val) *
        (GaIn0a + GaIn0b * T + (GaIn1a + GaIn1b * T) * (x_A.val - 1.0 * x_B.val)) *
        dealii::Utilities::fixed_power<-2>(
          constV(1.0) -
          1.0 * dealii::Utilities::fixed_power<2>(x_A.val - 1.0 * x_B.val)) +
      x_A.val * (-2.0 * PbGa1a - 2.0 * PbGa1b * T) *
        (-4.0 * x_A.val - 4.0 * x_B.val + constV(4.0)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(
                                -1.0 * x_A.val - 0.5 * x_B.val + constV(0.5))) +
      0.0625 * x_A.val *
        (PbGa0a + PbGa0b * T +
         (PbGa1a + PbGa1b * T) * (-2.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) *
        (-4.0 * x_A.val - 4.0 * x_B.val + constV(4.0)) *
        (8.0 * x_A.val + 4.0 * x_B.val + constV(-4.0)) *
        dealii::Utilities::fixed_power<-2>(
          constV(0.25) - 1.0 * dealii::Utilities::fixed_power<2>(
                                 -1.0 * x_A.val - 0.5 * x_B.val + constV(0.5))) -
      4.0 * x_A.val *
        (PbGa0a + PbGa0b * T +
         (PbGa1a + PbGa1b * T) * (-2.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(
                                -1.0 * x_A.val - 0.5 * x_B.val + constV(0.5))) +
      4.0 * x_B.val * (InPb1a + InPb1b * T) *
        (-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(0.5 * x_A.val + x_B.val +
                                                                constV(-0.5))) +
      4.0 * x_B.val *
        (GaIn0a + GaIn0b * T + (GaIn1a + GaIn1b * T) * (x_A.val - 1.0 * x_B.val)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) -
          1.0 * dealii::Utilities::fixed_power<2>(x_A.val - 1.0 * x_B.val)) +
      0.25 * x_B.val *
        (InPb0a + InPb0b * T +
         (InPb1a + InPb1b * T) * (x_A.val + 2.0 * x_B.val + constV(-1.0))) *
        (-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0)) *
        (2.0 * x_A.val + 4.0 * x_B.val + constV(-2.0)) *
        dealii::Utilities::fixed_power<-2>(
          constV(0.25) - 1.0 * dealii::Utilities::fixed_power<2>(0.5 * x_A.val + x_B.val +
                                                                 constV(-0.5))) -
      4.0 * x_B.val *
        (InPb0a + InPb0b * T +
         (InPb1a + InPb1b * T) * (x_A.val + 2.0 * x_B.val + constV(-1.0))) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(0.5 * x_A.val + x_B.val +
                                                                constV(-0.5))) +
      (PbGa0a + PbGa0b * T +
       (PbGa1a + PbGa1b * T) * (-2.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) *
        (-4.0 * x_A.val - 4.0 * x_B.val + constV(4.0)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(
                                -1.0 * x_A.val - 0.5 * x_B.val + constV(0.5)));

    dfdx_B.val =
      R * T * (log(x_B.val) - 1.0 * log(-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) +
      4.0 * x_A.val * x_B.val * (-1.0 * GaIn1a - 1.0 * GaIn1b * T) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) -
          1.0 * dealii::Utilities::fixed_power<2>(x_A.val - 1.0 * x_B.val)) +
      4.0 * x_A.val * x_B.val * (-2.0 * x_A.val + 2.0 * x_B.val) *
        (GaIn0a + GaIn0b * T + (GaIn1a + GaIn1b * T) * (x_A.val - 1.0 * x_B.val)) *
        dealii::Utilities::fixed_power<-2>(
          constV(1.0) -
          1.0 * dealii::Utilities::fixed_power<2>(x_A.val - 1.0 * x_B.val)) +
      x_A.val * (-1.0 * PbGa1a - 1.0 * PbGa1b * T) *
        (-4.0 * x_A.val - 4.0 * x_B.val + constV(4.0)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(
                                -1.0 * x_A.val - 0.5 * x_B.val + constV(0.5))) +
      4.0 * x_A.val *
        (GaIn0a + GaIn0b * T + (GaIn1a + GaIn1b * T) * (x_A.val - 1.0 * x_B.val)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) -
          1.0 * dealii::Utilities::fixed_power<2>(x_A.val - 1.0 * x_B.val)) +
      0.0625 * x_A.val *
        (PbGa0a + PbGa0b * T +
         (PbGa1a + PbGa1b * T) * (-2.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) *
        (-4.0 * x_A.val - 4.0 * x_B.val + constV(4.0)) *
        (4.0 * x_A.val + 2.0 * x_B.val + constV(-2.0)) *
        dealii::Utilities::fixed_power<-2>(
          constV(0.25) - 1.0 * dealii::Utilities::fixed_power<2>(
                                 -1.0 * x_A.val - 0.5 * x_B.val + constV(0.5))) -
      4.0 * x_A.val *
        (PbGa0a + PbGa0b * T +
         (PbGa1a + PbGa1b * T) * (-2.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(
                                -1.0 * x_A.val - 0.5 * x_B.val + constV(0.5))) +
      4.0 * x_B.val * (2.0 * InPb1a + 2.0 * InPb1b * T) *
        (-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(0.5 * x_A.val + x_B.val +
                                                                constV(-0.5))) +
      0.25 * x_B.val *
        (InPb0a + InPb0b * T +
         (InPb1a + InPb1b * T) * (x_A.val + 2.0 * x_B.val + constV(-1.0))) *
        (-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0)) *
        (4.0 * x_A.val + 8.0 * x_B.val + constV(-4.0)) *
        dealii::Utilities::fixed_power<-2>(
          constV(0.25) - 1.0 * dealii::Utilities::fixed_power<2>(0.5 * x_A.val + x_B.val +
                                                                 constV(-0.5))) -
      4.0 * x_B.val *
        (InPb0a + InPb0b * T +
         (InPb1a + InPb1b * T) * (x_A.val + 2.0 * x_B.val + constV(-1.0))) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(0.5 * x_A.val + x_B.val +
                                                                constV(-0.5))) +
      4.0 *
        (InPb0a + InPb0b * T +
         (InPb1a + InPb1b * T) * (x_A.val + 2.0 * x_B.val + constV(-1.0))) *
        (-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(0.5 * x_A.val + x_B.val +
                                                                constV(-0.5)));

    dfdx_A.grad = (2.0) * x_A.grad + (0.0) * x_B.grad;
    dfdx_B.grad = (0.0) * x_A.grad + (2.0) * x_B.grad;

    this->volumetrize_free_energy();
    this->nondimensionalize_free_energy();
  }
};

// Individual phases are derived classes of PhaseFieldContainer
template <int dim, int degree>
class SOLID : public PhaseFieldContainer<dim, degree>
{
public:
  SOLID(const IsothermalSystem                                          &isoSys,
        const std::string                                               &phase_name,
        const std::map<std::string, PhaseFieldContainer<dim, degree> *> &phase_fields)
    : PhaseFieldContainer<dim, degree>(isoSys, phase_name, phase_fields)
  {}

  dealii::VectorizedArray<double> GaIn0a = dealii::make_vectorized_array(0.);
  dealii::VectorizedArray<double> GaIn0b = dealii::make_vectorized_array(0.);
  dealii::VectorizedArray<double> GaIn1a = dealii::make_vectorized_array(0.);
  dealii::VectorizedArray<double> GaIn1b = dealii::make_vectorized_array(0.);

  dealii::VectorizedArray<double> InPb0a = dealii::make_vectorized_array(5069.);
  dealii::VectorizedArray<double> InPb0b = dealii::make_vectorized_array(-2.624);
  dealii::VectorizedArray<double> InPb1a = dealii::make_vectorized_array(456.);
  dealii::VectorizedArray<double> InPb1b = dealii::make_vectorized_array(0.875);

  dealii::VectorizedArray<double> PbGa0a = dealii::make_vectorized_array(0.);
  dealii::VectorizedArray<double> PbGa0b = dealii::make_vectorized_array(0.);
  dealii::VectorizedArray<double> PbGa1a = dealii::make_vectorized_array(0.);
  dealii::VectorizedArray<double> PbGa1b = dealii::make_vectorized_array(0.);

  inline void
  calculate_free_energy() override
  {
    const FieldContainer<dim> &x_A = this->comp_data["A"].x_data;
    const FieldContainer<dim> &x_B = this->comp_data["B"].x_data;

    FieldContainer<dim> &dfdx_A = this->comp_data["A"].dfdx;
    FieldContainer<dim> &dfdx_B = this->comp_data["B"].dfdx;

    this->phase_free_energy =
      R * T *
        (x_A.val * log(x_A.val) + x_B.val * log(x_B.val) +
         (-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0)) *
           log(-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) +
      4.0 * x_A.val * x_B.val *
        (GaIn0a + GaIn0b * T + (GaIn1a + GaIn1b * T) * (x_A.val - 1.0 * x_B.val)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) -
          1.0 * dealii::Utilities::fixed_power<2>(x_A.val - 1.0 * x_B.val)) +
      x_A.val *
        (PbGa0a + PbGa0b * T +
         (PbGa1a + PbGa1b * T) * (-2.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) *
        (-4.0 * x_A.val - 4.0 * x_B.val + constV(4.0)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(
                                -1.0 * x_A.val - 0.5 * x_B.val + constV(0.5))) +
      20000.0 * x_A.val + x_B.val * (-0.002 * T + constV(38.0)) +
      4.0 * x_B.val *
        (InPb0a + InPb0b * T +
         (InPb1a + InPb1b * T) * (x_A.val + 2.0 * x_B.val + constV(-1.0))) *
        (-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(0.5 * x_A.val + x_B.val +
                                                                constV(-0.5))) +
      (-1.183 * T + constV(967.0)) * (-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0));

    dfdx_A.val =
      R * T * (log(x_A.val) - 1.0 * log(-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) +
      1.183 * T +
      4.0 * x_A.val * x_B.val * (GaIn1a + GaIn1b * T) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) -
          1.0 * dealii::Utilities::fixed_power<2>(x_A.val - 1.0 * x_B.val)) +
      4.0 * x_A.val * x_B.val * (2.0 * x_A.val - 2.0 * x_B.val) *
        (GaIn0a + GaIn0b * T + (GaIn1a + GaIn1b * T) * (x_A.val - 1.0 * x_B.val)) *
        dealii::Utilities::fixed_power<-2>(
          constV(1.0) -
          1.0 * dealii::Utilities::fixed_power<2>(x_A.val - 1.0 * x_B.val)) +
      x_A.val * (-2.0 * PbGa1a - 2.0 * PbGa1b * T) *
        (-4.0 * x_A.val - 4.0 * x_B.val + constV(4.0)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(
                                -1.0 * x_A.val - 0.5 * x_B.val + constV(0.5))) +
      0.0625 * x_A.val *
        (PbGa0a + PbGa0b * T +
         (PbGa1a + PbGa1b * T) * (-2.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) *
        (-4.0 * x_A.val - 4.0 * x_B.val + constV(4.0)) *
        (8.0 * x_A.val + 4.0 * x_B.val + constV(-4.0)) *
        dealii::Utilities::fixed_power<-2>(
          constV(0.25) - 1.0 * dealii::Utilities::fixed_power<2>(
                                 -1.0 * x_A.val - 0.5 * x_B.val + constV(0.5))) -
      4.0 * x_A.val *
        (PbGa0a + PbGa0b * T +
         (PbGa1a + PbGa1b * T) * (-2.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(
                                -1.0 * x_A.val - 0.5 * x_B.val + constV(0.5))) +
      4.0 * x_B.val * (InPb1a + InPb1b * T) *
        (-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(0.5 * x_A.val + x_B.val +
                                                                constV(-0.5))) +
      4.0 * x_B.val *
        (GaIn0a + GaIn0b * T + (GaIn1a + GaIn1b * T) * (x_A.val - 1.0 * x_B.val)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) -
          1.0 * dealii::Utilities::fixed_power<2>(x_A.val - 1.0 * x_B.val)) +
      0.25 * x_B.val *
        (InPb0a + InPb0b * T +
         (InPb1a + InPb1b * T) * (x_A.val + 2.0 * x_B.val + constV(-1.0))) *
        (-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0)) *
        (2.0 * x_A.val + 4.0 * x_B.val + constV(-2.0)) *
        dealii::Utilities::fixed_power<-2>(
          constV(0.25) - 1.0 * dealii::Utilities::fixed_power<2>(0.5 * x_A.val + x_B.val +
                                                                 constV(-0.5))) -
      4.0 * x_B.val *
        (InPb0a + InPb0b * T +
         (InPb1a + InPb1b * T) * (x_A.val + 2.0 * x_B.val + constV(-1.0))) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(0.5 * x_A.val + x_B.val +
                                                                constV(-0.5))) +
      (PbGa0a + PbGa0b * T +
       (PbGa1a + PbGa1b * T) * (-2.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) *
        (-4.0 * x_A.val - 4.0 * x_B.val + constV(4.0)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(
                                -1.0 * x_A.val - 0.5 * x_B.val + constV(0.5))) +
      constV(19033.0);

    dfdx_B.val =
      R * T * (log(x_B.val) - 1.0 * log(-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) +
      1.181 * T +
      4.0 * x_A.val * x_B.val * (-1.0 * GaIn1a - 1.0 * GaIn1b * T) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) -
          1.0 * dealii::Utilities::fixed_power<2>(x_A.val - 1.0 * x_B.val)) +
      4.0 * x_A.val * x_B.val * (-2.0 * x_A.val + 2.0 * x_B.val) *
        (GaIn0a + GaIn0b * T + (GaIn1a + GaIn1b * T) * (x_A.val - 1.0 * x_B.val)) *
        dealii::Utilities::fixed_power<-2>(
          constV(1.0) -
          1.0 * dealii::Utilities::fixed_power<2>(x_A.val - 1.0 * x_B.val)) +
      x_A.val * (-1.0 * PbGa1a - 1.0 * PbGa1b * T) *
        (-4.0 * x_A.val - 4.0 * x_B.val + constV(4.0)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(
                                -1.0 * x_A.val - 0.5 * x_B.val + constV(0.5))) +
      4.0 * x_A.val *
        (GaIn0a + GaIn0b * T + (GaIn1a + GaIn1b * T) * (x_A.val - 1.0 * x_B.val)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) -
          1.0 * dealii::Utilities::fixed_power<2>(x_A.val - 1.0 * x_B.val)) +
      0.0625 * x_A.val *
        (PbGa0a + PbGa0b * T +
         (PbGa1a + PbGa1b * T) * (-2.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) *
        (-4.0 * x_A.val - 4.0 * x_B.val + constV(4.0)) *
        (4.0 * x_A.val + 2.0 * x_B.val + constV(-2.0)) *
        dealii::Utilities::fixed_power<-2>(
          constV(0.25) - 1.0 * dealii::Utilities::fixed_power<2>(
                                 -1.0 * x_A.val - 0.5 * x_B.val + constV(0.5))) -
      4.0 * x_A.val *
        (PbGa0a + PbGa0b * T +
         (PbGa1a + PbGa1b * T) * (-2.0 * x_A.val - 1.0 * x_B.val + constV(1.0))) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(
                                -1.0 * x_A.val - 0.5 * x_B.val + constV(0.5))) +
      4.0 * x_B.val * (2.0 * InPb1a + 2.0 * InPb1b * T) *
        (-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(0.5 * x_A.val + x_B.val +
                                                                constV(-0.5))) +
      0.25 * x_B.val *
        (InPb0a + InPb0b * T +
         (InPb1a + InPb1b * T) * (x_A.val + 2.0 * x_B.val + constV(-1.0))) *
        (-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0)) *
        (4.0 * x_A.val + 8.0 * x_B.val + constV(-4.0)) *
        dealii::Utilities::fixed_power<-2>(
          constV(0.25) - 1.0 * dealii::Utilities::fixed_power<2>(0.5 * x_A.val + x_B.val +
                                                                 constV(-0.5))) -
      4.0 * x_B.val *
        (InPb0a + InPb0b * T +
         (InPb1a + InPb1b * T) * (x_A.val + 2.0 * x_B.val + constV(-1.0))) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(0.5 * x_A.val + x_B.val +
                                                                constV(-0.5))) +
      4.0 *
        (InPb0a + InPb0b * T +
         (InPb1a + InPb1b * T) * (x_A.val + 2.0 * x_B.val + constV(-1.0))) *
        (-1.0 * x_A.val - 1.0 * x_B.val + constV(1.0)) *
        dealii::Utilities::fixed_power<-1>(
          constV(1.0) - 4.0 * dealii::Utilities::fixed_power<2>(0.5 * x_A.val + x_B.val +
                                                                constV(-0.5))) +
      constV(-929.0);

    dfdx_A.grad = (2.0) * x_A.grad + (0.0) * x_B.grad;
    dfdx_B.grad = (0.0) * x_A.grad + (2.0) * x_B.grad;

    this->volumetrize_free_energy();
    this->nondimensionalize_free_energy();
  }
};

//
template <int dim, int degree>
inline SystemContainer<dim, degree>::SystemContainer(
  const IsothermalSystem         &_isoSys,
  const userInputParameters<dim> &_userInputs)
  : isoSys(_isoSys)
  , userInputs(_userInputs)
{
  // For all phase names
  phase_fields.insert({"SOLID", new SOLID<dim, degree>(isoSys, "SOLID", phase_fields)});
  phase_fields.insert(
    {"LIQUID", new LIQUID<dim, degree>(isoSys, "LIQUID", phase_fields)});
}
template class SystemContainer<2, 1>;
template class SystemContainer<2, 2>;
template class SystemContainer<2, 3>;
template class SystemContainer<2, 4>;
template class SystemContainer<2, 5>;
template class SystemContainer<2, 6>;
template class SystemContainer<3, 1>;
template class SystemContainer<3, 2>;
template class SystemContainer<3, 3>;
template class SystemContainer<3, 4>;
template class SystemContainer<3, 5>;
template class SystemContainer<3, 6>;