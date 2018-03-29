/******************************************************************************/
/*                       toro, a MOOSE-based application                      */
/*                                                                            */
/*          Copyright (C) 2017 by Antoine B. Jacquey and Mauro Cacace         */
/*             GFZ Potsdam, German Research Centre for Geosciences            */
/*                                                                            */
/*    This program is free software: you can redistribute it and/or modify    */
/*    it under the terms of the GNU General Public License as published by    */
/*      the Free Software Foundation, either version 3 of the License, or     */
/*                     (at your option) any later version.                    */
/*                                                                            */
/*       This program is distributed in the hope that it will be useful,      */
/*       but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       */
/*                GNU General Public License for more details.                */
/*                                                                            */
/*      You should have received a copy of the GNU General Public License     */
/*    along with this program. If not, see <http://www.gnu.org/licenses/>     */
/******************************************************************************/
#include "toroAdvectionImplicit.h"
#include "MooseMesh.h"

registerMooseObject("ToroApp", toroAdvectionImplicit);

template <>
InputParameters
validParams<toroAdvectionImplicit>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription(
      "Advection term corrected by an artificial entropy viscosity stabilization term after "
      "Guermond et al. 2011. "
      "based on an impliciti discretization of the advection term -  relative velocities for ALE "
      "implementation.");
  params.addRequiredCoupledVar("velocities", "The string of velocities that advect the problem.");
  params.addCoupledVar("mesh_velocities",
                       "The string of velocities that describe the motion of the mesh.");
  params.addParam<bool>(
      "use_displaced_mesh", true, "Set the displaced mesh flag to true by default.");
  params.addParam<bool>("entropy_viscosity_stabilization",
                        false,
                        "Whether to perform entropy viscosity stabilization.");
  params.addParam<bool>(
      "is_temperature", true, "Whether we are advective the temperature field or not.");
  params.addParam<Real>("beta_stabilization", 0.026, "The beta local stabilization parameter.");
  params.addParam<Real>("cr_stabilization", 0.5, "The cr local stabilization parameter.");
  params.addParam<Real>("gamma_stabilization", 0.0, "The gamma local stabilization parameter.");
  params.addParam<Real>("model_length", 1.0, "The model scaling length.");
  params.addParam<MooseEnum>("ele_diameter",
                             toroAdvectionImplicit::elementDiameter() = "min",
                             "The diameter of a single cell.");
  params.addRequiredParam<PostprocessorName>(
      "pp_max_vel", "The postprocessor to retrieve the maximum velocity on the whole domain.");
  params.addRequiredParam<PostprocessorName>(
      "pp_max_var",
      "The postprocessor to retrieve the maximum advected variable on the whole domain.");
  params.addRequiredParam<PostprocessorName>(
      "pp_min_var",
      "The postprocessor to retrieve the minimum advected variable on the whole domain.");
  params.addRequiredParam<PostprocessorName>(
      "pp_avg_var",
      "The postprocessor to retrieve the average advected variable on the whole domain.");
  params.addRequiredParam<PostprocessorName>(
      "pp_max_entropy", "The postprocessor to retrieve the maximum entropy on the whole domain.");
  params.addRequiredParam<PostprocessorName>(
      "pp_min_entropy", "The postprocessor to retrieve the minimum entropy on the whole domain.");
  params.addRequiredParam<PostprocessorName>(
      "pp_avg_entropy", "The postprocessor to retrieve the average entropy on the whole domain.");
  params.addParam<Real>(
      "coeff_shear_heating", 0.0, "The coefficient in front of the shear heating generation.");
  return params;
}

toroAdvectionImplicit::toroAdvectionImplicit(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _entropy_viscosity_stabilization(getParam<bool>("entropy_viscosity_stabilization")),
    _is_temperature(getParam<bool>("is_temperature")),
    _beta_stabilization(getParam<Real>("beta_stabilization")),
    _cr_stabilization(getParam<Real>("cr_stabilization")),
    _gamma_stabilization(getParam<Real>("gamma_stabilization")),
    _model_length(getParam<Real>("model_length")),
    _ele_diameter(getParam<MooseEnum>("ele_diameter")),
    _n_vel(coupledComponents("velocities")),
    _vel_var(_n_vel),
    _vel(3),
    _vel_old(3),
    _vel_older(3),
    _has_mesh_velocity(isCoupled("mesh_velocities")),
    _n_vel_mesh(coupledComponents("mesh_velocities")),
    _vel_mesh(3),
    _vel_old_mesh(3),
    _vel_older_mesh(3),
    _value_old(_fe_problem.isTransient() ? valueOld() : _zero),
    _gradient_old(_fe_problem.isTransient() ? gradientOld() : _grad_zero),
    _second_old(_fe_problem.isTransient() ? secondOld() : _second_zero),
    _value_older(_fe_problem.isTransient() ? valueOlder() : _zero),
    _gradient_older(_fe_problem.isTransient() ? gradientOlder() : _grad_zero),
    _second_older(_fe_problem.isTransient() ? secondOlder() : _second_zero),
    _thermal_diff(getDefaultMaterialProperty<Real>("thermal_diffusivity")),
    _coeff_Hs(getParam<Real>("coeff_shear_heating")),
    _rhoC_b(getDefaultMaterialProperty<Real>("bulk_specific_heat")),
    _radiogenic_heat(getDefaultMaterialProperty<Real>("radiogenic_heat_production")),
    _inelastic_heat(getDefaultMaterialProperty<Real>("inelastic_heat")),
    _adiabatic_heat(getDefaultMaterialProperty<Real>("adiabatic_heat")),
    _pp_max_vel(getPostprocessorValue("pp_max_vel")),
    _pp_max_var(getPostprocessorValue("pp_max_var")),
    _pp_min_var(getPostprocessorValue("pp_min_var")),
    _pp_avg_var(getPostprocessorValue("pp_avg_var")),
    _pp_max_entropy(getPostprocessorValue("pp_max_entropy")),
    _pp_min_entropy(getPostprocessorValue("pp_min_entropy")),
    _pp_avg_entropy(getPostprocessorValue("pp_avg_entropy")),
    _residual(_fe_problem.getMaxQps())
{
  if (!_fe_problem.isTransient())
    mooseError("@ toroAdvectionImplicit: it does not make any sense to run this kernel on a "
               "steady state.");
  if ((_n_vel != _mesh.dimension()) || (_has_mesh_velocity && _n_vel_mesh != _mesh.dimension()))
    mooseError("@ toroAdvectionImplicit: The number of supplied velocities does not match the "
               "problem dimension.");
  if (_has_mesh_velocity && (_n_vel != _n_vel_mesh))
    mooseError("@ toroAdvectionImplicit: the number of components of the supplied velovities does "
               "not match.");

  for (unsigned i = 0; i < _n_vel; ++i)
  {
    _vel_var[i] = coupled("velocities", i);
    _vel[i] = &coupledValue("velocities", i);
    _vel_old[i] = &coupledValueOld("velocities", i);
    _vel_older[i] = &coupledValueOlder("velocities", i);
  }
  for (unsigned i = _n_vel; i < 3; ++i)
  {
    _vel[i] = &_zero;
    _vel_old[i] = &_zero;
    _vel_older[i] = &_zero;
  }

  if (_has_mesh_velocity)
  {
    for (unsigned i = 0; i < _n_vel; ++i)
    {
      _vel_mesh[i] = &coupledValue("mesh_velocities", i);
      _vel_old_mesh[i] = &coupledValueOld("mesh_velocities", i);
      _vel_older_mesh[i] = &coupledValueOlder("mesh_velocities", i);
    }
    for (unsigned i = _n_vel; i < 3; ++i)
    {
      _vel_mesh[i] = &_zero;
      _vel_old_mesh[i] = &_zero;
      _vel_older_mesh[i] = &_zero;
    }
  }
  else
    for (unsigned i = 0; i < 3; ++i)
    {
      _vel_mesh[i] = &_zero;
      _vel_old_mesh[i] = &_zero;
      _vel_older_mesh[i] = &_zero;
    }

  if (_cr_stabilization > 0.0)
  {
    _strain_rate_old = &getMaterialPropertyOld<RankTwoTensor>("total_strain_rate");
    _strain_rate_older = &getMaterialPropertyOlder<RankTwoTensor>("total_strain_rate");
  }
}

MooseEnum
toroAdvectionImplicit::elementDiameter()
{
  return MooseEnum("min=1 max=2 average=3");
}

Real
toroAdvectionImplicit::computeElementDiameter()
{
  Real diameter = 0.0;

  if (_current_elem->dim() == 1)
    diameter += _current_elem->volume();
  else
  {
    switch (_ele_diameter)
    {
      case 1: // min
        diameter += _current_elem->hmin();
        break;
      case 2: // max
        diameter += _current_elem->hmax();
        break;
      case 3: // average
        diameter += 0.5 * (_current_elem->hmin() + _current_elem->hmax());
        break;
    }
  }
  return diameter;
}

void
toroAdvectionImplicit::computeArtificialViscosity()
{
  Real diameter = computeElementDiameter();

  computeEntropyResidual();

  Real max_residual = 0.0;
  Real max_velocity = 0.0;

  Real max_diffusivity = _is_temperature ? 0.0 : 1.0;
  Real max_rho_cp = _is_temperature ? 0.0 : 1.0;

  for (unsigned int qp = 0; qp < _q_point.size(); ++qp)
  {
    // should one make use of the relative velocities?!
    Real u0 = 0.5 * ((*_vel_old[0])[qp] + (*_vel_older[0])[qp]);
    Real u1 = 0.5 * ((*_vel_old[1])[qp] + (*_vel_older[1])[qp]);
    Real u2 = 0.5 * ((*_vel_old[2])[qp] + (*_vel_older[2])[qp]);
    RealVectorValue u(u0, u1, u2);

    Real strain_rate = 0.0;
    if (_cr_stabilization > 0.0)
      strain_rate =
          (_t_step <= 1)
              ? (*_strain_rate_old)[qp].trace()
              : 0.5 * ((*_strain_rate_old)[qp].trace() + (*_strain_rate_older)[qp].trace());

    Real field = (_value_old[qp] + _value_older[qp]) * 0.5;
    _residual[qp] *= std::abs(field - _pp_avg_var);

    max_residual = std::max(_residual[qp], max_residual);
    max_velocity = std::max(u.norm() + _gamma_stabilization * strain_rate * diameter, max_velocity);

    if (_is_temperature)
    {
      max_diffusivity = std::max(_thermal_diff[qp], max_diffusivity);
      max_rho_cp = std::max(_rhoC_b[qp], max_rho_cp);
    }
  }

  // If velocity is null, assume a sensible value to get back
  // an artificial diffusion
  // here we do: velocity ~ _diffusivity / length_scale
  if (std::abs(_pp_max_vel) < 1e-50)
  {
    _artificial_viscosity =
        _beta_stabilization * (max_diffusivity / max_rho_cp) / _model_length * diameter;
    return;
  }

  Real max_viscosity = _beta_stabilization * max_rho_cp * max_velocity * diameter;
  Real entropy_variation =
      std::max(_pp_max_entropy - _pp_avg_entropy, _pp_avg_entropy - _pp_min_entropy);

  if (_t_step <= 1 || std::abs(entropy_variation) < 1e-50)
  {
    _artificial_viscosity = max_viscosity;
    return;
  }

  Real entropy_viscosity =
      _cr_stabilization * diameter * diameter * max_residual / entropy_variation;

  _artificial_viscosity = std::min(max_viscosity, entropy_viscosity);
}

void
toroAdvectionImplicit::computeEntropyResidual()
{
  for (unsigned int qp = 0; qp < _q_point.size(); ++qp)
  {
    // should one make use of the relative velocities?!
    Real u0 = 0.5 * ((*_vel_old[0])[qp] + (*_vel_older[0])[qp]);
    Real u1 = 0.5 * ((*_vel_old[1])[qp] + (*_vel_older[1])[qp]);
    Real u2 = 0.5 * ((*_vel_old[2])[qp] + (*_vel_older[2])[qp]);
    RealVectorValue u(u0, u1, u2);

    Real dvar_dt = (_dt_old == 0.0) ? 0.0 : (_value_old[qp] - _value_older[qp]) / _dt_old;
    Real u_grad_var = u * (_gradient_old[qp] + _gradient_older[qp]) * 0.5;

    if (_is_temperature)
    {
      Real laplace = 0.5 * (_second_old[qp].tr() + _second_older[qp].tr());
      Real kappa_laplace_var = _thermal_diff[qp] * laplace;
      Real Hr = _radiogenic_heat[qp];
      Real Hs = _coeff_Hs * _inelastic_heat[qp];
      Real Ha = _adiabatic_heat[qp];
      Real heat_sources = Hr + Hs + Ha;
      if (_rhoC_b[qp] != 0.0)
        heat_sources /= _rhoC_b[qp];

      _residual[qp] = std::abs(dvar_dt + u_grad_var - kappa_laplace_var - heat_sources);
      return;
    }
    _residual[qp] = std::abs(dvar_dt + u_grad_var);
  }
}

/****************************************************************************** /
/*                                  RESIDUALS                                 */
/******************************************************************************/

void
toroAdvectionImplicit::precalculateResidual()
{
  if (_entropy_viscosity_stabilization)
    computeArtificialViscosity();
}

Real
toroAdvectionImplicit::computeQpResidual()
{
  Real u0 = (*_vel[0])[_qp] - (*_vel_mesh[0])[_qp];
  Real u1 = (*_vel[1])[_qp] - (*_vel_mesh[1])[_qp];
  Real u2 = (*_vel[2])[_qp] - (*_vel_mesh[2])[_qp];
  RealVectorValue u(u0, u1, u2);

  Real res = u * _grad_u[_qp] * _test[_i][_qp];

  if (_entropy_viscosity_stabilization)
    res += _artificial_viscosity * _grad_test[_i][_qp] * _grad_u[_qp];

  return res;
}

/******************************************************************************/
/*                                  JACOBIAN                                  */
/******************************************************************************/

Real
toroAdvectionImplicit::computeQpJacobian()
{
  Real u0 = (*_vel[0])[_qp] - (*_vel_mesh[0])[_qp];
  Real u1 = (*_vel[1])[_qp] - (*_vel_mesh[1])[_qp];
  Real u2 = (*_vel[2])[_qp] - (*_vel_mesh[2])[_qp];
  RealVectorValue u(u0, u1, u2);

  Real jac = u * _grad_phi[_j][_qp] * _test[_i][_qp];

  if (_entropy_viscosity_stabilization)
    jac += _artificial_viscosity * _grad_test[_i][_qp] * _grad_phi[_j][_qp];

  return jac;
}

/******************************************************************************/
/*                              OFF-DIAG JACOBIAN                             */
/******************************************************************************/

Real
toroAdvectionImplicit::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real jac = 0.0;
  for (unsigned int coupled_component = 0; coupled_component < _n_vel; ++coupled_component)
    if (jvar == _vel_var[coupled_component])
      jac += _phi[_j][_qp] * _grad_u[_qp](coupled_component) * _test[_i][_qp];

  return jac;
}
