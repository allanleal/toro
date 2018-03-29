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

#include "toroSolidMomentum.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "Assembly.h"
#include "SystemBase.h"
#include "libmesh/quadrature.h"

registerMooseObject("ToroApp", toroSolidMomentum);

template <>
InputParameters
validParams<toroSolidMomentum>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Solid momentum kernel.");
  params.addRequiredCoupledVar(
      "displacements", "The string of displacements variables suitable for the problem statement.");
  params.addCoupledVar("temperature", "The temperature variable.");
  params.addCoupledVar("fluid_pressure", "The fluid pressure variable.");
  params.addCoupledVar("damage", "The damage variable.");
  params.set<bool>("use_displaced_mesh") = false;
  params.addParam<bool>(
      "use_finite_deform_jacobian", false, "Jacobian for corotational finite strain.");
  params.addParam<bool>("volumetric_locking_correction",
                        false,
                        "Set to false to turn off volumetric locking correction.");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction "
                                        "the variable this kernel acts in (0 for x, "
                                        "1 for y, 2 for z).");
  return params;
}

toroSolidMomentum::toroSolidMomentum(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _ndisp(coupledComponents("displacements")),
    _coupled_temp(isCoupled("temperature")),
    _coupled_pf(isCoupled("fluid_pressure")),
    _coupled_dam(isCoupled("damage")),
    _pf(_coupled_pf ? coupledValue("fluid_pressure") : _zero),
    _use_finite_deform_jacobian(getParam<bool>("use_finite_deform_jacobian")),
    _volumetric_locking_correction(getParam<bool>("volumetric_locking_correction")),
    _component(getParam<unsigned int>("component")),
    _biot(getDefaultMaterialProperty<Real>("biot_coefficient")),
    _gravity(getDefaultMaterialProperty<RealVectorValue>("gravity_vector")),
    _rho_b(getDefaultMaterialProperty<Real>("bulk_density")),
    _stress(getMaterialProperty<RankTwoTensor>("stress")),
    _tangent_modulus(getMaterialProperty<RankFourTensor>("tangent_modulus")),
    _dthermal_strain_dtemp(getDefaultMaterialProperty<RankTwoTensor>("dthermal_strain_dtemp")),
    _dstress_ddamage(getDefaultMaterialProperty<RankTwoTensor>("dstress_ddamage")),
    _disp_var(_ndisp),
    _temp_var(_coupled_temp ? coupled("temperature") : 0),
    _pf_var(_coupled_pf ? coupled("fluid_pressure") : 0),
    _dam_var(_coupled_dam ? coupled("damage") : 0),
    _assembly_undisplaced(_fe_problem.assembly(_tid)),
    _var_undisplaced(
        _fe_problem.getStandardVariable(_tid, parameters.get<NonlinearVariableName>("variable"))),
    _grad_phi_undisplaced(_assembly_undisplaced.gradPhi()),
    _grad_test_undisplaced(_var_undisplaced.gradPhi()),
    _avg_grad_test(_test.size(), std::vector<Real>(3, 0.0)),
    _avg_grad_phi(_phi.size(), std::vector<Real>(3, 0.0))
{
  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    mooseError("The number of displacement variables supplied must match the mesh dimension.");

  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);

  if (_use_finite_deform_jacobian)
  {
    _deformation_gradient = &getMaterialProperty<RankTwoTensor>("deformation_gradient");
    _deformation_gradient_old = &getMaterialPropertyOld<RankTwoTensor>("deformation_gradient");
    _rotation_increment = &getMaterialProperty<RankTwoTensor>("rotation_increment");
  }

  // Error if volumetic locking correction is turned on for 1D problems
  if (_ndisp == 1 && _volumetric_locking_correction)
    mooseError("Volumetric locking correction should be set to false for 1-D problems.");
}

/******************************************************************************/
/*                                  RESIDUALS                                 */
/******************************************************************************/

void
toroSolidMomentum::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  _local_re.resize(re.size());
  _local_re.zero();

  if (_volumetric_locking_correction)
    computeAverageGradientTest();

  precalculateResidual();
  for (_i = 0; _i < _test.size(); ++_i)
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      _local_re(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidual();

  re += _local_re;

  if (_has_save_in)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (const auto & var : _save_in)
      var->sys().solution().add_vector(_local_re, var->dofIndices());
  }
}

Real
toroSolidMomentum::computeQpResidual()
{
  RealVectorValue grav_term = -_rho_b[_qp] * _gravity[_qp];
  RealVectorValue poro_stress_row = _stress[_qp].row(_component);
  poro_stress_row(_component) -= _biot[_qp] * _pf[_qp];

  Real residual = poro_stress_row * _grad_test[_i][_qp] + grav_term(_component) * _test[_i][_qp];

  return residual;
}

/******************************************************************************/
/*                                  JACOBIAN                                  */
/******************************************************************************/

void
toroSolidMomentum::computeJacobian()
{
  if (_use_finite_deform_jacobian)
  {
    _finite_deform_tangent_modulus.resize(_qrule->n_points());

    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      computeFiniteDeformModulus();

    _fe_problem.prepareShapes(_var.number(), _tid);
    Kernel::computeJacobian();
  }
  else
  {
    if (_volumetric_locking_correction)
    {
      computeAverageGradientTest();
      computeAverageGradientPhi();
    }
    Kernel::computeJacobian();
  }
}

Real
toroSolidMomentum::computeQpJacobian()
{
  if (_use_finite_deform_jacobian)
    return elasticJacobian(_finite_deform_tangent_modulus[_qp],
                           _component,
                           _component,
                           _grad_test[_i][_qp],
                           _grad_phi_undisplaced[_j][_qp]);

  Real sum_C3x3 = _tangent_modulus[_qp].sum3x3();
  RealGradient sum_C3x1 = _tangent_modulus[_qp].sum3x1();

  Real jacobian = 0.0;
  jacobian += elasticJacobian(
      _tangent_modulus[_qp], _component, _component, _grad_test[_i][_qp], _grad_phi[_j][_qp]);

  if (_volumetric_locking_correction)
  {
    // jacobian = Bbar^T_i * C * Bbar_j where Bbar = B + Bvol
    // jacobian = B^T_i * C * B_j + Bvol^T_i * C * Bvol_j +  Bvol^T_i * C * B_j + B^T_i * C * Bvol_j

    // Bvol^T_i * C * Bvol_j
    jacobian += sum_C3x3 * (_avg_grad_test[_i][_component] - _grad_test[_i][_qp](_component)) *
                (_avg_grad_phi[_j][_component] - _grad_phi[_j][_qp](_component)) / 9.0;

    // B^T_i * C * Bvol_j
    jacobian += sum_C3x1(_component) * _grad_test[_i][_qp](_component) *
                (_avg_grad_phi[_j][_component] - _grad_phi[_j][_qp](_component)) / 3.0;

    // Bvol^T_i * C * B_j
    RankTwoTensor phi;
    if (_component == 0)
    {
      phi(0, 0) = _grad_phi[_j][_qp](0);
      phi(0, 1) = phi(1, 0) = _grad_phi[_j][_qp](1);
      phi(0, 2) = phi(2, 0) = _grad_phi[_j][_qp](2);
    }
    else if (_component == 1)
    {
      phi(1, 1) = _grad_phi[_j][_qp](1);
      phi(0, 1) = phi(1, 0) = _grad_phi[_j][_qp](0);
      phi(1, 2) = phi(2, 1) = _grad_phi[_j][_qp](2);
    }
    else if (_component == 2)
    {
      phi(2, 2) = _grad_phi[_j][_qp](2);
      phi(0, 2) = phi(2, 0) = _grad_phi[_j][_qp](0);
      phi(1, 2) = phi(2, 1) = _grad_phi[_j][_qp](1);
    }

    jacobian += (_tangent_modulus[_qp] * phi).trace() *
                (_avg_grad_test[_i][_component] - _grad_test[_i][_qp](_component)) / 3.0;
  }
  return jacobian;
}

/******************************************************************************/
/*                              OFF-DIAG JACOBIAN                             */
/******************************************************************************/

void
toroSolidMomentum::computeOffDiagJacobian(unsigned int jvar)
{
  if (_use_finite_deform_jacobian)
  {
    _finite_deform_tangent_modulus.resize(_qrule->n_points());

    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      computeFiniteDeformModulus();

    _fe_problem.prepareShapes(jvar, _tid);
    Kernel::computeOffDiagJacobian(jvar);
  }
  else
  {
    if (_volumetric_locking_correction)
    {
      computeAverageGradientPhi();
      computeAverageGradientTest();
    }
    Kernel::computeOffDiagJacobian(jvar);
  }
}

Real
toroSolidMomentum::computeQpOffDiagJacobian(unsigned int jvar)
{
  // off-diagonal Jacobian with respect to a coupled displacement component
  for (unsigned int coupled_component = 0; coupled_component < _ndisp; ++coupled_component)
    if (jvar == _disp_var[coupled_component])
    {
      if (_use_finite_deform_jacobian)
        return elasticJacobian(_finite_deform_tangent_modulus[_qp],
                               _component,
                               coupled_component,
                               _grad_test[_i][_qp],
                               _grad_phi_undisplaced[_j][_qp]);

      const Real sum_C3x3 = _tangent_modulus[_qp].sum3x3();
      const RealGradient sum_C3x1 = _tangent_modulus[_qp].sum3x1();

      Real jacobian = 0.0;
      jacobian += elasticJacobian(_tangent_modulus[_qp],
                                  _component,
                                  coupled_component,
                                  _grad_test[_i][_qp],
                                  _grad_phi[_j][_qp]);

      if (_volumetric_locking_correction)
      {
        // jacobian = Bbar^T_i * C * Bbar_j where Bbar = B + Bvol
        // jacobian = B^T_i * C * B_j + Bvol^T_i * C * Bvol_j +  Bvol^T_i * C * B_j + B^T_i * C *
        // Bvol_j

        // Bvol^T_i * C * Bvol_j
        jacobian += sum_C3x3 * (_avg_grad_test[_i][_component] - _grad_test[_i][_qp](_component)) *
                    (_avg_grad_phi[_j][coupled_component] - _grad_phi[_j][_qp](coupled_component)) /
                    9.0;

        // B^T_i * C * Bvol_j
        jacobian += sum_C3x1(_component) * _grad_test[_i][_qp](_component) *
                    (_avg_grad_phi[_j][coupled_component] - _grad_phi[_j][_qp](coupled_component)) /
                    3.0;

        // Bvol^T_i * C * B_i
        RankTwoTensor phi;
        for (unsigned int i = 0; i < 3; ++i)
          phi(coupled_component, i) = _grad_phi[_j][_qp](i);

        jacobian += (_tangent_modulus[_qp] * phi).trace() *
                    (_avg_grad_test[_i][_component] - _grad_test[_i][_qp](_component)) / 3.0;
      }

      return jacobian;
    }

  if (_coupled_temp && jvar == _temp_var)
    return -(_tangent_modulus[_qp] * _dthermal_strain_dtemp[_qp] *
             _grad_test[_i][_qp])(_component)*_phi[_j][_qp];

  if (_coupled_pf && jvar == _pf_var)
    return -_biot[_qp] * _phi[_j][_qp] * _grad_test[_i][_qp](_component);

  if (_coupled_dam && jvar == _dam_var)
    return -(_dstress_ddamage[_qp] * _grad_test[_i][_qp])(_component)*_phi[_j][_qp];

  return 0.0;
}

Real
toroSolidMomentum::elasticJacobian(const RankFourTensor & jacobian_r4t,
                                   unsigned int i,
                                   unsigned int k,
                                   const RealGradient & grad_test,
                                   const RealGradient & grad_phi)
{
  Real sum = 0.0;
  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    for (unsigned int l = 0; l < LIBMESH_DIM; ++l)
      sum += jacobian_r4t(i, j, k, l) * grad_phi(l) * grad_test(j);
  return sum;
}

void
toroSolidMomentum::computeFiniteDeformModulus()
{
  const RankTwoTensor I(RankTwoTensor::initIdentity);
  const RankFourTensor II_ijkl = I.mixedProductIkJl(I);

  // Bring back to unrotated config
  const RankTwoTensor unrotated_stress =
      (*_rotation_increment)[_qp].transpose() * _stress[_qp] * (*_rotation_increment)[_qp];

  // Incremental deformation gradient Fhat
  const RankTwoTensor Fhat =
      (*_deformation_gradient)[_qp] * (*_deformation_gradient_old)[_qp].inverse();
  const RankTwoTensor Fhatinv = Fhat.inverse();

  const RankTwoTensor rot_times_stress = (*_rotation_increment)[_qp] * unrotated_stress;
  const RankFourTensor dstress_drot =
      I.mixedProductIkJl(rot_times_stress) + I.mixedProductJkIl(rot_times_stress);
  const RankFourTensor rot_rank_four =
      (*_rotation_increment)[_qp].mixedProductIkJl((*_rotation_increment)[_qp]);
  const RankFourTensor drot_dUhatinv = Fhat.mixedProductIkJl(I);

  const RankTwoTensor A = I - Fhatinv;

  // Ctilde = Chat^-1 - I
  const RankTwoTensor Ctilde = A * A.transpose() - A - A.transpose();
  const RankFourTensor dCtilde_dFhatinv =
      -I.mixedProductIkJl(A) - I.mixedProductJkIl(A) + II_ijkl + I.mixedProductJkIl(I);

  // Second order approximation of Uhat - consistent with strain increment definition
  // const RankTwoTensor Uhat = I - 0.5 * Ctilde - 3.0/8.0 * Ctilde * Ctilde;

  RankFourTensor dUhatinv_dCtilde =
      0.5 * II_ijkl - 1.0 / 8.0 * (I.mixedProductIkJl(Ctilde) + Ctilde.mixedProductIkJl(I));
  RankFourTensor drot_dFhatinv = drot_dUhatinv * dUhatinv_dCtilde * dCtilde_dFhatinv;

  drot_dFhatinv -= Fhat.mixedProductIkJl((*_rotation_increment)[_qp].transpose());
  _finite_deform_tangent_modulus[_qp] = dstress_drot * drot_dFhatinv;

  const RankFourTensor dstrain_increment_dCtilde =
      -0.5 * II_ijkl + 0.25 * (I.mixedProductIkJl(Ctilde) + Ctilde.mixedProductIkJl(I));
  _finite_deform_tangent_modulus[_qp] +=
      rot_rank_four * _tangent_modulus[_qp] * dstrain_increment_dCtilde * dCtilde_dFhatinv;
  _finite_deform_tangent_modulus[_qp] += Fhat.mixedProductJkIl(_stress[_qp]);

  const RankFourTensor dFhat_dFhatinv = -Fhat.mixedProductIkJl(Fhat.transpose());
  const RankTwoTensor dJ_dFhatinv = dFhat_dFhatinv.innerProductTranspose(Fhat.ddet());

  // Component from Jacobian derivative
  _finite_deform_tangent_modulus[_qp] += _stress[_qp].outerProduct(dJ_dFhatinv);

  // Derivative of Fhatinv w.r.t. undisplaced coordinates
  const RankTwoTensor Finv = (*_deformation_gradient)[_qp].inverse();
  const RankFourTensor dFhatinv_dGradu = -Fhatinv.mixedProductIkJl(Finv.transpose());
  _finite_deform_tangent_modulus[_qp] = _finite_deform_tangent_modulus[_qp] * dFhatinv_dGradu;
}

void
toroSolidMomentum::computeAverageGradientTest()
{
  // Calculate volume averaged value of shape function derivative
  _avg_grad_test.resize(_test.size());
  for (_i = 0; _i < _test.size(); ++_i)
  {
    _avg_grad_test[_i].resize(3);
    _avg_grad_test[_i][_component] = 0.0;
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      _avg_grad_test[_i][_component] += _grad_test[_i][_qp](_component) * _JxW[_qp] * _coord[_qp];

    _avg_grad_test[_i][_component] /= _current_elem_volume;
  }
}

void
toroSolidMomentum::computeAverageGradientPhi()
{
  // Calculate volume average derivatives for phi
  _avg_grad_phi.resize(_phi.size());
  for (_i = 0; _i < _phi.size(); ++_i)
  {
    _avg_grad_phi[_i].resize(3);
    for (unsigned int component = 0; component < _mesh.dimension(); ++component)
    {
      _avg_grad_phi[_i][component] = 0.0;
      for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
        _avg_grad_phi[_i][component] += _grad_phi[_i][_qp](component) * _JxW[_qp] * _coord[_qp];

      _avg_grad_phi[_i][component] /= _current_elem_volume;
    }
  }
}
