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

#ifndef toroSOLIDMOMENTUM_H
#define toroSOLIDMOMENTUM_H

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"
#include "RankTwoTensor.h"

class toroSolidMomentum;

template <>
InputParameters validParams<toroSolidMomentum>();

class toroSolidMomentum : public DerivativeMaterialInterface<Kernel>
{
public:
  toroSolidMomentum(const InputParameters & parameters);

  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(unsigned int jvar) override;

protected:
  virtual void computeResidual() override;
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  virtual Real elasticJacobian(const RankFourTensor & jacobian_r4t,
                               unsigned int i,
                               unsigned int k,
                               const RealGradient & grad_test,
                               const RealGradient & grad_phi);
  virtual void computeFiniteDeformModulus();
  virtual void computeAverageGradientTest();
  virtual void computeAverageGradientPhi();

  unsigned int _ndisp;
  bool _coupled_temp;
  bool _coupled_pf;
  bool _coupled_dam;
  const VariableValue & _pf;
  bool _use_finite_deform_jacobian;
  bool _volumetric_locking_correction;
  const unsigned int _component;
  const MaterialProperty<Real> & _biot;
  const MaterialProperty<RealVectorValue> & _gravity;
  const MaterialProperty<Real> & _rho_b;
  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankFourTensor> & _tangent_modulus;
  const MaterialProperty<RankTwoTensor> & _dthermal_strain_dtemp;
  const MaterialProperty<RankTwoTensor> & _dstress_ddamage;
  std::vector<unsigned int> _disp_var;
  unsigned int _temp_var;
  unsigned int _pf_var;
  unsigned int _dam_var;
  std::vector<RankFourTensor> _finite_deform_tangent_modulus;
  Assembly & _assembly_undisplaced;
  MooseVariable & _var_undisplaced;
  const VariablePhiGradient & _grad_phi_undisplaced;
  const VariableTestGradient & _grad_test_undisplaced;
  const MaterialProperty<RankTwoTensor> * _deformation_gradient;
  const MaterialProperty<RankTwoTensor> * _deformation_gradient_old;
  const MaterialProperty<RankTwoTensor> * _rotation_increment;
  // Gradient of test function averaged over the element. Used in volumetric locking correction
  // calculation.
  std::vector<std::vector<Real>> _avg_grad_test;
  // Gradient of phi function averaged over the element. Used in volumetric locking correction
  // calculation.
  std::vector<std::vector<Real>> _avg_grad_phi;
};

#endif // toroSOLIDMOMENTUM_H
