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
#ifndef toroADVECTIONIMPLICIT_H
#define toroADVECTIONIMPLICIT_H

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"
#include "RankTwoTensor.h"

class toroAdvectionImplicit;

template <>
InputParameters validParams<toroAdvectionImplicit>();

class toroAdvectionImplicit : public DerivativeMaterialInterface<Kernel>
{
public:
  toroAdvectionImplicit(const InputParameters & parameters);
  static MooseEnum elementDiameter();

protected:
  virtual void precalculateResidual() override;
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int) override;

  virtual Real computeElementDiameter();
  virtual void computeArtificialViscosity();
  virtual void computeEntropyResidual();

  bool _entropy_viscosity_stabilization;
  bool _is_temperature;

  Real _beta_stabilization;
  Real _cr_stabilization;
  Real _gamma_stabilization;
  Real _model_length;

  MooseEnum _ele_diameter;

  // Eulerian velocities
  unsigned int _n_vel;
  std::vector<unsigned> _vel_var;
  std::vector<const VariableValue *> _vel;
  std::vector<const VariableValue *> _vel_old;
  std::vector<const VariableValue *> _vel_older;

  // mesh velocities - only set if the mesh is displaced by a different velocity than the solid one
  // or if running a porous problem
  bool _has_mesh_velocity;
  unsigned int _n_vel_mesh;
  std::vector<const VariableValue *> _vel_mesh;
  std::vector<const VariableValue *> _vel_old_mesh;
  std::vector<const VariableValue *> _vel_older_mesh;

  //  advected quantity - old and older values
  const VariableValue & _value_old;
  const VariableGradient & _gradient_old;
  const VariableSecond & _second_old;
  const VariableValue & _value_older;
  const VariableGradient & _gradient_older;
  const VariableSecond & _second_older;

  // diffusion residual - if advected quantity is temperature
  const MaterialProperty<Real> & _thermal_diff;
  // source term residual - if advected quantity is temperature
  Real _coeff_Hs;
  const MaterialProperty<Real> & _rhoC_b;
  const MaterialProperty<Real> & _radiogenic_heat;
  const MaterialProperty<Real> & _inelastic_heat;
  const MaterialProperty<Real> & _adiabatic_heat;

  // strain rate - still to be tested
  const MaterialProperty<RankTwoTensor> * _strain_rate_old;
  const MaterialProperty<RankTwoTensor> * _strain_rate_older;

  // postprocessor to get the global maxima and minima of all quantities
  const PostprocessorValue & _pp_max_vel;
  const PostprocessorValue & _pp_max_var;
  const PostprocessorValue & _pp_min_var;
  const PostprocessorValue & _pp_avg_var;
  const PostprocessorValue & _pp_max_entropy;
  const PostprocessorValue & _pp_min_entropy;
  const PostprocessorValue & _pp_avg_entropy;

private:
  Real _artificial_viscosity;
  std::vector<Real> _residual;
};

#endif // toroADVECTIONIMPLICIT_H
