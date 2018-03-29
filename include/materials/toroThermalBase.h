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

#ifndef toroTHERMALBASE_H
#define toroTHERMALBASE_H

#include "Material.h"
#include "RankTwoTensor.h"
#include "DerivativeMaterialInterface.h"

class toroThermalBase;

template <>
InputParameters validParams<toroThermalBase>();

class toroThermalBase : public DerivativeMaterialInterface<Material>
{
public:
  toroThermalBase(const InputParameters & parameters);
  virtual ~toroThermalBase() {}

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  virtual void computeQpSpecificHeat();
  virtual void computeQpThermalDiff();
  virtual void computeQpThermalExpansion();
  virtual void computeQpThermalStrain();
  virtual void computeQpThermalSource();
  virtual void computeQpHeatCap() = 0;
  virtual void computeQpThermalCond() = 0;
  virtual void computeQpThermalExp() = 0;
  virtual Real computeMixtureProperty(const Real fluid_prop, const Real solid_prop);

  const VariableValue & _temp_dot;
  const VariableValue & _porosity;
  Real _heat_source;

  const MaterialProperty<Real> & _rho_f;
  const MaterialProperty<Real> & _rho_s;
  MaterialProperty<Real> & _thermal_diff;
  MaterialProperty<Real> & _rhoC_b;
  MaterialProperty<Real> & _rhoC_f;
  MaterialProperty<Real> & _thermal_exp;
  MaterialProperty<RankTwoTensor> & _thermal_strain;
  const MaterialProperty<RankTwoTensor> & _thermal_strain_old;
  MaterialProperty<RankTwoTensor> & _dthermal_strain_dtemp;
  MaterialProperty<Real> & _radiogenic_heat;

  std::vector<Real> _c_f;
  std::vector<Real> _c_s;
  std::vector<Real> _lambda_f;
  std::vector<Real> _lambda_s;
  std::vector<Real> _beta_f;
  std::vector<Real> _beta_s;
};

#endif // toroTHERMALBASE_H
