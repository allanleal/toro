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

#include "toroThermalBase.h"

template <>
InputParameters
validParams<toroThermalBase>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("Base class for calculating the thermal properties.");
  params.addRequiredCoupledVar("temperature", "The temperature variable.");
  params.addCoupledVar("porosity", "The porosity auxiliary variable.");
  params.addParam<Real>("radiogenic_heat_production", 0.0, "The radiogenic heat production value.");
  return params;
}

toroThermalBase::toroThermalBase(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _temp_dot(_fe_problem.isTransient() ? coupledDot("temperature") : _zero),
    _porosity(isCoupled("porosity") ? coupledValue("porosity") : _zero),
    _heat_source(getParam<Real>("radiogenic_heat_production")),
    _rho_f(getDefaultMaterialProperty<Real>("fluid_density")),
    _rho_s(getDefaultMaterialProperty<Real>("solid_density")),
    _thermal_diff(declareProperty<Real>("thermal_diffusivity")),
    _rhoC_b(declareProperty<Real>("bulk_specific_heat")),
    _rhoC_f(declareProperty<Real>("fluid_specific_heat")),
    _thermal_exp(declareProperty<Real>("thermal_expansion_coefficient")),
    _thermal_strain(declareProperty<RankTwoTensor>("thermal_strain")),
    _thermal_strain_old(getMaterialPropertyOld<RankTwoTensor>("thermal_strain")),
    _dthermal_strain_dtemp(declareProperty<RankTwoTensor>("dthermal_strain_dtemp")),
    _radiogenic_heat(declareProperty<Real>("radiogenic_heat_production")),
    _c_f(_fe_problem.getMaxQps()),
    _c_s(_fe_problem.getMaxQps()),
    _lambda_f(_fe_problem.getMaxQps()),
    _lambda_s(_fe_problem.getMaxQps()),
    _beta_f(_fe_problem.getMaxQps()),
    _beta_s(_fe_problem.getMaxQps())
{
}

void
toroThermalBase::initQpStatefulProperties()
{
  _thermal_strain[_qp].zero();
}

void
toroThermalBase::computeQpProperties()
{
  computeQpSpecificHeat();
  computeQpThermalDiff();
  computeQpThermalExpansion();
  computeQpThermalStrain();
  computeQpThermalSource();
}

void
toroThermalBase::computeQpSpecificHeat()
{
  computeQpHeatCap();

  _rhoC_f[_qp] = _rho_f[_qp] * _c_f[_qp];
  _rhoC_b[_qp] = computeMixtureProperty(_rhoC_f[_qp], _rho_s[_qp] * _c_s[_qp]);
}

void
toroThermalBase::computeQpThermalDiff()
{
  computeQpThermalCond();

  _thermal_diff[_qp] = computeMixtureProperty(_lambda_f[_qp], _lambda_s[_qp]);
  if (_rhoC_b[_qp] != 0.0)
    _thermal_diff[_qp] /= _rhoC_b[_qp];
}

void
toroThermalBase::computeQpThermalExpansion()
{
  computeQpThermalExp();

  _thermal_exp[_qp] = computeMixtureProperty(_beta_f[_qp], _beta_s[_qp]);
}

void
toroThermalBase::computeQpThermalStrain()
{
  RankTwoTensor thermal_strain_incr = RankTwoTensor();
  thermal_strain_incr.addIa(_thermal_exp[_qp] * _temp_dot[_qp] / 3.0);

  _thermal_strain[_qp] = _thermal_strain_old[_qp] + thermal_strain_incr * _dt;
  _dthermal_strain_dtemp[_qp].zero();
  _dthermal_strain_dtemp[_qp].addIa(_thermal_exp[_qp] / 3.0);
  // + dthermal_exp_dtemp * _temp_dot / 3
}

void
toroThermalBase::computeQpThermalSource()
{
  // In this material, we compute only radiogenic heat production
  // In toroElasticRheology, we compute shear heating and adiabatic heating
  _radiogenic_heat[_qp] = _heat_source;
}

Real
toroThermalBase::computeMixtureProperty(const Real fluid_prop, const Real solid_prop)
{
  return _porosity[_qp] * fluid_prop + (1.0 - _porosity[_qp]) * solid_prop;
}
