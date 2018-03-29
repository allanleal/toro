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

#include "toroThermalConstant.h"

registerMooseObject("ToroApp", toroThermalConstant);

template <>
InputParameters
validParams<toroThermalConstant>()
{
  InputParameters params = validParams<toroThermalBase>();
  params.addClassDescription("Constant thermal properties.");
  params.addParam<Real>("fluid_thermal_conductivity", 0.0, "The fluid thermal conductivity.");
  params.addRequiredParam<Real>("solid_thermal_conductivity", "The solid thermal conductivity.");
  params.addParam<Real>("fluid_heat_capacity", 0.0, "The fluid heat capacity.");
  params.addParam<Real>("solid_heat_capacity", 0.0, "The solid heat capacity.");
  params.addParam<Real>(
      "fluid_thermal_expansion", 0.0, "The fluid volumetric thermal expansion coefficient.");
  params.addParam<Real>(
      "solid_thermal_expansion", 0.0, "The solid volumetric thermal expansion coefficient.");
  return params;
}

toroThermalConstant::toroThermalConstant(const InputParameters & parameters)
  : toroThermalBase(parameters),
    _fluid_thermal_cond(getParam<Real>("fluid_thermal_conductivity")),
    _solid_thermal_cond(getParam<Real>("solid_thermal_conductivity")),
    _fluid_heat_cap(getParam<Real>("fluid_heat_capacity")),
    _solid_heat_cap(getParam<Real>("solid_heat_capacity")),
    _fluid_thermal_exp(getParam<Real>("fluid_thermal_expansion")),
    _solid_thermal_exp(getParam<Real>("solid_thermal_expansion"))
{
}

void
toroThermalConstant::computeQpHeatCap()
{
  _c_f[_qp] = _fluid_heat_cap;
  _c_s[_qp] = _solid_heat_cap;
}

void
toroThermalConstant::computeQpThermalCond()
{
  _lambda_f[_qp] = _fluid_thermal_cond;
  _lambda_s[_qp] = _solid_thermal_cond;
}

void
toroThermalConstant::computeQpThermalExp()
{
  _beta_f[_qp] = _fluid_thermal_exp;
  _beta_s[_qp] = _solid_thermal_exp;
}
