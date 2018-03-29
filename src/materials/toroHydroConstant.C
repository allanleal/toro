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

#include "toroHydroConstant.h"

registerMooseObject("ToroApp", toroHydroConstant);

template <>
InputParameters
validParams<toroHydroConstant>()
{
  InputParameters params = validParams<toroHydroBase>();
  params.addClassDescription("Constant thermal properties.");
  params.addRequiredParam<Real>("permeability", "The permeability of the matrix.");
  params.addRequiredParam<Real>("fluid_viscosity", "The viscosity of the fluid.");
  params.addParam<Real>("fluid_modulus", "The bulk modulus of the fluid phase.");
  params.addParam<Real>("solid_modulus", "The bulk modulus of the solid phase.");
  return params;
}

toroHydroConstant::toroHydroConstant(const InputParameters & parameters)
  : toroHydroBase(parameters),
    _perm(getParam<Real>("permeability")),
    _fluid_viscosity(getParam<Real>("fluid_viscosity"))
{
  if (isParamValid("fluid_modulus"))
    _fluid_compr = 1.0 / getParam<Real>("fluid_modulus");
  else
    _fluid_compr = 0.0;

  if (isParamValid("solid_modulus"))
    _solid_compr = 1.0 / getParam<Real>("solid_modulus");
  else
    _solid_compr = 0.0;
}

void
toroHydroConstant::computeQpFluidCompressibility()
{
  _C_f[_qp] = _fluid_compr;
}

void
toroHydroConstant::computeQpSolidCompressibility()
{
  _C_s[_qp] = _solid_compr;
}

void
toroHydroConstant::computeQpPermeability()
{
  _k[_qp] = _perm;
}

void
toroHydroConstant::computeQpFluidViscosity()
{
  _eta_f[_qp] = _fluid_viscosity;
}
