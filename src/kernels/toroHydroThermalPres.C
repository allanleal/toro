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

#include "toroHydroThermalPres.h"

registerMooseObject("ToroApp", toroHydroThermalPres);

template <>
InputParameters
validParams<toroHydroThermalPres>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Thermal pressurisation coupling kernel.");
  params.addRequiredCoupledVar("temperature", "The temperature variable.");
  return params;
}

toroHydroThermalPres::toroHydroThermalPres(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _temp_dot(coupledDot("temperature")),
    _dtemp_dot_dtemp(coupledDotDu("temperature")),
    _biot(getDefaultMaterialProperty<Real>("biot_coefficient")),
    _thermal_exp(getDefaultMaterialProperty<Real>("thermal_expansion_coefficient"))
{
}

/******************************************************************************/
/*                                  RESIDUALS                                 */
/******************************************************************************/

Real
toroHydroThermalPres::computeQpResidual()
{
  return -_biot[_qp] * _thermal_exp[_qp] * _temp_dot[_qp] * _test[_i][_qp];
}

/******************************************************************************/
/*                                  JACOBIAN                                  */
/******************************************************************************/

Real
toroHydroThermalPres::computeQpJacobian()
{
  return 0.0;
}

/******************************************************************************/
/*                              OFF-DIAG JACOBIAN                             */
/******************************************************************************/

Real
toroHydroThermalPres::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _temp_var)
    return -_biot[_qp] * _thermal_exp[_qp] * _dtemp_dot_dtemp[_qp] * _phi[_j][_qp] * _test[_i][_qp];

  return 0.0;
}
