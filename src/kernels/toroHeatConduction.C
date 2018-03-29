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

#include "toroHeatConduction.h"

registerMooseObject("ToroApp", toroHeatConduction);

template <>
InputParameters
validParams<toroHeatConduction>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Heat conduction kernel.");
  return params;
}

toroHeatConduction::toroHeatConduction(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _thermal_diff(getMaterialProperty<Real>("thermal_diffusivity"))
{
}

/******************************************************************************/
/*                                  RESIDUALS                                 */
/******************************************************************************/

Real
toroHeatConduction::computeQpResidual()
{
  return _thermal_diff[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];
}

/******************************************************************************/
/*                                  JACOBIAN                                  */
/******************************************************************************/

Real
toroHeatConduction::computeQpJacobian()
{
  return _thermal_diff[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}

/******************************************************************************/
/*                              OFF-DIAG JACOBIAN                             */
/******************************************************************************/

// Real
// toroHeatConduction::computeQpOffDiagJacobian(unsigned int jvar)
// {
//   return 0.0;
// }
