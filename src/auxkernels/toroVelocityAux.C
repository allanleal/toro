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

#include "toroVelocityAux.h"

registerMooseObject("ToroApp", toroVelocityAux);

template <>
InputParameters
validParams<toroVelocityAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("displacement",
                               "The displacement variable to calculate the velocity.");
  params.addClassDescription("Calculates a component of the solid velocity based on displacement.");
  return params;
}

toroVelocityAux::toroVelocityAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _disp(coupledValue("displacement")),
    _disp_old(_is_transient ? coupledValueOld("displacement") : _zero)
{
}

Real
toroVelocityAux::computeValue()
{
  Real inv_dt = 1.0;
  if (_is_transient)
    inv_dt /= _dt;

  return (_disp[_qp] - _disp_old[_qp]) * inv_dt;
}
