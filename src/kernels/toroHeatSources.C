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

#include "toroHeatSources.h"

registerMooseObject("ToroApp", toroHeatSources);

template <>
InputParameters
validParams<toroHeatSources>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription(
      "Heat generation by radiogenic, shear heating and adiabatic processes.");
  params.addParam<Real>(
      "coeff_shear_heating", 0.0, "The coefficient in front of the shear heating generation.");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  return params;
}

toroHeatSources::toroHeatSources(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _coeff_Hs(getParam<Real>("coeff_shear_heating")),
    _rhoC_b(getDefaultMaterialProperty<Real>("bulk_specific_heat")),
    _radiogenic_heat(getDefaultMaterialProperty<Real>("radiogenic_heat_production")),
    _inelastic_heat(getDefaultMaterialProperty<Real>("inelastic_heat")),
    _adiabatic_heat(getDefaultMaterialProperty<Real>("adiabatic_heat")),
    _dinelastic_heat_dstrain(getDefaultMaterialProperty<RankTwoTensor>("dinelastic_heat_dstrain")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);

  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    mooseError("toroHeatSources: The number of displacement variables supplied must match the "
               "mesh dimension.");
}

/******************************************************************************/
/*                                  RESIDUALS                                 */
/******************************************************************************/

Real
toroHeatSources::computeQpResidual()
{
  Real Hr = _radiogenic_heat[_qp];
  Real Hs = _coeff_Hs * _inelastic_heat[_qp];
  Real Ha = _adiabatic_heat[_qp];
  Real heat_sources = Hr + Hs + Ha;
  if (_rhoC_b[_qp] != 0.0)
    heat_sources /= _rhoC_b[_qp];

  return -heat_sources * _test[_i][_qp];
}

/******************************************************************************/
/*                                  JACOBIAN                                  */
/******************************************************************************/

Real
toroHeatSources::computeQpJacobian()
{
  return 0.0;
}

/******************************************************************************/
/*                              OFF-DIAG JACOBIAN                             */
/******************************************************************************/

Real
toroHeatSources::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    if (jvar == _disp_var[i])
      return -_test[_i][_qp] * _coeff_Hs * (_dinelastic_heat_dstrain[_qp] * _grad_phi[_j][_qp])(i);

  return 0.0;
}
