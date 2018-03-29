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

#include "toroSmallStrain.h"
#include "Assembly.h"
#include "libmesh/quadrature.h"

registerMooseObject("ToroApp", toroSmallStrain);

template <>
InputParameters
validParams<toroSmallStrain>()
{
  InputParameters params = validParams<toroStrainBase>();
  params.addClassDescription(
      "Compute a strain increment and rotation increment for small strains.");
  return params;
}

toroSmallStrain::toroSmallStrain(const InputParameters & parameters) : toroStrainBase(parameters) {}

void
toroSmallStrain::computeProperties()
{
  Real volumetric_strain = 0.0;
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    RankTwoTensor total_strain_increment;
    computeTotalStrainIncrement(total_strain_increment);

    _strain_increment[_qp] = total_strain_increment;

    if (_volumetric_locking_correction)
      volumetric_strain += total_strain_increment.trace() * _JxW[_qp] * _coord[_qp];
  }
  if (_volumetric_locking_correction)
    volumetric_strain /= _current_elem_volume;

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    Real trace = _strain_increment[_qp].trace();
    if (_volumetric_locking_correction)
    {
      _strain_increment[_qp](0, 0) += (volumetric_strain - trace) / 3.0;
      _strain_increment[_qp](1, 1) += (volumetric_strain - trace) / 3.0;
      _strain_increment[_qp](2, 2) += (volumetric_strain - trace) / 3.0;
    }

    // Substract EigenStrain increment
    subtractEigenstrainIncrement(_strain_increment[_qp]);

    // Update strain in intermediate configuration: rotations are not needed
    _total_strain[_qp] = _total_strain_old[_qp] + _strain_increment[_qp];
  }
}

void
toroSmallStrain::computeTotalStrainIncrement(RankTwoTensor & total_strain_increment)
{
  // Deformation gradient
  RankTwoTensor A(
      (*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]); // Deformation gradient
  RankTwoTensor Fbar((*_grad_disp_old[0])[_qp],
                     (*_grad_disp_old[1])[_qp],
                     (*_grad_disp_old[2])[_qp]); // Old Deformation gradient

  _deformation_gradient[_qp] = A;
  _deformation_gradient[_qp].addIa(1.0);

  A -= Fbar; // A = grad_disp - grad_disp_old

  total_strain_increment = 0.5 * (A + A.transpose());
}
