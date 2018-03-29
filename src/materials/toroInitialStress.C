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

#include "toroInitialStress.h"
#include "Function.h"

registerMooseObject("ToroApp", toroInitialStress);

template <>
InputParameters
validParams<toroInitialStress>()
{
  InputParameters params = validParams<toroEigenStrainBase>();
  params.addClassDescription("Computes an eigenstrain from an initial stress");
  params.addRequiredRangeCheckedParam<Real>(
      "bulk_modulus", "bulk_modulus > 0.0", "The drained bulk modulus of the material.");
  params.addRequiredRangeCheckedParam<Real>(
      "shear_modulus", "shear_modulus > 0.0", "The shear modulus of the material.");
  params.addRequiredParam<std::vector<FunctionName>>(
      "initial_stress",
      "A list of functions describing the initial stress.  If provided, there must be 3 of these, "
      "corresponding to the xx, yy, zz components respectively (principle stresses). To compute "
      "the eigenstrain correctly, your elasticity tensor should not be time-varying in the first "
      "timestep");
  return params;
}

toroInitialStress::toroInitialStress(const InputParameters & parameters)
  : toroEigenStrainBase(parameters),
    _bulk_modulus(getParam<Real>("bulk_modulus")),
    _shear_modulus(getParam<Real>("shear_modulus")),
    _eigenstrain_old(getMaterialPropertyOld<RankTwoTensor>(_eigenstrain_name))
{
  const std::vector<FunctionName> & fcn_names(
      getParam<std::vector<FunctionName>>("initial_stress"));
  const unsigned num = fcn_names.size();

  if (num != 2.0 * LIBMESH_DIM)
    mooseError("toroInitialStress: ",
               2.0 * LIBMESH_DIM,
               " initial stress functions must be provided. You supplied ",
               num,
               "\n");

  _initial_stress_fcn.resize(num);
  for (unsigned i = 0; i < num; ++i)
    _initial_stress_fcn[i] = &getFunctionByName(fcn_names[i]);
}

void
toroInitialStress::computeQpEigenstrain()
{
  if (_t_step == 1)
  {
    RankFourTensor Eijkl0 = RankFourTensor();
    std::vector<Real> iso_const(2);
    iso_const[0] = _bulk_modulus - 2.0 / 3.0 * _shear_modulus;
    iso_const[1] = _shear_modulus;
    Eijkl0.fillFromInputVector(iso_const, RankFourTensor::symmetric_isotropic);

    RankTwoTensor initial_stress = RankTwoTensor();
    for (unsigned i = 0; i < LIBMESH_DIM; ++i)
      for (unsigned j = i; j < LIBMESH_DIM; ++j)
        initial_stress(i, j) = initial_stress(j, i) =
            _initial_stress_fcn[i * LIBMESH_DIM + j - 1 * (i > 0) - 2 * (i > 1)]->value(
                _t, _q_point[_qp]);

    _eigenstrain[_qp] = -Eijkl0.invSymm() * initial_stress;
  }
  else
    _eigenstrain[_qp] = _eigenstrain_old[_qp];
}
