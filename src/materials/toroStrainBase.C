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

#include "toroStrainBase.h"
#include "MooseMesh.h"
#include "Assembly.h"

template <>
InputParameters
validParams<toroStrainBase>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("Base class for calculating the strain tensor.");
  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system.");
  params.addParam<bool>(
      "volumetric_locking_correction", false, "Flag to correct volumetric locking.");
  params.addParam<std::vector<MaterialPropertyName>>(
      "eigenstrain_names", "List of eigenstrains to be applied in this strain calculation.");
  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

toroStrainBase::toroStrainBase(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _ndisp(coupledComponents("displacements")),
    _grad_disp(3),
    _grad_disp_old(3),
    _total_strain(declareProperty<RankTwoTensor>("total_strain")),
    _total_strain_old(getMaterialPropertyOld<RankTwoTensor>("total_strain")),
    _eigenstrain_names(getParam<std::vector<MaterialPropertyName>>("eigenstrain_names")),
    _eigenstrains(_eigenstrain_names.size()),
    _eigenstrains_old(_eigenstrain_names.size()),
    _strain_increment(declareProperty<RankTwoTensor>("strain_increment")),
    _rotation_increment(declareProperty<RankTwoTensor>("rotation_increment")),
    _deformation_gradient(declareProperty<RankTwoTensor>("deformation_gradient")),
    _volumetric_locking_correction(getParam<bool>("volumetric_locking_correction")),
    _current_elem_volume(_assembly.elemVolume())
{
  for (unsigned int i = 0; i < _eigenstrains.size(); ++i)
  {
    _eigenstrains[i] = &getMaterialProperty<RankTwoTensor>(_eigenstrain_names[i]);
    _eigenstrains_old[i] = &getMaterialPropertyOld<RankTwoTensor>(_eigenstrain_names[i]);
  }

  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    mooseError(
        "The number of variables supplied in 'displacements' must match the mesh dimension.");

  // fetch coupled variables and gradients (as stateful properties if necessary)
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _grad_disp[i] = &coupledGradient("displacements", i);
    if (_fe_problem.isTransient())
      _grad_disp_old[i] = &coupledGradientOld("displacements", i);
    else
      _grad_disp_old[i] = &_grad_zero;
  }

  // set unused dimensions to zero
  for (unsigned i = _ndisp; i < 3; ++i)
  {
    _grad_disp[i] = &_grad_zero;
    _grad_disp_old[i] = &_grad_zero;
  }

  if (getParam<bool>("use_displaced_mesh"))
    mooseError("The strain calculator needs to run on the undisplaced mesh.");
}

void
toroStrainBase::initQpStatefulProperties()
{
  _total_strain[_qp].zero();
  _deformation_gradient[_qp].zero();
  _deformation_gradient[_qp].addIa(1.0);

  _rotation_increment[_qp].zero();
  _rotation_increment[_qp].addIa(1.0);
}

void
toroStrainBase::subtractEigenstrainIncrement(RankTwoTensor & strain)
{
  for (unsigned int i = 0; i < _eigenstrains.size(); ++i)
  {
    strain -= (*_eigenstrains[i])[_qp];
    strain += (*_eigenstrains_old[i])[_qp];
  }
}
