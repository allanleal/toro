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

#include "toroElasticRheology.h"

registerMooseObject("ToroApp", toroElasticRheology);

template <>
InputParameters
validParams<toroElasticRheology>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("Class calculating stress based on an elastic rheology.");
  params.addCoupledVar("temperature", "The temperature variable.");
  params.addRequiredRangeCheckedParam<Real>(
      "bulk_modulus", "bulk_modulus > 0.0", "The drained bulk modulus of the material.");
  params.addRequiredRangeCheckedParam<Real>(
      "shear_modulus", "shear_modulus > 0.0", "The shear modulus of the material.");
  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

toroElasticRheology::toroElasticRheology(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _coupled_temp(isCoupled("temperature")),
    _temp(_coupled_temp ? coupledValue("temperature") : _zero),
    _bulk_modulus(getParam<Real>("bulk_modulus")),
    _shear_modulus(getParam<Real>("shear_modulus")),
    _strain_increment(getMaterialProperty<RankTwoTensor>("strain_increment")),
    _rotation_increment(getMaterialProperty<RankTwoTensor>("rotation_increment")),
    _elastic_strain(declareProperty<RankTwoTensor>("elastic_strain")),
    _elastic_strain_old(getMaterialPropertyOld<RankTwoTensor>("elastic_strain")),
    _elasticity_tensor(declareProperty<RankFourTensor>("elasticity_tensor")),
    _tangent_modulus(declareProperty<RankFourTensor>("tangent_modulus")),
    _thermal_strain(getDefaultMaterialProperty<RankTwoTensor>("thermal_strain")),
    _thermal_strain_old(getMaterialPropertyOld<RankTwoTensor>("thermal_strain")),
    _stress(declareProperty<RankTwoTensor>("stress")),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>("stress"))
{
  if (getParam<bool>("use_displaced_mesh"))
    mooseError("The stress calculator needs to run on the undisplaced mesh.");
}

void
toroElasticRheology::initQpStatefulProperties()
{
  _stress[_qp].zero();
  _elastic_strain[_qp].zero();
}

void
toroElasticRheology::computeQpProperties()
{
  computeQpRheologyProperties();
  computeQpStress();
  rotateIncrement();
}

void
toroElasticRheology::computeQpRheologyProperties()
{
  std::vector<Real> iso_const(2);
  iso_const[0] = _bulk_modulus - 2.0 / 3.0 * _shear_modulus;
  iso_const[1] = _shear_modulus;
  _elasticity_tensor[_qp].fillFromInputVector(iso_const, RankFourTensor::symmetric_isotropic);
}

void
toroElasticRheology::computeQpStress()
{
  RankTwoTensor elastic_strain_increment = RankTwoTensor();

  updateQpElasticState(elastic_strain_increment);
  _elastic_strain[_qp] = _elastic_strain_old[_qp] + elastic_strain_increment;
}

void
toroElasticRheology::rotateIncrement()
{
  // Perform rotations
  _elastic_strain[_qp] =
      _rotation_increment[_qp] * _elastic_strain[_qp] * _rotation_increment[_qp].transpose();
  _stress[_qp] = _rotation_increment[_qp] * _stress[_qp] * _rotation_increment[_qp].transpose();
}

void
toroElasticRheology::updateQpElasticState(RankTwoTensor & elastic_strain_increment)
{
  elastic_strain_increment = _strain_increment[_qp];
  // Remove the thermal strain increment
  subtractEigenstrainIncrementFromStrain(elastic_strain_increment);
    _stress[_qp] = _stress_old[_qp] + _elasticity_tensor[_qp] * elastic_strain_increment;

  _tangent_modulus[_qp] = _elasticity_tensor[_qp];
}

void
toroElasticRheology::subtractEigenstrainIncrementFromStrain(RankTwoTensor & strain)
{
  strain -= _thermal_strain[_qp];
  strain += _thermal_strain_old[_qp];
}
