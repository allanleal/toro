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

#ifndef toroELASTICRHEOLOGY_H
#define toroELASTICRHEOLOGY_H

#include "Material.h"
#include "RankTwoTensor.h"
#include "DerivativeMaterialInterface.h"

class toroElasticRheology;

template <>
InputParameters validParams<toroElasticRheology>();

class toroElasticRheology : public DerivativeMaterialInterface<Material>
{
public:
  toroElasticRheology(const InputParameters & parameters);
  virtual ~toroElasticRheology() {}

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  virtual void computeQpRheologyProperties();
  virtual void computeQpStress();
  virtual void rotateIncrement();
  virtual void updateQpElasticState(RankTwoTensor & elastic_strain_increment);
  virtual void subtractEigenstrainIncrementFromStrain(RankTwoTensor & strain);

  bool _coupled_temp;
  const VariableValue & _temp;
  const Real _bulk_modulus;
  const Real _shear_modulus;

  const MaterialProperty<RankTwoTensor> & _strain_increment;
  const MaterialProperty<RankTwoTensor> & _rotation_increment;
  MaterialProperty<RankTwoTensor> & _elastic_strain;
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;
  MaterialProperty<RankFourTensor> & _elasticity_tensor;
  MaterialProperty<RankFourTensor> & _tangent_modulus;
  const MaterialProperty<RankTwoTensor> & _thermal_strain;
  const MaterialProperty<RankTwoTensor> & _thermal_strain_old;
  MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankTwoTensor> & _stress_old;
};

#endif // toroELASTICRHEOLOGY_H
