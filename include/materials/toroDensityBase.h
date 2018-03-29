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

#ifndef toroDENSITYBASE_H
#define toroDENSITYBASE_H

#include "Material.h"
#include "DerivativeMaterialInterface.h"

class toroDensityBase;

template <>
InputParameters validParams<toroDensityBase>();

class toroDensityBase : public DerivativeMaterialInterface<Material>
{
public:
  toroDensityBase(const InputParameters & parameters);
  virtual ~toroDensityBase() {}

protected:
  virtual void computeQpGravity();

  const VariableValue & _porosity;
  bool _has_gravity;
  Real _g;
  Real _fluid_density;
  Real _solid_density;

  MaterialProperty<RealVectorValue> & _gravity;
  MaterialProperty<Real> & _rho_f;
  MaterialProperty<Real> & _rho_s;
  MaterialProperty<Real> & _rho_b;
};

#endif // toroDENSITYBASE_H
