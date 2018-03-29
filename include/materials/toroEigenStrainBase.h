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

#ifndef toroEIGENSTRAINBASE_H
#define toroEIGENSTRAINBASE_H

#include "Material.h"
#include "RankTwoTensor.h"

class toroEigenStrainBase;

template <>
InputParameters validParams<toroEigenStrainBase>();

class toroEigenStrainBase : public Material
{
public:
  toroEigenStrainBase(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();
  // Compute the eigenstrain and store in _eigenstrain
  virtual void computeQpEigenstrain() = 0;

  // Material property name for the eigenstrain tensor
  std::string _eigenstrain_name;

  // Stores the current total eigenstrain
  MaterialProperty<RankTwoTensor> & _eigenstrain;

  // Restartable data to check for the zeroth and first time steps for thermal calculations
  bool & _step_zero;
};

#endif // toroEIGENSTRAINBASE_H
