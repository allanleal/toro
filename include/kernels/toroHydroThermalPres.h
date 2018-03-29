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

#ifndef toroHYDROTHERMALPRES_H
#define toroHYDROTHERMALPRES_H

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"

class toroHydroThermalPres;

template <>
InputParameters validParams<toroHydroThermalPres>();

class toroHydroThermalPres : public DerivativeMaterialInterface<Kernel>
{
public:
  toroHydroThermalPres(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const VariableValue & _temp_dot;
  const VariableValue & _dtemp_dot_dtemp;
  const MaterialProperty<Real> & _biot;
  const MaterialProperty<Real> & _thermal_exp;
  unsigned int _temp_var;
};

#endif // toroHYDROTHERMALPRES_H
