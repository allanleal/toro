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

#ifndef toroHYDROPOROMECH_H
#define toroHYDROPOROMECH_H

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"

class toroHydroPoroMech;

template <>
InputParameters validParams<toroHydroPoroMech>();

class toroHydroPoroMech : public DerivativeMaterialInterface<Kernel>
{
public:
  toroHydroPoroMech(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  unsigned int _ndisp;
  const MaterialProperty<Real> & _poro_mech;
  const MaterialProperty<Real> & _poro_mech_jac;
  std::vector<unsigned int> _disp_var;
};

#endif // toroHYDROPOROMECH_H
