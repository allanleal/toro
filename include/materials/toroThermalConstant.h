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

#ifndef toroTHERMALCONSTANT_H
#define toroTHERMALCONSTANT_H

#include "toroThermalBase.h"

class toroThermalConstant;

template <>
InputParameters validParams<toroThermalConstant>();

class toroThermalConstant : public toroThermalBase
{
public:
  toroThermalConstant(const InputParameters & parameters);

protected:
  virtual void computeQpHeatCap() override;
  virtual void computeQpThermalCond() override;
  virtual void computeQpThermalExp() override;

  Real _fluid_thermal_cond;
  Real _solid_thermal_cond;
  Real _fluid_heat_cap;
  Real _solid_heat_cap;
  Real _fluid_thermal_exp;
  Real _solid_thermal_exp;
};

#endif // toroTHERMALCONSTANT_H
