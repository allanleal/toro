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

#ifndef toroHEATSOURCES_H
#define toroHEATSOURCES_H

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"
#include "RankTwoTensor.h"

class toroHeatSources;

template <>
InputParameters validParams<toroHeatSources>();

class toroHeatSources : public DerivativeMaterialInterface<Kernel>
{
public:
  toroHeatSources(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  Real _coeff_Hs;

  const MaterialProperty<Real> & _rhoC_b;
  const MaterialProperty<Real> & _radiogenic_heat;
  const MaterialProperty<Real> & _inelastic_heat;
  const MaterialProperty<Real> & _adiabatic_heat;
  // const MaterialProperty<Real> & _dinelastic_heat_dtemp;
  const MaterialProperty<RankTwoTensor> & _dinelastic_heat_dstrain;
  unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;
};

#endif // toroHEATSOURCES_H
