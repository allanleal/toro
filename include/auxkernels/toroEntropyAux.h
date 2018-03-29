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
#ifndef toroENTROPYAUX_H
#define toroENTROPYAUX_H

#include "AuxKernel.h"
#include "DerivativeMaterialInterface.h"

class toroEntropyAux;
template <>
InputParameters validParams<toroEntropyAux>();

class toroEntropyAux : public DerivativeMaterialInterface<AuxKernel>
{
public:
  toroEntropyAux(const InputParameters & parameters);
  virtual ~toroEntropyAux() {}

protected:
  virtual Real computeValue();
  const VariableValue & _var_old;
  const VariableValue & _var_older;

  const PostprocessorValue & _pp_max_var;
  const PostprocessorValue & _pp_min_var;
};

#endif // toroENTROPY_H
