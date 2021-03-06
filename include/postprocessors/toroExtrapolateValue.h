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

#ifndef toroEXTRAPOLATEVALUE_H
#define toroEXTRAPOLATEVALUE_H

#include "ElementExtremeValue.h"
#include "DerivativeMaterialInterface.h"

class toroExtrapolateValue;

template <>
InputParameters validParams<toroExtrapolateValue>();

class toroExtrapolateValue : public DerivativeMaterialInterface<ElementExtremeValue>
{
public:
  toroExtrapolateValue(const InputParameters & parameters);

protected:
  virtual void computeQpValue() override;

  const VariableValue & _value_old;
  const VariableValue & _value_older;
};

#endif // toroEXTRAPOLATEVALUE_H
