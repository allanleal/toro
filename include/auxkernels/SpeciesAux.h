//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECIESAUX_H
#define SPECIESAUX_H

#include "ReaktoroProblemUserObject.h"

#include "AuxKernel.h"
#include "Reaktoro/Reaktoro.hpp"

class SpeciesAux;

template <>
InputParameters validParams<SpeciesAux>();

/**
 * Aux kernel that updated values of coupled variables
 */
class SpeciesAux : public AuxKernel
{
public:
  SpeciesAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  unsigned int _n_elements;
  unsigned int _n_species;
  std::vector<MooseVariable *> _nonlinear_element_vars;
  std::vector<MooseVariable *> _element_vars;
  std::vector<MooseVariable *> _species_vars;

  const ReaktoroProblemUserObject & _reaktoro_problem;

  const VariableValue & _temp;
  const VariableValue & _pressure;

  Reaktoro::ChemicalState _state;
};

#endif /* SPECIESAUX_H */
