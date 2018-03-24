//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ADDSPECIESACTION_H
#define ADDSPECIESACTION_H

#include "InputParameters.h"
#include "Action.h"

#include "Reaktoro/Reaktoro.hpp"

// Forward declaration
class AddSpeciesAction;

/**
 * Automatically generates all variables to model a polycrystal with op_num orderparameters
 */
class AddSpeciesAction : public Action
{
public:
  AddSpeciesAction(const InputParameters & params);

  virtual void act();

private:
  void createChemicalSystem();

  std::vector<VariableName> _nonlinear_element_names;
  std::vector<VariableName> _element_names;
  std::vector<VariableName> _species_names;
};

template <>
InputParameters validParams<AddSpeciesAction>();

#endif // ADDSPECIESACTION_H
