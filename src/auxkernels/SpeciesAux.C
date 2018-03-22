//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SpeciesAux.h"
#include "Reaktoro/Reaktoro.hpp"

registerMooseObject("ToroApp", SpeciesAux);

template <>
InputParameters
validParams<SpeciesAux>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addRequiredCoupledVar("species", "unknown (nl-variable)");
  params.addPrivateParam<Reaktoro::ChemicalSystem *>("reaktoro_system");

  return params;
}

SpeciesAux::SpeciesAux(const InputParameters & parameters)
    : AuxKernel(parameters), _n_vars(coupledComponents("species")),
      _reaktoro_system(*getCheckedPointerParam<Reaktoro::ChemicalSystem *>("reaktoro_system"))
{
  /*
  EquilibriumProblem problem_bc(_system);
  problem_bc.setTemperature(25, "celsius");
  problem_bc.setPressure(1, "bar");
  problem_bc.add("H2O", 1, "kg");
  problem_bc.add("NaCl", 0.1, "mol");
  */

  Reaktoro::EquilibriumProblem reaktoro_problem(_reaktoro_system);
  reaktoro_problem.setTemperature(25, "celsius");
  reaktoro_problem.setPressure(1, "bar");
  reaktoro_problem.add("H2O", 1, "kg");

  _state = equilibrate(reaktoro_problem);

  for (unsigned int i = 0; i < _n_vars; i++)
    _vars.push_back(dynamic_cast<MooseVariable *>(getVar("species", i)));
}

void
SpeciesAux::compute()
{
  std::vector<Real> values(_n_vars);

  computeVarValues(values);

  for (unsigned int i = 0; i < _n_vars; i++)
    _vars[i]->setNodalValue(values[i]);

  _var.setNodalValue(0.0);
}

Real
SpeciesAux::computeValue()
{
  return 0.0;
}

void
SpeciesAux::computeVarValues(std::vector<Real> & values)
{
  // _state_bc.setSpeciesAmounts(i, val);
  /*
  for (auto i = beginIndex(species_amounts); i < species_amounts.size(); i++)
    _state.setSpeciesAmount(i, species_amounts[i]);
  */

  equilibrate(_state);

  const auto & species_amounts = _state.speciesAmounts();

  for (auto i = beginIndex(species_amounts); i < species_amounts.size(); i++)
    values[i] = species_amounts[i];
}
