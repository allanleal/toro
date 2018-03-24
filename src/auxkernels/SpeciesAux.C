//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SpeciesAux.h"
#include "ReaktoroProblemUserObject.h"

registerMooseObject("ToroApp", SpeciesAux);

template <>
InputParameters
validParams<SpeciesAux>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addRequiredCoupledVar("species", "unknown (nl-variable)");
  params.addPrivateParam<Reaktoro::ChemicalSystem *>("reaktoro_system");

  params.addCoupledVar("temperature", 300, "The temperature.");
  params.addCoupledVar("pressure", 1e5, "The pressure.");

  return params;
}

SpeciesAux::SpeciesAux(const InputParameters & parameters)
    : AuxKernel(parameters), _n_vars(coupledComponents("species")),
      _reaktoro_problem(getUserObject<ReaktoroProblemUserObject>("reaktoro_problem")),
      _temp(coupledValue("temperature")),
      _pressure(coupledValue("pressure")),
      _state(_reaktoro_problem.getInitialState())
{
  /*
  EquilibriumProblem problem_bc(_system);
  problem_bc.setTemperature(25, "celsius");
  problem_bc.setPressure(1, "bar");
  problem_bc.add("H2O", 1, "kg");
  problem_bc.add("NaCl", 0.1, "mol");
  */

  for (unsigned int i = 0; i < _n_vars; i++)
    _vars.push_back(dynamic_cast<MooseVariable *>(getVar("species", i)));
}

Real
SpeciesAux::computeValue()
{
  std::vector<Real> values(_n_vars);

  computeVarValues(values);

  for (unsigned int i = 0; i < _n_vars; i++)
    _vars[i]->setNodalValue(values[i]);

  return _vars[0]->nodalValue()[0];
}

void
SpeciesAux::computeVarValues(std::vector<Real> & values)
{
  // _state_bc.setSpeciesAmounts(i, val);
  /*
  for (auto i = beginIndex(species_amounts); i < species_amounts.size(); i++)
    _state.setSpeciesAmount(i, species_amounts[i]);
  */

  _state.setTemperature(_temp[_qp]);
  _state.setPressure(_pressure[_qp]);

  equilibrate(_state);

  const auto & species_amounts = _state.speciesAmounts();

  for (auto i = beginIndex(species_amounts); i < species_amounts.size(); i++)
    values[i] = species_amounts[i];

}
