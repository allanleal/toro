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

  params.addParam<std::string>("pressure_units", "Pa", "Units pressure is in (Pa)");
  params.addParam<std::string>("temperature_units", "K", "Units temperature is in (K)");

  params.addCoupledVar("temperature", 300, "The temperature.");
  params.addCoupledVar("pressure", 1e5, "The pressure.");

  params.addRequiredParam<std::vector<std::string>>("substance_names", "Names of the substances, these go with substance_amounts and substance_units");
  params.addRequiredParam<std::vector<Real>>("substance_amounts", "The amount of the substance_names");
  params.addRequiredParam<std::vector<std::string>>("substance_units", "The units to use for each amount (kg, mol)");

  return params;
}

SpeciesAux::SpeciesAux(const InputParameters & parameters)
    : AuxKernel(parameters), _n_vars(coupledComponents("species")),
      _reaktoro_system(*getCheckedPointerParam<Reaktoro::ChemicalSystem *>("reaktoro_system")),
      _substance_names(getParam<std::vector<std::string>>("substance_names")),
      _substance_amounts(getParam<std::vector<Real>>("substance_amounts")),
      _substance_units(getParam<std::vector<std::string>>("substance_units")),
      _temp(coupledValue("temperature")),
      _pressure(coupledValue("pressure"))
{
  /*
  EquilibriumProblem problem_bc(_system);
  problem_bc.setTemperature(25, "celsius");
  problem_bc.setPressure(1, "bar");
  problem_bc.add("H2O", 1, "kg");
  problem_bc.add("NaCl", 0.1, "mol");
  */

  std::cout << "n_vars: "<<_n_vars<<std::endl;


  Reaktoro::EquilibriumProblem reaktoro_problem(_reaktoro_system);

  reaktoro_problem.setTemperature(300, getParam<std::string>("temperature_units"));
  reaktoro_problem.setPressure(1e5, getParam<std::string>("pressure_units"));

  for (unsigned int i = 0; i < _substance_names.size(); i++)
    reaktoro_problem.add(_substance_names[i], _substance_amounts[i], _substance_units[i]);

  _state = equilibrate(reaktoro_problem);

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

  if (_temp[_qp] > 399)
    _state.output("results.txt");
}
