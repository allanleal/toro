//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ReaktoroProblemUserObject.h"
#include "libmesh/perf_log.h"

#include "Reaktoro/Reaktoro.hpp"

#include <ostream>

registerMooseObject("ToroApp", ReaktoroProblemUserObject);

template <>
InputParameters
validParams<ReaktoroProblemUserObject>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.set<ExecFlagEnum>("execute_on") = EXEC_FINAL;
  params.addParam<std::string>("outfile", "perflog.csv", "name of perf log output file");
  params.addClassDescription("Dumps perlog information to a csv file for further analysis.");

  params.addParam<Real>("initial_temperature",
                        300,
                        "The temperature to use in the calculation of the initial equilibrium");
  params.addParam<Real>(
      "initial_pressure", 1e5, "The pressure to use in the calculation of the initial equilibrium");

  params.addParam<std::string>("pressure_units", "Pa", "Units pressure is in (Pa)");
  params.addParam<std::string>("temperature_units", "K", "Units temperature is in (K)");

  params.addRequiredParam<std::vector<std::string>>(
      "substance_names",
      "Names of the substances, these go with substance_amounts and substance_units");
  params.addRequiredParam<std::vector<Real>>("substance_amounts",
                                             "The amount of the substance_names");
  params.addRequiredParam<std::vector<std::string>>("substance_units",
                                                    "The units to use for each amount (kg, mol)");

  return params;
}

ReaktoroProblemUserObject::ReaktoroProblemUserObject(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _substance_names(getParam<std::vector<std::string>>("substance_names")),
    _substance_amounts(getParam<std::vector<Real>>("substance_amounts")),
    _substance_units(getParam<std::vector<std::string>>("substance_units")),
    _initial_temperature(getParam<Real>("initial_temperature")),
    _initial_pressure(getParam<Real>("initial_pressure"))
{
  Reaktoro::Database database("supcrt98.xml");

  Reaktoro::ChemicalEditor editor(database);
  editor.addAqueousPhase({"H2O(l)", "H+", "OH-", "Na+", "Cl-"});

  _system = Reaktoro::ChemicalSystem(editor);

  _reaktoro_problem = libmesh_make_unique<Reaktoro::EquilibriumProblem>(_system);

  _reaktoro_problem->setTemperature(_initial_temperature,
                                    getParam<std::string>("temperature_units"));
  _reaktoro_problem->setPressure(_initial_pressure, getParam<std::string>("pressure_units"));

  for (unsigned int i = 0; i < _substance_names.size(); i++)
    _reaktoro_problem->add(_substance_names[i], _substance_amounts[i], _substance_units[i]);
}

std::vector<std::string>
ReaktoroProblemUserObject::getElementNames() const
{
  return names(_system.elements());
}

std::vector<std::string>
ReaktoroProblemUserObject::getSpeciesNames() const
{
  return names(_system.species());
}

Reaktoro::ChemicalState
ReaktoroProblemUserObject::getInitialState() const
{
  return equilibrate(*_reaktoro_problem);
}

void
ReaktoroProblemUserObject::execute()
{
}
