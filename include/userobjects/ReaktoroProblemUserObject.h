//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef REAKTOROPROBLEMUSEROBJECT_H
#define REAKTOROPROBLEMUSEROBJECT_H

#include "GeneralUserObject.h"

#include "Reaktoro/Reaktoro.hpp"

class ReaktoroProblemUserObject;

template <>
InputParameters validParams<ReaktoroProblemUserObject>();

/// Records all post processor data in a CSV file.
class ReaktoroProblemUserObject : public GeneralUserObject
{
public:
  ReaktoroProblemUserObject(const InputParameters & parameters);

  std::vector<std::string> getElementNames() const;
  std::vector<std::string> getSpeciesNames() const;
  Reaktoro::ChemicalState getInitialState() const;

  virtual void initialize() override{};
  virtual void execute() override;
  virtual void finalize() override{};

  Reaktoro::ChemicalSystem _system;
  std::unique_ptr<Reaktoro::EquilibriumProblem> _reaktoro_problem;
  Reaktoro::ChemicalState _state;

  const std::vector<std::string> & _substance_names;
  const std::vector<Real> & _substance_amounts;
  const std::vector<std::string> & _substance_units;

  const Real & _initial_temperature;
  const Real & _initial_pressure;
};

#endif
