//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ReaktoroMULTIAPP_H
#define ReaktoroMULTIAPP_H

#include "TransientMultiApp.h"
#include "BoundaryRestrictable.h"

#include "Reaktoro/Reaktoro.hpp"

class ReaktoroMultiApp;

template <>
InputParameters validParams<ReaktoroMultiApp>();

/**
 * Automatically generates Sub-App positions from positions in the master app's mesh.
 */
class ReaktoroMultiApp : public TransientMultiApp, public BoundaryRestrictable
{
public:
  ReaktoroMultiApp(const InputParameters & parameters);

  std::vector<std::string> getSpecies();

protected:
  virtual void fillPositions() override;

  virtual void initialSetup() override;

  virtual bool solveStep(Real dt, Real target_time, bool auto_advance = true) override;

private:
  Reaktoro::ChemicalSystem _system;
  Reaktoro::ChemicalState _state_bc;
  Reaktoro::ChemicalState _state_ic;
};

#endif // ReaktoroMULTIAPP_H
