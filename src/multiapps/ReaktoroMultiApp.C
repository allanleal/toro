//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "ReaktoroMultiApp.h"
#include "MooseMesh.h"
#include "FEProblem.h"
#include "Reaktoro.hpp"

registerMooseObject("MooseApp", ReaktoroMultiApp);

template <>
InputParameters
validParams<ReaktoroMultiApp>()
{
  InputParameters params = validParams<TransientMultiApp>();

  params += validParams<BoundaryRestrictable>();

  params.suppressParameter<std::vector<Point>>("positions");
  params.suppressParameter<std::vector<FileName>>("positions_file");

  // Turn off the base class position parameter
  params.set<bool>("use_positions") = false;

  return params;
}

ReaktoroMultiApp::ReaktoroMultiApp(const InputParameters & parameters)
  : TransientMultiApp(parameters), BoundaryRestrictable(this, true) // true for applying to nodesets
{
  //fillPositions();
  using namespace Reaktoro;

  Database database("supcrt98.xml");

  ChemicalEditor editor(database);
  editor.addAqueousPhase({"H2O(l)", "H+", "OH-", "Na+", "Cl-"});

  ChemicalSystem system(editor);

  EquilibriumProblem problem_bc(system);
  problem_bc.setTemperature(25, "celsius");
  problem_bc.setPressure(1, "bar");
  problem_bc.add("H2O", 1, "kg");
  problem_bc.add("NaCl", 0.1, "mol");

  EquilibriumProblem problem_ic(system);
  problem_ic.setTemperature(25, "celsius");
  problem_ic.setPressure(1, "bar");
  problem_ic.add("H2O", 1, "kg");

  ChemicalState state_ic = equilibrate(problem_ic);
  ChemicalState state_bc = equilibrate(problem_bc);

  state_bc.output("state_bc.txt");
  state_ic.output("state_ic.txt");

  fillPositions();
}

void
ReaktoroMultiApp::initialSetup()
{
  init(_positions.size());

  MultiApp::initialSetup();
}

void
ReaktoroMultiApp::fillPositions()
{
  MooseMesh & master_mesh = _fe_problem.mesh();

  const std::set<BoundaryID> & bids = boundaryIDs();

  for (const auto & boundary_id : bids)
  {
    // Grab the nodes on the boundary ID and add a Sub-App at each one.
    const std::vector<dof_id_type> & boundary_node_ids = master_mesh.getNodeList(boundary_id);

    for (const auto & boundary_node_id : boundary_node_ids)
    {
      Node & node = master_mesh.nodeRef(boundary_node_id);

      _positions.push_back(node);
    }
  }
}
