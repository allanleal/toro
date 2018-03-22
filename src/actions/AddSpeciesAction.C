//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AddSpeciesAction.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "Factory.h"
#include "FEProblem.h"
#include "ReaktoroMultiApp.h"

#include "libmesh/string_to_enum.h"

using namespace Reaktoro;

template <>
InputParameters
validParams<AddSpeciesAction>()
{
  InputParameters params = validParams<Action>();
  params.addClassDescription("Set up order parameter variables for a polycrystal simulation");
  // Get MooseEnums for the possible order/family options for this variable
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());
  params.addParam<MooseEnum>("family",
                             families,
                             "Specifies the family of FE "
                             "shape function to use for the order parameters");
  params.addParam<MooseEnum>("order",
                             orders,
                             "Specifies the order of the FE "
                             "shape function to use for the order parameters");
  params.addRequiredParam<MultiAppName>("multiapp", "The name of the multiapp");
  params.addParam<Real>("scaling", 1.0, "Specifies a scaling factor to apply to this variable");
  params.addRequiredParam<std::string>("species_prefix", "specifies the base name of the variables");
  return params;
}

AddSpeciesAction::AddSpeciesAction(const InputParameters & params)
  : Action(params),
    _multiapp_name(getParam<MultiAppName>("multiapp")),
    _species_prefix(getParam<std::string>("species_prefix"))
{
}

void
AddSpeciesAction::act()
{
  if (_current_task == "add_aux_kernel")
  {
    InputParameters params = _factory.getValidParams("SpeciesAux");
    params.set<AuxVariableName>("variable") = _species_names[0];
    params.set<std::vector<VariableName>>("species") = _species_names;
    params.set<Reaktoro::ChemicalSystem *>("reaktoro_system") = &_system;

    _problem->addAuxKernel("SpeciesAux", "species_aux", params);
  }
  else if (_current_task == "add_aux_variable")
  {
    createChemicalSystem();

    auto species = names(_system.species());
    for (const auto & species : species)
    {
      std::cout << species << std::endl;

      // Add the variable
      _problem->addAuxVariable(species,
                               FEType(Utility::string_to_enum<Order>(getParam<MooseEnum>("order")),
                                      Utility::string_to_enum<FEFamily>(getParam<MooseEnum>("family"))));

      _species_names.push_back(species);
    }
  }
}

void
AddSpeciesAction::createChemicalSystem()
{
  Database database("supcrt98.xml");

  ChemicalEditor editor(database);
  editor.addAqueousPhase({"H2O(l)", "H+", "OH-", "Na+", "Cl-"});

  _system = ChemicalSystem(editor);



//  _state_bc.output("state_bc.txt");
//  _state_ic.output("state_ic.txt");
}
