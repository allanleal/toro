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
#include "ReaktoroProblemUserObject.h"

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
  params.addParam<Real>("scaling", 1.0, "Specifies a scaling factor to apply to this variable");

  params.addParam<std::string>("pressure_units", "Pa", "Units pressure is in (Pa)");
  params.addParam<std::string>("temperature_units", "K", "Units temperature is in (K)");

  params.addParam<std::vector<VariableName>>("temperature", "The temperature.");
  params.addParam<std::vector<VariableName>>("pressure", "The pressure.");

  params.addRequiredParam<std::vector<std::string>>("substance_names", "Names of the substances, these go with substance_amounts and substance_units");
  params.addRequiredParam<std::vector<Real>>("substance_amounts", "The amount of the substance_names");
  params.addRequiredParam<std::vector<std::string>>("substance_units", "The units to use for each amount (kg, mol)");

  return params;
}

AddSpeciesAction::AddSpeciesAction(const InputParameters & params)
  : Action(params)
{
}

void
AddSpeciesAction::act()
{
  std::cout << "AddSpeciesAction::act(): " << _current_task << std::endl;



  if (_current_task == "add_reaktoro_aux_kernels")
  {
    std::cout<< "Adding Kernels!" <<std::endl;

    InputParameters params = _factory.getValidParams("SpeciesAux");

    params.set<std::vector<VariableName>>("temperature") = getParam<std::vector<VariableName>>("temperature");
    params.set<std::vector<VariableName>>("pressure") = getParam<std::vector<VariableName>>("pressure");

    params.set<AuxVariableName>("variable") = _species_names[0];
    params.set<std::vector<VariableName>>("species") = _species_names;

    params.set<UserObjectName>("reaktoro_problem") = "reaktoro_problem";

    _problem->addAuxKernel("SpeciesAux", "species_aux", params);
  }
  else if (_current_task == "add_reaktoro_aux_variables")
  {
    std::cout<< "Adding Aux Variables!" <<std::endl;

    const auto & reaktoro_problem = _problem->getUserObject<ReaktoroProblemUserObject>("reaktoro_problem");

    auto species = reaktoro_problem.getSpeciesNames();

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
  else if (_current_task == "add_user_object")
  {
    std::cout<< "Adding UserObject!" <<std::endl;

    InputParameters params = _factory.getValidParams("ReaktoroProblemUserObject");

    params.set<std::vector<VariableName>>("species") = _species_names;

    params.applyParameters(_pars);

    _problem->addUserObject("ReaktoroProblemUserObject", "reaktoro_problem", params);
  }
}
