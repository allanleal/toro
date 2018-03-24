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

  params.addParam<std::vector<SubdomainName>>(
    "block", "The list of block ids (SubdomainID) that this object will be applied");

  params.addParam<std::vector<BoundaryName>>(
      "boundary", "The list of boundary IDs from the mesh where this boundary condition applies");

  ExecFlagEnum execute_options = MooseUtils::getDefaultExecFlagEnum();
  execute_options = EXEC_TIMESTEP_END;
  params.addParam<ExecFlagEnum>("execute_on", execute_options, execute_options.getDocString());

  return params;
}

AddSpeciesAction::AddSpeciesAction(const InputParameters & params)
  : Action(params)
{
}

void
AddSpeciesAction::act()
{
  if (_current_task == "add_reaktoro_aux_kernels")
  {
    std::cout<< "Adding Kernels!" <<std::endl;

    InputParameters params = _factory.getValidParams("SpeciesAux");

    if (isParamValid("block"))
      params.set<std::vector<SubdomainName>>("block") = getParam<std::vector<SubdomainName>>("block");
    else if (isParamValid("boundary"))
      params.set<std::vector<BoundaryName>>("boundary") = getParam<std::vector<BoundaryName>>("boundary");

    params.set<std::vector<VariableName>>("temperature") = getParam<std::vector<VariableName>>("temperature");
    params.set<std::vector<VariableName>>("pressure") = getParam<std::vector<VariableName>>("pressure");

    params.set<AuxVariableName>("variable") = _species_names[0];
    params.set<std::vector<VariableName>>("species") = _species_names;

    params.set<std::vector<VariableName>>("elements") = _element_names;

    params.set<std::vector<VariableName>>("nonlinear_elements") = _nonlinear_element_names;

    params.set<UserObjectName>("reaktoro_problem") = _name + "_reaktoro_problem";

    params.set<ExecFlagEnum>("execute_on") = getParam<ExecFlagEnum>("execute_on");

    _problem->addAuxKernel("SpeciesAux", _name + "_species_aux", params);
  }
  else if (_current_task == "add_reaktoro_aux_variables")
  {
    std::cout<< "Adding Aux Variables!" <<std::endl;

    const auto & reaktoro_problem = _problem->getUserObject<ReaktoroProblemUserObject>(_name + "_reaktoro_problem");

    auto elements = reaktoro_problem.getElementNames();

    for (const auto & element : elements)
    {
      std::cout << element << std::endl;

      // Add the variable
      _problem->addAuxVariable("rt_" + element,
                               FEType(Utility::string_to_enum<Order>(getParam<MooseEnum>("order")),
                                      Utility::string_to_enum<FEFamily>(getParam<MooseEnum>("family"))));

      _problem->addVariable(element,
                            FEType(Utility::string_to_enum<Order>(getParam<MooseEnum>("order")),
                                   Utility::string_to_enum<FEFamily>(getParam<MooseEnum>("family"))), 1.0);

      _nonlinear_element_names.push_back(element);
      _element_names.push_back("rt_" + element);
    }

    auto species = reaktoro_problem.getSpeciesNames();

    for (const auto & species : species)
    {
      std::cout << species << std::endl;

      // Add the variable
      _problem->addAuxVariable("rt_" + species,
                               FEType(Utility::string_to_enum<Order>(getParam<MooseEnum>("order")),
                                      Utility::string_to_enum<FEFamily>(getParam<MooseEnum>("family"))));

      _species_names.push_back("rt_" + species);
    }
  }
  else if (_current_task == "add_user_object")
  {
    std::cout<< "Adding UserObject!" <<std::endl;

    InputParameters params = _factory.getValidParams("ReaktoroProblemUserObject");

    params.set<std::vector<VariableName>>("species") = _species_names;

    params.applyParameters(_pars);

    _problem->addUserObject("ReaktoroProblemUserObject", _name + "_reaktoro_problem", params);
  }
}
