//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "ToroTestApp.h"
#include "ToroApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

template <>
InputParameters
validParams<ToroTestApp>()
{
  InputParameters params = validParams<ToroApp>();
  return params;
}

ToroTestApp::ToroTestApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  ToroApp::registerObjectDepends(_factory);
  ToroApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  ToroApp::associateSyntaxDepends(_syntax, _action_factory);
  ToroApp::associateSyntax(_syntax, _action_factory);

  Moose::registerExecFlags(_factory);
  ModulesApp::registerExecFlags(_factory);
  ToroApp::registerExecFlags(_factory);

  bool use_test_objs = getParam<bool>("allow_test_objects");
  if (use_test_objs)
  {
    ToroTestApp::registerObjects(_factory);
    ToroTestApp::associateSyntax(_syntax, _action_factory);
    ToroTestApp::registerExecFlags(_factory);
  }
}

ToroTestApp::~ToroTestApp() {}

void
ToroTestApp::registerApps()
{
  registerApp(ToroApp);
  registerApp(ToroTestApp);
}

void
ToroTestApp::registerObjects(Factory & /*factory*/)
{
  /* Uncomment Factory parameter and register your new test objects here! */
}

void
ToroTestApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
  /* Uncomment Syntax and ActionFactory parameters and register your new test objects here! */
}

void
ToroTestApp::registerExecFlags(Factory & /*factory*/)
{
  /* Uncomment Factory parameter and register your new execute flags here! */
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
ToroTestApp__registerApps()
{
  ToroTestApp::registerApps();
}

// External entry point for dynamic object registration
extern "C" void
ToroTestApp__registerObjects(Factory & factory)
{
  ToroTestApp::registerObjects(factory);
}

// External entry point for dynamic syntax association
extern "C" void
ToroTestApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  ToroTestApp::associateSyntax(syntax, action_factory);
}

// External entry point for dynamic execute flag loading
extern "C" void
ToroTestApp__registerExecFlags(Factory & factory)
{
  ToroTestApp::registerExecFlags(factory);
}
