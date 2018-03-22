#include "ToroApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

template <>
InputParameters
validParams<ToroApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

ToroApp::ToroApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  ToroApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  ToroApp::associateSyntax(_syntax, _action_factory);

  Moose::registerExecFlags(_factory);
  ModulesApp::registerExecFlags(_factory);
  ToroApp::registerExecFlags(_factory);
}

ToroApp::~ToroApp() {}

void
ToroApp::registerApps()
{
  registerApp(ToroApp);
}

void
ToroApp::registerObjects(Factory & factory)
{
    Registry::registerObjectsTo(factory, {"ToroApp"});
}

void
ToroApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & action_factory)
{
  Registry::registerActionsTo(action_factory, {"ToroApp"});

  /* Uncomment Syntax parameter and register your new production objects here! */
}

void
ToroApp::registerObjectDepends(Factory & /*factory*/)
{
}

void
ToroApp::associateSyntaxDepends(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}

void
ToroApp::registerExecFlags(Factory & /*factory*/)
{
  /* Uncomment Factory parameter and register your new execution flags here! */
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
ToroApp__registerApps()
{
  ToroApp::registerApps();
}

extern "C" void
ToroApp__registerObjects(Factory & factory)
{
  ToroApp::registerObjects(factory);
}

extern "C" void
ToroApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  ToroApp::associateSyntax(syntax, action_factory);
}

extern "C" void
ToroApp__registerExecFlags(Factory & factory)
{
  ToroApp::registerExecFlags(factory);
}
