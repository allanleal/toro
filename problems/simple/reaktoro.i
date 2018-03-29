[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 20
  ny = 20
  nz = 20
[]

[Variables]
  [./pressure]
    initial_condition = 1e5
  [../]
  [./temperature]
    initial_condition = 300
  [../]
[]

[Kernels]
  [./pressure_storage]
    type =TimeDerivative
    variable = pressure
  [../]
  [./DarcyFlow]
    type = toroHydroDarcy
    variable = pressure
  [../]
  [./temperature_storage]
    type =TimeDerivative
    variable = temperature
  [../]
  [./TemperatureConduction]
    type = toroHeatConduction
    variable = temperature
  [../]
[]

[Materials]
  [./material_darcy]
    type = toroHydroConstant
    block = 0
   fluid_modulus = 1e9
   fluid_viscosity = 1e-3
   porosity = 0.2
   permeability = 1e-15
  [../]
  [./temperature_material]
    type = toroThermalConstant
    block = 0
    solid_heat_capacity = 2000
    fluid_heat_capacity = 1000
    solid_thermal_conductivity = 2.4
    fluid_thermal_conductivity = 0.65
    temperature = temperature
  [../]
  [./density_materia]
    type = toroDensityConstant
    block = 0
    solid_density = 2000
    fluid_density = 1000
  [../]
[]

[BCs]
  [./leftP]
    type = DirichletBC
    variable = pressure
    boundary = 'left'
    value = 2e5
  [../]
  [./leftT]
    type = DirichletBC
    variable = temperature
    boundary = 'left'
    value = 400
  [../]
  [./rightP]
    type = DirichletBC
    variable = pressure
    boundary = 'right'
    value = 1e5
  [../]
  [./rightT]
    type = DirichletBC
    variable = temperature
    boundary = 'right'
    value = 300
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 10
  dt = 10
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
[]

[Modules]
  [Reaktoro]
    # temperature = 'u'
    family = LAGRANGE
    order = FIRST
    substance_names = 'H2O NaCl'
    substance_amounts = '1 0.1'
    substance_units = 'kg mol'
    temperature = 'temperature'
    pressure = 'pressure'
  []
[]
