[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 5
  ny = 3
  nz = 2
[]

[Variables]
  [temp]
    initial_condition = 300
  []
[]

[AuxVariables]
  [pressure]
    initial_condition = 5e6
  []
[]

[Functions]
  [bc_func]
    type = ParsedFunction
    value = '300+(100*y)'
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = temp
  []
  [time]
    type = TimeDerivative
    variable = temp
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = temp
    boundary = 'left'
    value = 300
  []
  [right]
    type = FunctionDirichletBC
    variable = temp
    boundary = 'right'
    function = bc_func
  []
[]

[Executioner]
  type = Transient
  num_steps = 5
  dt = 0.1
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
[]

[Modules]
  [Reaktoro]
    [Problems]
      # active = right
      [bulk]
        family = LAGRANGE
        order = FIRST
        substance_names = 'H2O NaCl'
        substance_amounts = '1 0.1'
        substance_units = 'kg mol'
        temperature = 'temp'
        pressure = 'pressure'
      []
      [right]
        family = LAGRANGE
        order = FIRST
        substance_names = 'H2O'
        substance_amounts = '1.00000001'
        substance_units = 'kg'
        temperature = 'temp'
        pressure = 'pressure'
        boundary = 'left'
      []
    []
  []
[]