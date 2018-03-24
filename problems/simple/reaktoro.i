[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 5
  ny = 3
  nz = 2
[]

[Variables]
  [u]
    initial_condition = 300
  []
[]

[AuxVariables]
  [from_sub]
    initial_condition = 300
  []
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
    variable = u
  []
  [time]
    type = TimeDerivative
    variable = u
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = 'left'
    value = 300
  []
  [right]
    type = FunctionDirichletBC
    variable = u
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
    [./Problems]
      [./interior]
        family = LAGRANGE
	order = FIRST

        substance_names = 'H2O NaCl'
        substance_amounts = '1 0.1'
        substance_units = 'kg mol'

        temperature = 'u'
        pressure = 'pressure'
      []
    []
  []
[]
