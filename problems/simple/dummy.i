[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 5
  ny = 3
  nz = 2
[]

[Variables]
  [./u]
  [../]
[]

[AuxVariables]
  [./from_sub]
  [../]
[]

[Functions]
  [./bc_func]
    type = ParsedFunction
    value = y+1
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
  [./time]
    type = TimeDerivative
    variable = u
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  [../]
  [./right]
    type = FunctionDirichletBC
    variable = u
    boundary = right
    function = bc_func
  [../]
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