project('Elasticity', ['cpp'],
  default_options : ['cpp_std=c++17'])

cmake = import('cmake')
opt_var = cmake.subproject_options()
opt_var.add_cmake_defines({'EIGEN_BUILD_DOC': false, 'BUILD_TESTING': false,
                           'EIGEN_LEAVE_TEST_IN_ALL_TARGET': false})
opt_var.set_override_option('cpp_std', 'c++17')
eigen = cmake.subproject('Eigen', options: opt_var).dependency('eigen')

galerkin = subproject('Galerkin').get_variable('galerkin_dep')
smv = subproject('SmallVector').get_variable('smallvector_dep')
tetmesh = subproject('TetMesh').get_variable('tetmesh_dep')
symsparse = subproject('SymSparse').get_variable('symsparse_dep')

libelasticity = library('elasticity',
  files(['integrate_error.cpp',
         'integrate_volume_force_1.cpp',
         'integrate_volume_force_2.cpp',
         'integrate_volume_force_3.cpp',
         'integrate_volume_force_4.cpp',
         'assemble_stiffness_1_double.cpp',
         'assemble_stiffness_1_vector.cpp',
         'assemble_stiffness_2_double.cpp',
         'assemble_stiffness_2_vector.cpp',
         'assemble_stiffness_3_double.cpp',
         'assemble_stiffness_3_vector.cpp',
         'assemble_stiffness_4_double.cpp',
         'assemble_stiffness_4_vector.cpp',
         'impose_dirichlet_condition_1.cpp',
         'impose_dirichlet_condition_2.cpp',
         'impose_dirichlet_condition_3.cpp',
         'impose_dirichlet_condition_4.cpp',
         'stiffness_matrix_1.cpp',
         'stiffness_matrix_2.cpp',
         'stiffness_matrix_3.cpp',
         'stiffness_matrix_4.cpp',
         'mass_matrix_1.cpp',
         'mass_matrix_2.cpp',
         'mass_matrix_3.cpp',
         'mass_matrix_4.cpp',
         'integrate_traction_force_1.cpp',
         'integrate_traction_force_2.cpp',
         'integrate_traction_force_3.cpp',
         'integrate_traction_force_4.cpp',
         'add_traction_force.cpp',
         'read_mesh.cpp',
         'mesh_tools.cpp',
         'fill_internal_coordinates.cpp',
         'impose_homogeneous.cpp'
         ]),
  dependencies : [eigen, galerkin, tetmesh, symsparse])

elasticity_dep = declare_dependency(link_with : libelasticity,
  include_directories : include_directories('.'),
  dependencies : [eigen, galerkin, tetmesh, symsparse])

subdir('validation')

