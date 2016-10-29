# Generated automatically using the command :
# c++2py.py ../c++/solver_core.hpp -p -mpytriqs.applications.impurity_solvers.alps_cthyb -o alps_cthyb --moduledoc "The ALPS cthyb solver" -I /opt/local/include/openmpi-clang38
from wrap_generator import *

# The module
module = module_(full_name = "pytriqs.applications.impurity_solvers.alps_cthyb", doc = "The ALPS cthyb solver", app_name = "pytriqs.applications.impurity_solvers.alps_cthyb")

# All the triqs C++/Python modules
module.use_module('gf', 'triqs')
module.use_module('operators', 'triqs')

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("solver_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/python_tools/converters/pair.hpp>
#include <triqs/python_tools/converters/map.hpp>
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/variant.hpp>
#include <triqs/python_tools/converters/arrays.hpp>
using namespace triqs::gfs;
using triqs::operators::many_body_operator;
using namespace alps_cthyb;
#include "./pytriqs.applications.impurity_solvers.alps_cthyb_converters.hxx"
""")

# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "solver_core",   # name of the C++ class
        doc = r"DOC OF SOLVER CORE",   # doc of the C++ class
)

c.add_constructor("""(double beta, std::map<std::string,indices_type> gf_struct, int n_iw = 1025, int n_tau = 10001, int n_l = 50)""",
                  doc = """ """)

c.add_method("""void solve (**alps_cthyb::solve_parameters_t)""",
             doc = """+----------------+----------+---------------------------+--------------------------------------------------------------------------------+
| Parameter Name | Type     | Default                   | Documentation                                                                  |
+================+==========+===========================+================================================================================+
| h_int          | Operator |                           | Interacting part of the atomic Hamiltonian                                     |
+----------------+----------+---------------------------+--------------------------------------------------------------------------------+
| random_seed    | int      | 34788 + 928374 * MPI.rank | Seed for random number generator                                               |
+----------------+----------+---------------------------+--------------------------------------------------------------------------------+
| max_time       | int      | -1 = infinite             | Maximum runtime in seconds, use -1 to set infinite                             |
+----------------+----------+---------------------------+--------------------------------------------------------------------------------+
| verbosity      | int      | 0                         |                                                                                |
+----------------+----------+---------------------------+--------------------------------------------------------------------------------+
| imag_threshold | double   | 1.e-15                    | Threshold below which imaginary components of Delta and h_loc are set to zero  |
+----------------+----------+---------------------------+--------------------------------------------------------------------------------+ """)

c.add_property(name = "h_loc",
               getter = cfunction("many_body_op_t h_loc ()"),
               doc = """The local Hamiltonian of the problem : H_loc used in the last call to solve. """)

c.add_property(name = "last_solve_parameters",
               getter = cfunction("alps_cthyb::solve_parameters_t last_solve_parameters ()"),
               doc = """Set of parameters used in the last call to solve """)

c.add_property(name = "G0_iw",
               getter = cfunction("block_gf_view<imfreq> G0_iw ()"),
               doc = """G0(iw) in imaginary frequencies """)

c.add_property(name = "Delta_tau",
               getter = cfunction("block_gf_view<imtime> Delta_tau ()"),
               doc = """Delta(tau) in imaginary time """)

c.add_property(name = "G_tau",
               getter = cfunction("block_gf_view<imtime> G_tau ()"),
               doc = """G(tau) in imaginary time """)

c.add_property(name = "G_l",
               getter = cfunction("block_gf_view<legendre> G_l ()"),
               doc = """G_l in Legendre polynomials representation """)

c.add_property(name = "density_matrix",
               getter = cfunction("std::vector<matrix_t> density_matrix ()"),
               doc = """Density matrix """)

c.add_property(name = "average_sign",
               getter = cfunction("mc_weight_t average_sign ()"),
               doc = """Monte Carlo average sign """)

c.add_property(name = "solve_status",
               getter = cfunction("int solve_status ()"),
               doc = """Status of the solve on exit """)

module.add_class(c)

module.generate_code()