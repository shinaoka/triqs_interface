#pragma once
#include "./config.hpp"

namespace alps_cthyb {

using namespace triqs::operators;
using indices_map_t = std::map<triqs::operators::indices_t,triqs::operators::indices_t>;

// All the arguments of the solve function
struct solve_parameters_t {

 /// If assume_real == true, the real-number solver will be used. Otherwise, the complex version will called.
 //bool assume_real;
 
 /// Interacting part of the atomic Hamiltonian
 /// type: Operator
 many_body_op_t h_int;

 /// Number of QMC cycles
 //int n_cycles;

 /// Partition method
 /// type: str
 //std::string partition_method = "autopartition";

 /// Quantum numbers
 /// type: list(Operator)
 /// default: []
 //std::vector<many_body_op_t> quantum_numbers = std::vector<many_body_op_t>{};

 /// Length of a single QMC cycle
 /// default: 50
 //int length_cycle = 50;

 /// Number of cycles for thermalization
 /// default: 5000
 //int n_warmup_cycles = 5000;

 /// Seed for random number generator
 /// default: 34788 + 928374 * MPI.rank
 int random_seed = 34788 + 928374 * triqs::mpi::communicator().rank();

 /// Name of random number generator
 /// type: str
 //std::string random_name = "";

 /// Maximum runtime in seconds, Maximum runtime in seconds, use -1 to use a default value
 int max_time = -1;

 /// Thermalization time in seconds
 /// default: -1 = 10 % of total simulation time
 int thermalization_time = -1;

 /// Verbosity level
 /// default: 3 on MPI rank 0, 0 otherwise.
 //int verbosity = ((triqs::mpi::communicator().rank() == 0) ? 3 : 0); // silence the slave nodes
 int verbosity = 0;

 /// Add shifting a move as a move?
 //bool move_shift = true;

 /// Add double insertions as a move?
 //bool move_double = false;

 /// Calculate the full trace or use an estimate?
 //bool use_trace_estimator = false;

 /// Measure G(tau)?
 //bool measure_g_tau = true;

 /// Measure G_l (Legendre)?
 //bool measure_g_l = false;

 /// Measure perturbation order?
 //bool measure_pert_order = false;

 /// Measure the contribution of each atomic state to the trace?
 //bool measure_density_matrix= false;

 /// Use the norm of the density matrix in the weight if true, otherwise use Trace
 //bool use_norm_as_weight = false;

 /// Analyse performance of trace computation with histograms (developers only)?
 //bool performance_analysis = false;

 /// Operator insertion/removal probabilities for different blocks
 /// type: dict(str:float)
 /// default: {}
 //std::map<std::string,double> proposal_prob = (std::map<std::string,double>{});

 /// List of global moves (with their names).
 /// Each move is specified with an index substitution dictionary
 /// type: dict(str : dict(indices : indices))
 /// default: {}
 //std::map<std::string,indices_map_t> move_global = (std::map<std::string,indices_map_t>{});

 /// Overall probability of the global moves
 //double move_global_prob = 0.05;

 /// Threshold below which imaginary components of Delta and h_loc are set to zero
 double imag_threshold = 1.e-15;

 /// 0 : no rotation, 1 : diagonalize the local transfer matrix
 int basis_rotation = 0;

 /// Dump parameters passed to ALPSCore/CT-HYB into a file
 //std::string params_dump_file = "";

 solve_parameters_t() {}

 //solve_parameters_t(many_body_op_t h_int, int n_cycles) : h_int(h_int), n_cycles(n_cycles) {}

};
}
