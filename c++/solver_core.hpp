/*******************************************************************************
 *
 * C++ interface to TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016, H. Shinaoka, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
 *
 * This file is licenced under the GNU General Public License, version 3.
 *
 ******************************************************************************/
#pragma once
#include <complex>

#include <triqs/mc_tools.hpp>
#include <triqs/utility/callbacks.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/statistics/histograms.hpp>
#include <triqs/gfs/functions/functions.hpp>
#include <triqs/arrays/linalg/eigenelements.hpp>

#include <alps/params.hpp>

#include <alps/cthyb/solver.hpp>

#include <boost/multi_array.hpp>

#include "solve_parameters.hpp"
#include "config.hpp"
//#include "atom_diag.hpp"
//#include "atom_diag_functions.hpp"

namespace alps_cthyb {

using namespace triqs::utility;
using namespace triqs::statistics;
using namespace triqs::gfs;
using histo_map_t = std::map<std::string, histogram>;
using indices_type = triqs::operators::indices_t;

/**  DOC OF SOLVER CORE*/
class solver_core {
 int num_flavors;
 int n_iw_;
 int n_tau_;
 int n_l_;
 double beta;                                   // inverse temperature
 std::map<std::string, indices_type> gf_struct; // Block structure of the Green function FIXME
 bool assume_real_;                             // Assume real Hamiltonian
 many_body_op_t _h_loc;                         // The local Hamiltonian = h_int + h0
 block_gf<imfreq> _G0_iw;                       // Green's function containers: imaginary-freq Green's functions
 block_gf<imtime> _Delta_tau, _G_tau;           // Green's function containers: imaginary-time Green's functions
 block_gf<imtime, g_target_t> _G_tau_accum;     // Intermediate object to accumulate g(tau), either real or complex
 block_gf<legendre> _G_l;                       // Green's function containers: Legendre coefficients
 boost::multi_array<double,3> delta_tau_Re_, delta_tau_Im_;
 histogram _pert_order_total;                   // Histogram of the total perturbation order
 histo_map_t _pert_order;                       // Histograms of the perturbation order for each block
 std::vector<matrix_t> _density_matrix;         // density matrix, when used in Norm mode
 triqs::mpi::communicator _comm;                // define the communicator, here MPI_COMM_WORLD
 solve_parameters_t _last_solve_parameters;     // parameters of the last call to solve
 histo_map_t _performance_analysis;             // Histograms used for performance analysis
 mc_weight_t _average_sign;                     // average sign of the QMC
 int _solve_status;                             // Status of the solve upon exit: 0 for clean termination, > 0 otherwise.

 public:
 solver_core(double beta, std::map<std::string, indices_type> const & gf_struct, bool assume_real, int n_iw=1025, int n_tau=10001, int n_l=50);

 /// Solve the impurity problem for the given Hamiltonian h_loc and with specified parameters params.
 TRIQS_WRAP_ARG_AS_DICT // Wrap the solver parameters as a dictionary in python with the clang tool
 void solve(solve_parameters_t const & p);

 /// The local Hamiltonian of the problem : H_loc used in the last call to solve.
 many_body_op_t const & h_loc() const { return _h_loc; }

 /// Set of parameters used in the last call to solve
 solve_parameters_t last_solve_parameters() const {return _last_solve_parameters;}

 /// G0(iw) in imaginary frequencies
 block_gf_view<imfreq> G0_iw() { return _G0_iw; }

 /// Delta(tau) in imaginary time
 block_gf_view<imtime> Delta_tau() { return _Delta_tau; }

 /// G(tau) in imaginary time
 block_gf_view<imtime> G_tau() { return _G_tau; }

 /// G_l in Legendre polynomials representation
 block_gf_view<legendre> G_l() { return _G_l; }

 /// Atomic G(tau) in imaginary time
 //block_gf_view<imtime> atomic_gf() const { return ::cthyb::atomic_gf(h_diag, beta, gf_struct, _G_tau[0].mesh().size()); }

 /// Density matrix
 std::vector<matrix_t> const & density_matrix() const { return _density_matrix;}

 /// Diagonalization of h_loc
 //atom_diag const & h_loc_diagonalization() const { return h_diag;}

 /// Histogram of the total perturbation order
 //histogram const& get_perturbation_order_total() const { return _pert_order_total; }

 /// Histograms of the perturbation order for each block
 //histo_map_t const& get_perturbation_order() const { return _pert_order; }

 /// Histograms related to the performance analysis
 //histo_map_t const& get_performance_analysis() const { return _performance_analysis; }

 /// Monte Carlo average sign
 mc_weight_t average_sign() const { return _average_sign; }

 /// Status of the solve on exit
 int solve_status() const { return _solve_status; }

};

namespace detail {
  template<typename GA, typename GT>
  void copy_from_alps_to_triqs_gf(const GA &ga, GT &gt, bool assume_real=false) {
    const auto num_blocks = gt.data().size();
    const auto n_tau = gt.data()[0].data().shape()[0];

    //Count the number of flavors
    int num_flavors = 0;
    for (int b : range(num_blocks)) {
      num_flavors += gt.data()[b].data().shape()[1];
    }

    int offset = 0;
    for (int b : range(num_blocks)) {
      const auto num_flavors_block = gt.data()[b].data().shape()[1];
      for (auto itau = 0; itau < n_tau; ++itau) {
        for (auto f1 = 0; f1 < num_flavors_block; ++f1) {
          for (auto f2 = 0; f2 < num_flavors_block; ++f2) {
            if (assume_real) {
              gt.data()[b].data()(itau,f1,f2) =
                  ga(alps::gf::itime_index(itau), alps::gf::index(f1+offset), alps::gf::index(f2+offset)).real();
            } else {
              gt.data()[b].data()(itau,f1,f2) =
                  ga(alps::gf::itime_index(itau), alps::gf::index(f1+offset), alps::gf::index(f2+offset));
            }
          }
        }
      }
      offset += num_flavors_block;
      gt[b].singularity()(1) = 1.0;
    }
  }

  template<typename GA, typename GT>
  void copy_Gl(const GA &ga, GT &gt, bool assume_real=false) {
    const auto num_blocks = gt.data().size();
    const auto n_tau = gt.data()[0].data().shape()[0];

    //Count the number of flavors
    int num_flavors = 0;
    for (int b : range(num_blocks)) {
      num_flavors += gt.data()[b].data().shape()[1];
    }

    int offset = 0;
    for (int b : range(num_blocks)) {
      const auto num_flavors_block = gt.data()[b].data().shape()[1];
      for (auto il = 0; il < n_tau; ++il) {
        for (auto f1 = 0; f1 < num_flavors_block; ++f1) {
          for (auto f2 = 0; f2 < num_flavors_block; ++f2) {
            if (assume_real) {
              gt.data()[b].data()(il,f1,f2) = ga[f1+offset][f2+offset][il].real();
            } else {
              gt.data()[b].data()(il,f1,f2) = ga[f1+offset][f2+offset][il];
            }
          }
        }
      }
      offset += num_flavors_block;
      //gt[b].singularity()(1) = 1.0;
      triqs::arrays::matrix<double> id(get_target_shape(gt[b]));
      id() = 1.0;
      auto v = triqs::gfs::gf_view<triqs::gfs::legendre>(gt[b]);
      triqs::gfs::enforce_discontinuity(v, id);
    }
  }

  template <typename... T>
  void mpi_broadcast(block_gf<T...> &g, triqs::mpi::communicator c = {}, int root = 0) {
    const auto num_blocks = g.data().size();
    for (int b : range(num_blocks)) {
      triqs::arrays::mpi_broadcast(g[b].data(), c, root);
      triqs::arrays::mpi_broadcast(g[b].singularity().data(), c, root);
    }
  }
}

}
