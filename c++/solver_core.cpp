#include "./solver_core.hpp"
#include <fstream>

namespace alps_cthyb {

struct block_visitor {
  int num_flavors;
  void operator()(const indices_type &bl) { num_flavors += bl.size(); }
};

struct index_visitor {
  std::vector<std::string> indices;
  void operator()(int i) { indices.push_back(std::to_string(i)); }
  void operator()(std::string s) { indices.push_back(s); }
};

solver_core::solver_core(double beta_,
                         std::map<std::string, indices_type> const &gf_struct_,
                         int n_iw,
                         int n_tau,
                         int n_l) :
    beta(beta_), gf_struct(gf_struct_), num_flavors(0), n_iw_(n_iw), n_tau_(n_tau), n_l_(n_l) {

  if (n_tau < 2 * n_iw) {
    TRIQS_RUNTIME_ERROR << "Must use as least twice as many tau points as Matsubara frequencies: n_iw = " << n_iw
                        << " but n_tau = " << n_tau << ".";
  }

  std::vector<std::string> block_names;
  std::vector<gf<imfreq>> g0_iw_blocks;
  std::vector<gf<imtime>> g_tau_blocks;
  std::vector<gf<legendre>> g_l_blocks;
  std::vector<gf<imtime>> delta_tau_blocks;
  std::vector<gf<imtime, delta_target_t>>
      g_tau_accum_blocks; //  Local real or complex (if complex mode) quantities for accumulation

  for (auto const &bl : gf_struct) {
    block_names.push_back(bl.first);
    int n = bl.second.size();
    num_flavors += bl.second.size();

    index_visitor iv;
    for (auto &ind: bl.second) { apply_visitor(iv, ind); }
    std::vector<std::vector<std::string>> indices{{iv.indices, iv.indices}};

    g0_iw_blocks.push_back(gf<imfreq>{{beta, Fermion, n_iw}, {n, n}, indices});
    g_tau_blocks.push_back(gf<imtime>{{beta, Fermion, n_tau}, {n, n}, indices});
    g_l_blocks.push_back(gf<legendre>{{beta, Fermion, static_cast<size_t>(n_l)}, {n, n},
                                      indices}); // FIXME: cast is ugly
    delta_tau_blocks.push_back(gf<imtime>{{beta, Fermion, n_tau}, {n, n}, indices});
    g_tau_accum_blocks.push_back(gf<imtime, delta_target_t>{{beta, Fermion, n_tau}, {n, n}});
  }

  _G0_iw = make_block_gf(block_names, g0_iw_blocks);
  _G_tau = make_block_gf(block_names, g_tau_blocks);
  _G_l = make_block_gf(block_names, g_l_blocks);
  _Delta_tau = make_block_gf(block_names, delta_tau_blocks);
  _G_tau_accum = make_block_gf(block_names, g_tau_accum_blocks);

}

/// -------------------------------------------------------------------------------------------

//void solver_core::solve(solve_parameters_t const & params) {
//
//}

void solver_core::solve(solve_parameters_t const &params) {
  using variant_int_string = triqs::utility::variant_int_string;

  _last_solve_parameters = params;

  // determine basis of operators to use
  fundamental_operator_set fops;
  for (auto const &bl: gf_struct) {
    for (auto const &a: bl.second) {
      std::cout << "debug : " << bl.first << " : " << bl.second << " : " << a << std::endl;
      fops.insert(bl.first, a);
    }
  }

  std::cout << "num_flavors " << num_flavors << std::endl;

  // setup the linear index map
  std::map<std::pair<int, int>, int> linindex;
  std::map<std::pair<variant_int_string, variant_int_string>, int> linindex2;
  int block_index = 0;
  int flavor_index = 0;
  for (auto const &bl: gf_struct) {
    int inner_index = 0;
    for (auto const &a: bl.second) {
      linindex[std::make_pair(block_index, inner_index)] = fops[{bl.first, a}];
      linindex2[std::make_pair(variant_int_string(bl.first), a)] = flavor_index;
      ++ flavor_index;
      ++ inner_index;
    }
    ++ block_index;
  }

  // Make list of block sizes
  std::vector<int> n_inner;
  for (auto const &bl : gf_struct) {
    n_inner.push_back(bl.second.size());
  }

  // Calculate imfreq quantities
  auto G0_iw_inv = map([](gf_const_view<imfreq> x) { return triqs::gfs::inverse(x); }, _G0_iw);
  auto Delta_iw = G0_iw_inv;

  //Compute Coulomb tensor
  std::cout << "h_loc " << params.h_int << std::endl;
  boost::multi_array<std::complex<double>, 4> Uijkl(boost::extents[num_flavors][num_flavors][num_flavors][num_flavors]);
  std::fill(Uijkl.origin(), Uijkl.origin() + Uijkl.num_elements(), 0.0);
  for (auto it = _h_loc.cbegin(); it != _h_loc.cend(); ++it) {
    if (it->coef == 0.0) {
      continue;
    }
    auto num_ops = std::distance(it->monomial.cbegin(), it->monomial.cend());
    if (num_ops != 4) {
      std::stringstream ss;
      ss << "Unsupported interaction term: " << it->monomial;
      TRIQS_RUNTIME_ERROR << ss.str();
    }
    auto iop = 0;
    boost::array<int,4> indices;
    for (auto it_op = it->monomial.cbegin(); it_op != it->monomial.cend(); ++it_op) {
      if ((iop < num_ops/2 && !it_op->dagger)
          || (iop >= num_ops/2 && it_op->dagger)
          || std::distance(it_op->indices.cbegin(), it_op->indices.cend()) != 2) {
        std::stringstream ss;
        ss << "Unsupported interaction term: " << it->monomial;
        TRIQS_RUNTIME_ERROR << ss.str();
      }
      indices[iop] = linindex2[std::make_pair(it_op->indices[0], it_op->indices[1])];
      ++ iop;
    }
    Uijkl[indices[0]][indices[1]][indices[2]][indices[3]] = it->coef;
  }
  {
    std::vector<std::complex<double> > Uijkl_vec(Uijkl.num_elements());
    std::copy(Uijkl.origin(), Uijkl.origin() + Uijkl.num_elements(), Uijkl_vec.begin());
  }

  // Do I have imaginary components in my local Hamiltonian?
  auto max_imag = 0.0;
  for (int b : range(gf_struct.size())) {
    max_imag = std::max(max_imag, max_element(abs(real(_G0_iw[b].singularity()(2)))));
  }

  // Add quadratic terms to h_loc
  std::vector<h_scalar_t> h_loc_vec(num_flavors * num_flavors);
  int b = 0;
  int offset = 0;
  for (auto const &bl: gf_struct) {

    int n1 = 0;
    for (auto const &a1: bl.second) {
      int n2 = 0;
      for (auto const &a2 : bl.second) {
        dcomplex e_ij;
        if (max_imag > params.imag_threshold) {
          e_ij = _G0_iw[b].singularity()(2)(n1, n2);
        } else {
          e_ij = _G0_iw[b].singularity()(2)(n1, n2).real();
        }
        _h_loc = _h_loc + e_ij * c_dag<h_scalar_t>(bl.first, a1) * c<h_scalar_t>(bl.first, a2);
        //h_loc_vec[(n1 + offset) * num_flavors + (n2 + offset)] = ;
        std::cout << linindex[std::make_pair(b, n1)] << "  "  << linindex[std::make_pair(b,n2)] << " " << e_ij << std::endl;
        n2++;
      }
      n1++;
    }
    b++;
    offset += bl.second.size();
  }

  // Determine terms Delta_iw from G0_iw and ensure that the 1/iw behaviour of G0_iw is correct
  b = 0;
  range _;
  triqs::clef::placeholder<0> iw_;
  for (auto const &bl : gf_struct) {
    Delta_iw[b](iw_) << G0_iw_inv[b].singularity()(-1) * iw_ + G0_iw_inv[b].singularity()(0);
    Delta_iw[b] = Delta_iw[b] - G0_iw_inv[b];
    _Delta_tau[b]() = inverse_fourier(Delta_iw[b]);
    // Force all diagonal elements to be real
    for (int i : range(bl.second.size())) _Delta_tau[b].data()(_, i, i) = real(_Delta_tau[b].data()(_, i, i));
    // If off-diagonal elements are below threshold, set to real
    if (max_element(abs(imag(_Delta_tau[b].data()))) < params.imag_threshold)
      _Delta_tau[b].data() = real(_Delta_tau[b].data());
    b++;
  }

  // Determine which solver to be used the real-number solver or the complex-number solver
  alps::params par;
  boost::shared_ptr<alps::cthyb::Solver> p_solver;
  alps::cthyb::MatrixSolver<std::complex<double> >::define_parameters(par);
  //if (par.help_requested(std::cout)) { exit(0); } //If help message is requested, print it and exit normally.

  // Set parameters for ALPS CT-HYB solver
  par["timelimit"] = params.max_time > 0 ? params.max_time : 1E+30;
  par["verbose"] = params.verbosity == 0 ? 0 : 1;

  par["model.sites"] = num_flavors;
  par["model.spins"] = 1;
  par["model.beta"] = beta;
  par["model.n_tau_hyb"] = n_tau_;
  par["model.coulomb_tensor"] = 0;
  par["model.hopping_matrix"] = 0;
  par["model.delta"] = 0;

  par["measurement.G1.n_legendre"] = n_l_;
  par["measurement.G1.n_tau"] = n_tau_;
  par["measurement.G1.n_iw"] = n_iw_;

  //int n_iw,
  //int n_tau,
  //int n_l) :

  p_solver.reset(new alps::cthyb::MatrixSolver<std::complex<double> >(par));

  // Call the ALPS CT-HYB solver

  //if (params.verbosity >= 2) std::cout << "Average sign: " << _average_sign << std::endl;

  // Copy local (real or complex) G_tau back to complex G_tau
  //if (params.measure_g_tau) _G_tau = _G_tau_accum;
}

}
