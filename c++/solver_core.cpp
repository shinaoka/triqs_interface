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
    std::cout << "block_name " << bl.first << std::endl;
    block_names.push_back(bl.first);
    std::cout << "block_name " << bl.first << std::endl;
    int n = bl.second.size();
    std::cout << "block_name " << bl.first << std::endl;
    num_flavors += bl.second.size();

    std::cout << "A " << bl.first << std::endl;
    index_visitor iv;
    std::cout << "B " << bl.first << std::endl;
    for (auto &ind: bl.second) { apply_visitor(iv, ind); }
    std::cout << "C " << bl.first << std::endl;
    std::vector<std::vector<std::string>> indices{{iv.indices, iv.indices}};
    std::cout << "D " << bl.first << std::endl;

    g0_iw_blocks.push_back(gf<imfreq>{{beta, Fermion, n_iw}, {n, n}, indices});
    g_tau_blocks.push_back(gf<imtime>{{beta, Fermion, n_tau}, {n, n}, indices});
    g_l_blocks.push_back(gf<legendre>{{beta, Fermion, static_cast<size_t>(n_l)}, {n, n}, indices});
    delta_tau_blocks.push_back(gf<imtime>{{beta, Fermion, n_tau}, {n, n}, indices});
    g_tau_accum_blocks.push_back(gf<imtime, delta_target_t>{{beta, Fermion, n_tau}, {n, n}});
    std::cout << "E " << bl.first << std::endl;
  }

  auto tmp = delta_tau_blocks[0].indices()[0];
  auto tmp2 = delta_tau_blocks[0].indices().ind_vec;
  for (auto it = tmp2.begin(); it != tmp2.end(); ++it) {
    int i = 0;
    for (auto it2 = it->begin(); it2 != it->end(); ++it2) {
      std::cout << "debug3   " << i << " " << *it2 << std::endl;
      ++ i;
    }
  }


  //n_tau_hyb_ = delta_tau_blocks[0].data().shape(0);
  //delta_tau_Re_.resize(boost::extents[n_tau_hyb_][num_flavors][num_flavors]);
  //delta_tau_Im_.resize(boost::extents[n_tau_hyb_][num_flavors][num_flavors]);
  //for (auto delta_bl = delta_tau_blocks.cbegin(); delta_bl != delta_tau_blocks.cend(); ++ delta_bl) {
  //}

  auto tmp3 = delta_tau_blocks[0];
  std::cout << "rank " << tmp3.data().num_elements() << std::endl;
  std::cout << "shape " << tmp3.data().shape(0) << std::endl;
  std::cout << "shape " << tmp3.data().shape(1) << std::endl;
  std::cout << "shape " << tmp3.data().shape(2) << std::endl;
  std::cout << "data " << tmp3.data()(0,0,0) << std::endl;

  //std::cout << "indices " << delta_tau_blocks[0].indices()[0] << std::endl;

  int idx = 0;
  for (auto it = tmp.cbegin(); it != tmp.cend(); ++it) {
    std::cout << "debug2   " << idx << " " << *it << std::endl;
    ++ idx;
  }

  _G0_iw = make_block_gf(block_names, g0_iw_blocks);
  _G_tau = make_block_gf(block_names, g_tau_blocks);
  _G_l = make_block_gf(block_names, g_l_blocks);
  _Delta_tau = make_block_gf(block_names, delta_tau_blocks);
  _G_tau_accum = make_block_gf(block_names, g_tau_accum_blocks);

  //std::cout << "block Gf _rank " << typeid(_G_tau.data()).name() << std::endl;
  //std::cout << "block Gf _rank " << _G_tau.data().data()->data().num_elements() << std::endl;
  //std::cout << "block Gf _rank " << _G_tau.data().data()->mesh().size() << std::endl;
  std::cout << "block Gf _rank " << _G_tau.data().size() << std::endl;
  std::cout << "block Gf _rank " << _G_tau.data()[0].data().rank << std::endl;
  std::cout << "block Gf _rank " << _G_tau.data()[0].data().shape() << std::endl;
  std::cout << "block Gf _rank " << _G_tau.data()[0].data()(0,0,0) << std::endl;
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
      fops.insert(bl.first, a);
    }
  }

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
  boost::multi_array<double, 4> Uijkl_Re(boost::extents[num_flavors][num_flavors][num_flavors][num_flavors]);
  boost::multi_array<double, 4> Uijkl_Im(boost::extents[num_flavors][num_flavors][num_flavors][num_flavors]);
  std::fill(Uijkl_Re.origin(), Uijkl_Re.origin() + Uijkl_Re.num_elements(), 0.0);
  std::fill(Uijkl_Im.origin(), Uijkl_Im.origin() + Uijkl_Im.num_elements(), 0.0);
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
    Uijkl_Re[indices[0]][indices[1]][indices[2]][indices[3]] = (it->coef).real();
    Uijkl_Im[indices[0]][indices[1]][indices[2]][indices[3]] = (it->coef).imag();
  }
  //{
    //std::vector<std::complex<double> > Uijkl_vec(Uijkl.num_elements());
    //std::copy(Uijkl.origin(), Uijkl.origin() + Uijkl.num_elements(), Uijkl_vec.begin());
  //}

  // Do I have imaginary components in my local Hamiltonian?
  auto max_imag = 0.0;
  for (int b : range(gf_struct.size())) {
    max_imag = std::max(max_imag, max_element(abs(real(_G0_iw[b].singularity()(2)))));
  }

  // Add quadratic terms to h_loc
  std::vector<double> h_loc_vec_Re(num_flavors * num_flavors), h_loc_vec_Im(num_flavors * num_flavors);
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
        std::cout << linindex[std::make_pair(b, n1)] << "  "  << linindex[std::make_pair(b,n2)] << " " << e_ij << std::endl;
        const int flavor1 = linindex[std::make_pair(b, n1)];//CORRECT???? FIX ME!
        const int flavor2 = linindex[std::make_pair(b, n2)];
        h_loc_vec_Re[flavor1*num_flavors + flavor2] = e_ij.real();
        h_loc_vec_Im[flavor1*num_flavors + flavor2] = e_ij.imag();

        n2++;
      }
      n1++;
    }
    b++;
    offset += bl.second.size();
  }

  // Determine terms Delta_iw from G0_iw and ensure that the 1/iw behaviour of G0_iw is correct
  {
    b = 0;
    range _;
    triqs::clef::placeholder<0> iw_;
    int offset = 0;
    delta_tau_Re_.resize(boost::extents[n_tau_][num_flavors][num_flavors]);
    delta_tau_Im_.resize(boost::extents[n_tau_][num_flavors][num_flavors]);
    for (auto const &bl : gf_struct) {
      Delta_iw[b](iw_) << G0_iw_inv[b].singularity()(-1) * iw_ + G0_iw_inv[b].singularity()(0);
      Delta_iw[b] = Delta_iw[b] - G0_iw_inv[b];
      _Delta_tau[b]() = inverse_fourier(Delta_iw[b]);
      // Force all diagonal elements to be real
      for (int i : range(bl.second.size())) _Delta_tau[b].data()(_, i, i) = real(_Delta_tau[b].data()(_, i, i));
      // If off-diagonal elements are below threshold, set to real
      if (max_element(abs(imag(_Delta_tau[b].data()))) < params.imag_threshold) {
        _Delta_tau[b].data() = real(_Delta_tau[b].data());
      }
      const auto num_flavors_block = bl.second.size();
      for (auto itau = 0; itau < n_tau_; ++itau) {
        for (auto f1 = 0; f1 < num_flavors_block; ++f1) {
          for (auto f2 = 0; f2 < num_flavors_block; ++f2) {
            delta_tau_Re_[itau][f1+offset][f2+offset] = _Delta_tau[b].data()(itau,f1,f2).real();
            delta_tau_Im_[itau][f1+offset][f2+offset] = _Delta_tau[b].data()(itau,f1,f2).imag();
          }
        }
      }
      b++;
      offset += num_flavors_block;
    }
  }

  // Determine which solver to be used the real-number solver or the complex-number solver
  alps::params par;
  boost::shared_ptr<alps::cthyb::Solver> p_solver;
  alps::cthyb::MatrixSolver<std::complex<double> >::define_parameters(par);
  //if (par.help_requested(std::cout)) { exit(0); } //If help message is requested, print it and exit normally.

  // Set parameters for ALPS CT-HYB solver
  par["timelimit"] = static_cast<unsigned long>(
      (params.max_time > 0 && params.max_time < ULONG_MAX)
      ? params.max_time : ULONG_MAX
  );
  par["verbose"] = params.verbosity == 0 ? 0 : 1;
  par["SEED"] = params.random_seed;

  par["model.sites"] = num_flavors;
  par["model.spins"] = 1;
  par["model.beta"] = beta;
  par["model.command_line_mode"] = true;
  par["model.coulomb_tensor_Re"] = std::vector<double>(Uijkl_Re.origin(), Uijkl_Re.origin() + Uijkl_Re.num_elements());
  par["model.coulomb_tensor_Im"] = std::vector<double>(Uijkl_Im.origin(), Uijkl_Im.origin() + Uijkl_Im.num_elements());
  par["model.hopping_matrix_Re"] = h_loc_vec_Re;
  par["model.hopping_matrix_Im"] = h_loc_vec_Im;
  par["model.n_tau_hyb"] = n_tau_ - 1;
  par["model.delta_Re"] = std::vector<double>(delta_tau_Re_.origin(), delta_tau_Re_.origin()+delta_tau_Re_.num_elements());
  par["model.delta_Im"] = std::vector<double>(delta_tau_Im_.origin(), delta_tau_Im_.origin()+delta_tau_Im_.num_elements());

  par["measurement.G1.n_legendre"] = n_l_;
  par["measurement.G1.n_tau"] = n_tau_ - 1;//Note: the minus 1
  par["measurement.G1.n_matsubara"] = n_iw_;

  if (params.verbosity >= 1) {
    std::cout << "Parameters passed to ALPS CT-HYB solver are the following." << std::endl;
    std::cout << par << std::endl;
  }

  // Call the ALPS CT-HYB solver
  p_solver.reset(new alps::cthyb::MatrixSolver<std::complex<double> >(par));
  p_solver->solve();
  const auto &alps_results = p_solver->get_results();

  if (params.verbosity >= 2) {
    std::cout << "Average sign: " << boost::any_cast<double>(alps_results.at("Sign")) << std::endl;
  }

  // Copy local (real or complex) G_tau back to complex G_tau
  using alps_gtau_t = alps::cthyb::MatrixSolver<std::complex<double> >::G1_tau_t;
  detail::copy_from_alps_to_triqs_gf(
      boost::any_cast<const alps_gtau_t&>(alps_results.at("gf")),
      _G_tau
  );

  //_G_tau = _G_tau_accum;
}

}
