// DO NOT EDIT
// Generated automatically using libclang using the command :
// c++2py.py ../c++/solver_core.hpp -p -mpytriqs.applications.impurity_solvers.alps_cthyb -o alps_cthyb --moduledoc "The ALPS cthyb solver" -I /opt/local/include/openmpi-clang38 -I /opt/ALPSCore/include


// --- C++ Python converter for solve_parameters_t
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/string.hpp>
#include <algorithm>

namespace triqs { namespace py_tools {

template <> struct py_converter<solve_parameters_t> {
 static PyObject *c2py(solve_parameters_t const & x) {
  PyObject * d = PyDict_New();
  PyDict_SetItemString( d, "assume_real"   , convert_to_python(x.assume_real));
  PyDict_SetItemString( d, "h_int"         , convert_to_python(x.h_int));
  PyDict_SetItemString( d, "random_seed"   , convert_to_python(x.random_seed));
  PyDict_SetItemString( d, "max_time"      , convert_to_python(x.max_time));
  PyDict_SetItemString( d, "verbosity"     , convert_to_python(x.verbosity));
  PyDict_SetItemString( d, "imag_threshold", convert_to_python(x.imag_threshold));
  return d;
 }

 template <typename T, typename U> static void _get_optional(PyObject *dic, const char *name, T &r, U const &init_default) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = init_default;
 }

 template <typename T> static void _get_optional(PyObject *dic, const char *name, T &r) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = T{};
 }

 static solve_parameters_t py2c(PyObject *dic) {
  solve_parameters_t res;
  res.assume_real = convert_from_python<bool>(PyDict_GetItemString(dic, "assume_real"));
  res.h_int = convert_from_python<many_body_op_t>(PyDict_GetItemString(dic, "h_int"));
  _get_optional(dic, "random_seed"   , res.random_seed      ,34788+928374*triqs::mpi::communicator().rank());
  _get_optional(dic, "max_time"      , res.max_time         ,-1);
  _get_optional(dic, "verbosity"     , res.verbosity        ,0);
  _get_optional(dic, "imag_threshold", res.imag_threshold   ,1.e-15);
  return res;
 }

 template <typename T>
 static void _check(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!convertible_from_python<T>(PyDict_GetItemString(dic, name), false))
   fs << "\n" << ++err << " The parameter " << name << " does not have the right type : expecting " << tname
      << " in C++, but got '" << PyDict_GetItemString(dic, name)->ob_type->tp_name << "' in Python.";
 }

 template <typename T>
 static void _check_mandatory(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!PyDict_Contains(dic, pyref::string(name)))
   fs << "\n" << ++err << " Mandatory parameter " << name << " is missing.";
  else _check<T>(dic,fs,err,name,tname);
 }

 template <typename T>
 static void _check_optional(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (PyDict_Contains(dic, pyref::string(name))) _check<T>(dic, fs, err, name, tname);
 }

 static bool is_convertible(PyObject *dic, bool raise_exception) {
  if (dic == nullptr or !PyDict_Check(dic)) {
   if (raise_exception) { PyErr_SetString(PyExc_TypeError, "The function must be called with named arguments");}
   return false;
  }
  std::stringstream fs, fs2; int err=0;

#ifndef TRIQS_ALLOW_UNUSED_PARAMETERS
  std::vector<std::string> ks, all_keys = {"assume_real","h_int","random_seed","max_time","verbosity","imag_threshold"};
  pyref keys = PyDict_Keys(dic);
  if (!convertible_from_python<std::vector<std::string>>(keys, true)) {
   fs << "\nThe dict keys are not strings";
   goto _error;
  }
  ks = convert_from_python<std::vector<std::string>>(keys);
  for (auto & k : ks)
   if (std::find(all_keys.begin(), all_keys.end(), k) == all_keys.end())
    fs << "\n"<< ++err << " The parameter '" << k << "' is not recognized.";
#endif

  _check_mandatory<bool          >(dic, fs, err, "assume_real"   , "bool");
  _check_mandatory<many_body_op_t>(dic, fs, err, "h_int"         , "many_body_op_t");
  _check_optional <int           >(dic, fs, err, "random_seed"   , "int");
  _check_optional <int           >(dic, fs, err, "max_time"      , "int");
  _check_optional <int           >(dic, fs, err, "verbosity"     , "int");
  _check_optional <double        >(dic, fs, err, "imag_threshold", "double");
  if (err) goto _error;
  return true;

 _error:
   fs2 << "\n---- There " << (err > 1 ? "are " : "is ") << err<< " error"<<(err >1 ?"s" : "")<< " in Python -> C++ transcription for the class solve_parameters_t\n" <<fs.str();
   if (raise_exception) PyErr_SetString(PyExc_TypeError, fs2.str().c_str());
  return false;
 }
};

}}