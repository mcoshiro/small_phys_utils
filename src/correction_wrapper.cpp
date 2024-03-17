/**
 * Class to wrap correctionlib corrections since Cling (ROOT 24) has conflicts with C++17
 */
#include "correction_wrapper.hpp"
#include "correction.hpp"
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using std::string;
using std::vector;
using std::unique_ptr;
using correction_mod::CorrectionSet;
using correction_mod::Variable;
using Ref = correction_mod::Correction::Ref;

CorrectionWrapper::CorrectionWrapper(string filename, string correction_name) {
  unique_ptr<CorrectionSet>* cs = new unique_ptr<CorrectionSet>(CorrectionSet::from_file(filename));
  Ref* correction_map = new Ref((*cs)->at(correction_name.c_str()));
  cs_voidptr = static_cast<void*>(cs);
  correction_map_voidptr = static_cast<void*>(correction_map);
}

CorrectionWrapper::~CorrectionWrapper() {
  unique_ptr<CorrectionSet>* cs = static_cast<unique_ptr<CorrectionSet>*>(cs_voidptr);
  Ref* correction_map = static_cast<Ref*>(correction_map_voidptr);
  delete cs;
  delete correction_map;
}

double CorrectionWrapper::evaluate(string value_name, vector<double> values) const {
  Ref* correction_map = static_cast<Ref*>(correction_map_voidptr);
  std::vector<Variable::Type> evaluate_argument;
  evaluate_argument.push_back(value_name);
  for (double value : values) {
    evaluate_argument.push_back(value);
  }
  return (*correction_map)->evaluate(evaluate_argument);
}

double CorrectionWrapper::evaluate(vector<double> values) const {
  Ref* correction_map = static_cast<Ref*>(correction_map_voidptr);
  std::vector<Variable::Type> evaluate_argument;
  for (double value : values) {
    evaluate_argument.push_back(value);
  }
  return (*correction_map)->evaluate(evaluate_argument);
}
