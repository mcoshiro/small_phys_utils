/**
 * Class to wrap correctionlib corrections since Cling (ROOT 24) has conflicts with C++17
 */
#ifndef H_CORRECTION_LIB
#define H_CORRECTION_LIB

#include <memory>
#include <string>
#include <vector>

class CorrectionWrapper {
  public:
    CorrectionWrapper(std::string filename, std::string correction_name);
    ~CorrectionWrapper();
    double evaluate(std::string value_name, std::vector<double> values) const;
    double evaluate(std::vector<double> values) const;

  private:
    //kludge to avoid including C++17 dependencies in this file
    void* cs_voidptr;
    void* correction_map_voidptr;
};

#endif 
