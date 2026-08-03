#pragma once

#include <utils/types.hh>

#include <vector>
#include <string>

namespace fb::utils
{

  // Enum to understand if a segment is of type forward (1) or backward (2) or transitions (3) or unknown (0)
  enum class SegmentType
  {
    FORWARD = 1,
    BACKWARD = 2,
    TRANSITION = 3,
    FORWARD_NAN = 4,
    BACKWARD_NAN = 5,
    TRANSITION_NAN = 6,
    UNKNOWN = 0
  };

  real clip(real x, real a, real b);

  real piramid(real x, real y);

  real square_conversion(real a, real amin, real amax);

  real signed_distance(real ax, real axmin, real axmax, real ay, real aymin, real aymax);

class LinearInterpolator {
public:
  LinearInterpolator() = default;
  LinearInterpolator(const std::vector<real>& X, const std::vector<real>& Y);
  void setup( const std::vector<real>& X, const std::vector<real>& Y);

  [[nodiscard]] real eval(real xi) const ;

private:
  std::vector<real> m_x, m_y;
  integer m_nx{0};

};

class LinearInterpolatorSet {
public:
  LinearInterpolatorSet() = default;
  LinearInterpolatorSet(const std::vector<real>& X, const std::vector<std::vector<real>>& Y, const std::vector<std::string> &headers);
  void setup( const std::vector<real>& X, const std::vector<std::vector<real>>& Y, const std::vector<std::string> &headers);

  [[nodiscard]] real eval(const std::string & name, real xi) const ;

private:
  std::vector<real> m_x;
  std::vector<std::vector<real>> m_y;
  std::vector<std::string> m_headers;
  integer m_nx{0};
  integer m_ny{0};

};

size_type find_interval(const std::vector<real>& vec, real value) ;

size_type find_interval_with_guess(const std::vector<real>& vec, real value, size_type guess) ;

size_type find_interval_binary_search(const std::vector<real>& vec, real value) ;

size_type find_interval_by_equispaced(const std::vector<real>& vec, real value) ;

class BilinearInterpolator {
public:
  BilinearInterpolator() = default;
  BilinearInterpolator(const std::vector<real>& X,
                        const std::vector<real>& Y,
                        const std::vector<real>& Z);
  void setup( const std::vector<real>& X,
              const std::vector<real>& Y,
              const std::vector<real>& Z,
              bool transpose = false);

  [[nodiscard]] real eval(real xi, real yi) const ;

private:
  std::vector<real> m_x, m_y, m_z;
  integer m_nx{0}, m_ny{0};

};

void compute_finite_difference(const std::vector<double>& X, const std::vector<double>& Y, std::vector<real> & dY_dx);
std::vector<real> compute_finite_difference(const std::vector<double>& X, const std::vector<double>& Y);

void linspace(std::vector<real> & vec, real x1, real x2, integer N);

} // namespace fb::utils
