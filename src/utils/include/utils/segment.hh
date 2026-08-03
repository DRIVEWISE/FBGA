#pragma once

#include <utils/types.hh>
#include <string>
#include <cmath>

#include <utils/gg_utils.hh>

namespace fb::utils
{

  template<typename T>
  constexpr T eval_v(const T& S, const T& A, const T& V0) {
      return std::sqrt(static_cast<T>(2) * A * S + std::pow(V0, 2));
  }

  template<typename T>
  constexpr T eval_v2(const T& S, const T& A, const T& V0) {
      return static_cast<T>(2) * A * S + std::pow(V0, 2);
  }
  class Segment
  {
  private:
    real m_s0{0};
    real m_s1{0};
    real m_v0{0};
    real m_a{0};
    real m_L{0};
    real m_k0{0};
    real m_k1{0};
    real m_T{-1.0}; // Total time
    real m_T0{-1.0}; // Total time
    real m_T1{-1.0}; // Total time
    SegmentType m_type{SegmentType::UNKNOWN};
  public:
    // constructors
    Segment() = default;
    //
    Segment(const real L, const real v0, const real k0, const real k1)
    : m_s1(L), m_v0(v0), m_L(L), m_k0(k0), m_k1(k1) {}
    //
    Segment(const real s0, const real L, const real v0, const real k0, const real k1)
    : m_s0(s0), m_s1(s0 + L), m_v0(v0), m_L(L), m_k0(k0), m_k1(k1) {}
    // setters
    void set_s0(const real s0){ this->m_s0 = s0; }
    void set_s1(const real s1){ this->m_s1 = s1; }
    void set_a(const real a){ this->m_a = a; }
    void set_v0(const real v0){ this->m_v0 = v0; }
    void set_L(const real L){ this->m_L = L; }
    void set_k0(const real k0){ this->m_k0 = k0; }
    void set_k1(const real k1){ this->m_k1 = k1; }
    void set_type(const SegmentType &type){ this->m_type = type; }
    // getters
    [[nodiscard]] real s0() const { return this->m_s0; }
    [[nodiscard]] real s1() const { return this->m_s1; }
    [[nodiscard]] real v0() const { return this->m_v0; }
    [[nodiscard]] real a() const { return this->m_a; }
    [[nodiscard]] real L() const { return this->m_L; }
    [[nodiscard]] real k0() const { return this->m_k0; }
    [[nodiscard]] real k1() const { return this->m_k1; }
    [[nodiscard]] SegmentType type() const { return this->m_type; }
    [[nodiscard]] real V(const real s) const { return eval_v(s, this->m_a, this->m_v0); }
    // evaluate v at L (final velocity)
    [[nodiscard]] real VF() const { return eval_v(this->m_L, this->m_a, this->m_v0); }
    // evaluate v (final) given a (forward propagation)
    [[nodiscard]] real VA(const real a) const { return eval_v(this->m_L, a, this->m_v0); }
    // evaluate v (initial) given b and vf (backward propagation)
    [[nodiscard]] real VB(const real b, const real vf) const { return sqrt(-static_cast<real>(2)*b*this->m_L+std::pow(vf,2)); }
    // evaluate t at s
    [[nodiscard]] real t(const real s) const { return static_cast<real>(2)*s / (this->m_v0 + eval_v(s,this->m_a,this->m_v0)); }
    // evaluate Total time
    [[nodiscard]] real T() const { return static_cast<real>(2)*this->m_L / (this->m_v0 + eval_v(this->m_L,this->m_a,this->m_v0)); }
    // evaluate curvature at s
    [[nodiscard]] real K(const real s) const { return this->m_k0 + s*(this->m_k1-this->m_k0)/this->m_L; }
    // evaluate lateral acceleration at s (ay = k*v^2)
    [[nodiscard]] real AY(const real s) const { return this->K(s)*std::pow(this->V(s),2); }
    // evaluate lateral acceleration at s (constant)
    [[nodiscard]] real AX(const real s) const { return this->m_a + 0*s; }
    // evaluate lateral acceleration at final point (forward)
    [[nodiscard]] real AYA(const real a) const { return this->m_k1 * std::pow(this->VA(a),2); }
    // evaluate lateral acceleration at initial point (backward)
    [[nodiscard]] real AYB(const real b, const real vf) const { return this->m_k0 * std::pow(this->VB(b,vf),2); }
    // retrieve initial lateral acceleration
    [[nodiscard]] real AY0() const { return this->m_k0 * std::pow(this->m_v0,2); }
    // retrieve final lateral acceleration
    [[nodiscard]] real AYF() const { return this->m_k1 * std::pow(this->VF(),2); }
    // eval t (generic)
    static real eval_t(const real s, const real a, const real v0) { return static_cast<real>(2)*s / (v0 + eval_v(s,a,v0)); }
    [[nodiscard]] real S(const real t_tot) const
    {
      real const t = t_tot - this->m_T0;
      return 1.0/2.0*this->m_a*t*t + this->m_v0*t + this->m_s0;
    }
    void set_T0(const real T0) { this->m_T0 = T0; }
    void set_T1(const real T1) { this->m_T1 = T1;}
    void set_T (const real T) { this->m_T = T; }
    void set_times(const real T0)
    {
      real const T = this->T();
      this->set_T0(T0);
      this->set_T1(T0 + T);
      this->set_T(T);
    }
    [[nodiscard]] real getT() const { return this->m_T; }
    [[nodiscard]] real getT0() const { return this->m_T0; }
    [[nodiscard]] real getT1() const { return this->m_T1; }

  };

}
