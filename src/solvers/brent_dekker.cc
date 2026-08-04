#include "solvers/brent_dekker.hh"

namespace fb::solvers
{

BrentDekker::BrentDekker(const real tol) : tol(tol), verbose("iter"){ }

BrentDekker::BrentDekker(const real tol, const int max_iter) : tol(tol), max_iter(max_iter), verbose("iter"){}

BrentDekker::BrentDekker(const real tol, const int max_iter, const std::string &verbose) : tol(tol), max_iter(max_iter), verbose(verbose){}

void BrentDekker::setup(const real tol, const int max_iter, const std::string &verbose)
{
  this->set_tolerance(tol);
  this->set_max_iter(max_iter);
  this->set_verbose(verbose);
}

} // namespace GG
