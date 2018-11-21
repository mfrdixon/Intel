/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Peter Caspers

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file sabrzc.hpp
    \brief zero correlation sabr

    Reference: Antonov, Konikov, Spector, SABR spreads its wings,
               Risk, August 2013
*/

#ifndef quantlib_zcsabr_hpp
#define quantlib_zcsabr_hpp

namespace QuantLib {

namespace detail {

namespace ZcSabr {

Real G(const Real t, const Real s);

Real integrand(const Real expiryTime, const Real forward, const Real alpha,
               const Real beta, const Real nu, const Real shift) {
    const Real eta = 1.0 / (2.0 * (beta - 1.0));
    const Real phi = 2.0 * std::atan(
};

} // namesapce ZcSabr
} // namespace detail

class ZcSabr {
  public:
    ZcSabr(const Real expiryTime, const Real forward, const Real alpha,
           const Real beta, const Real nu);

  private:
    std::pair<Real,Real> phi_psi(const Real s) const;

    const Real expiryTime_, forward_, alpha_, beta_, nu_;
    Real q_, q0_, smin_, smax_;
};

} // namesapce QuantLib

#endif
