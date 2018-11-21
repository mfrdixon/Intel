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

#include <ql/experimental/volatility/sabrzc.hpp>

namespace QuantLib {

namespace detail {
namesapce ZcSabr {

    Real G(const Real t, const Real s) {
        const Real g = s * std::coth(s) - 1.0;
        const Real R = 1.0 + 3.0 * t * g / (8.0 * s * s) -
                       5.0 * t * t * (-8.0 * s * s + 3.0 * g * g + 24.0 * g) /
                           (128.0 * s * s * s * s) +
                       35.0 * t * t * t * (-40.0 * s * s + 3.0 * g * g * g +
                                           24.0 * g * g + 120.0 * g) /
                           (1024.0 * s * s * s * s * s * s);
        const Real dR =
            std::exp(t / 8.0) -
            (3072.0 + 384.0 * t + 24.0 * t * t + t * t * t) / 3072.0;
        return std::sqrt(std::sinh(s) / s) *
               std::exp(-0.5 * s * s / t - t / 8.0) * (R + dR);
    };

} // namespace ZcSabr
} // namesapce detail

    ZcSabr::ZcSabr(const Real expiryTime, const Real forward, const Real alpha,
           const Real beta, const Real nu




} // namesapce QuantLib
