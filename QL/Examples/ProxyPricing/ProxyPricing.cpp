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

// Notes: This file has been adapted from an implementation by Peter Caspers. That implementation is for LSMC proxy pricing for a Bermudan swaption. This file implements CVA estimation of a portfolio of Bermudan swaptions using LSMC proxy pricing for each instrument. 

// The model assumes the following:
// (i) rates are modeled with a one factor short rate model
// (ii) each rate corresponding to each xIBOR curve are uncorrelated. 
// (iii) the default model is Poisson with a fixed hazard rate 
//
// Output: the code compares the CVA accuracy with LSMC and under a closed form solution. The timings of the LSMC are also given. 
// Further details of the model are given in docs/CVA_technical_note.pdf
//
// Author: Matthew Dixon
// Date:   December the 2nd, 2018

#include <ql/quantlib.hpp>
#include <boost/timer.hpp>
#include <ql/termstructures/credit/piecewisedefaultcurve.hpp>
#include <ql/termstructures/credit/defaultprobabilityhelpers.hpp>
#include <ql/termstructures/credit/flathazardrate.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/quotes/simplequote.hpp>

using namespace QuantLib;

#if defined(QL_ENABLE_SESSIONS)
namespace QuantLib {
Integer sessionId() { return 0; }
}
#endif

class Timer {
    boost::timer timer_;
    double elapsed_;

  public:
    void start() { timer_ = boost::timer(); }
    void stop() { elapsed_ = timer_.elapsed(); }
    double elapsed() const { return elapsed_; }
};

// reusable code snipped to perform one pricing step
#define PROXY_PRICING \
        npvRef = swaptionRef2->NPV(); \
        timer.start(); \
        npvProxy = swaption2->NPV(); \
        /*underlyingProxy = 0.0; swaption2->result<Real>("exerciseValue");*/ \
        timer.stop(); \
        npvProxyTiming = timer.elapsed(); \
        std::clog << "\nPricing results on " \
                  << Settings::instance().evaluationDate() \
                  << ", reference rate " << rateLevelRefQuote->value() \
                  << " with maturity " << maturityRefQuote->value() << "\n"; \
        std::clog << "Integral engine npv = " << npvRef << "\n"; \
        std::clog << "Proxy    engine npv = " << npvProxy \
                  << " (timing: " << npvProxyTiming*1000000.0 \
        << "mus)" /*<< ", underlying npv = " << underlyingProxy*/ << "\n";

// here the main part of the code starts

int main(int argc, char *argv[]) {

    try {

        std::clog << "Bermudan swaption proxy pricing example\n";

        // the index for the swaption's underlying.
        // one should again use the fixed yts ...
        std::list<boost::shared_ptr<IborIndex>> indices; //({AUDLibor, CADLibor, Cdor, CHFLibor, DKKLibor, Eonia, Euribor, EURLibor, Fedfunds, GBPLibor, Jibar, JPYLibor, Libor, NZDLibor, SEKLibor,Shibor,Sonia,Tibor,TRLibor,USDLibor,Zibor});
        
        // but for the reference pricings in the integral engine
        // we need a floating version as well.
        std::list<boost::shared_ptr<IborIndex>> indicesRef;


        // Timer seutp

        Timer timer;

        // Original evaluation date and rate level

        Real rateLevelOrig = 0.02;
        Date refDateOrig(12, January, 2015);
        Real recovery_rate = 0.5;

        Calendar calendar = TARGET();
        Date todaysDate(12, June, 2015);
    // must be a business day
        todaysDate = calendar.adjust(todaysDate);


        Settings::instance().evaluationDate() = refDateOrig;
        boost::shared_ptr<Quote> flatRate(new SimpleQuote(0.01));

        Handle<YieldTermStructure> tsCurve(
        boost::make_shared<FlatForward>(
            todaysDate, Handle<Quote>(flatRate), Actual365Fixed()));

    /*
      In Lehmans Brothers "guide to exotic credit derivatives"
      p. 32 there's a simple case, zero flat curve with a flat CDS
      curve with constant market spreads of 150 bp and RR = 50%
      corresponds to a flat 3% hazard rate. The implied 1-year
      survival probability is 97.04% and the 2-years is 94.18%
    */

    // market
    
        Real quoted_spread = 0.0150;
        std::vector<Period> tenors;
        for (Size i=0;i<10;i++){
            tenors.push_back(i * Years);
        }
        

        /*std::vector<boost::shared_ptr<DefaultProbabilityHelper> > instruments;
        for (Size i = 0; i < 10; i++) {
            instruments.push_back(boost::shared_ptr<DefaultProbabilityHelper>(
                new SpreadCdsHelper(Handle<Quote>(boost::shared_ptr<Quote>(
                                    new SimpleQuote(quoted_spread))),
                                tenors[i], 0, calendar, Quarterly, Following,
                                DateGeneration::TwentiethIMM, Actual365Fixed(),
                                recovery_rate, tsCurve)));

        }*/

    // Bootstrap hazard rates
        /*boost::shared_ptr<PiecewiseDefaultCurve<HazardRate, BackwardFlat> >
            hazardRateStructure(new PiecewiseDefaultCurve<HazardRate, BackwardFlat>(
        todaysDate, instruments, Actual365Fixed()));*/


        // the yield term structure for the original pricing
        // this must _not_ be floating, see the warning in
        // the proxy engine's header.

        Handle<YieldTermStructure> ytsOrig(boost::make_shared<FlatForward>(
            refDateOrig, rateLevelOrig, Actual365Fixed()));


        // the yield term structure for reference pricings
        // (integral engine) on future dates, this _floating_

        boost::shared_ptr<SimpleQuote> rateLevelRefQuote =
            boost::make_shared<SimpleQuote>(0.02);
        Handle<Quote> rateLevelRef(rateLevelRefQuote);
        Handle<YieldTermStructure> ytsRef(boost::make_shared<FlatForward>(
            0, TARGET(), rateLevelRef, Actual365Fixed()));
        // Construct the target xIBOR portfolio. To do: This should be placed in a for loop 

        indices.push_back(boost::make_shared<USDLibor>(6 * Months, ytsOrig));
        indicesRef.push_back(boost::make_shared<USDLibor>(6 * Months, ytsOrig));
        
        indices.push_back(boost::make_shared<AUDLibor>(6 * Months, ytsOrig));
        indicesRef.push_back(boost::make_shared<AUDLibor>(6 * Months, ytsOrig));

        indices.push_back(boost::make_shared<CADLibor>(6 * Months, ytsOrig));
        indicesRef.push_back(boost::make_shared<CADLibor>(6 * Months, ytsOrig));
        
        indices.push_back(boost::make_shared<SEKLibor>(6 * Months, ytsOrig));
        indicesRef.push_back(boost::make_shared<SEKLibor>(6 * Months, ytsOrig));
        
        indices.push_back(boost::make_shared<CHFLibor>(6 * Months, ytsOrig));
        indicesRef.push_back(boost::make_shared<CHFLibor>(6 * Months, ytsOrig));
        
        indices.push_back(boost::make_shared<DKKLibor>(6 * Months, ytsOrig));
        indicesRef.push_back(boost::make_shared<DKKLibor>(6 * Months, ytsOrig));
        
        indices.push_back(boost::make_shared<EURLibor>(6 * Months, ytsOrig));
        indicesRef.push_back(boost::make_shared<EURLibor>(6 * Months, ytsOrig));

        indices.push_back(boost::make_shared<JPYLibor>(6 * Months, ytsOrig));
        indicesRef.push_back(boost::make_shared<JPYLibor>(6 * Months, ytsOrig));

        indices.push_back(boost::make_shared<NZDLibor>(6 * Months, ytsOrig));
        indicesRef.push_back(boost::make_shared<NZDLibor>(6 * Months, ytsOrig));


        indices.push_back(boost::make_shared<TRLibor>(6 * Months, ytsOrig));
        indicesRef.push_back(boost::make_shared<TRLibor>(6 * Months, ytsOrig));

        indices.push_back(boost::make_shared<Shibor>(6 * Months, ytsOrig));
        indicesRef.push_back(boost::make_shared<Shibor>(6 * Months, ytsOrig));

        indices.push_back(boost::make_shared<GBPLibor>(6 * Months, ytsOrig));
        indicesRef.push_back(boost::make_shared<GBPLibor>(6 * Months, ytsOrig));
        //indices.push_back(boost::make_shared<Eonia>(6 * Months, ytsOrig));
        indices.push_back(boost::make_shared<Cdor>(6 * Months, ytsOrig));
        indicesRef.push_back(boost::make_shared<Cdor>(6 * Months, ytsOrig));


        //indices.push_back(boost::make_shared<FedFunds>(6 * Months, ytsOrig));
        indices.push_back(boost::make_shared<Jibar>(6 * Months, ytsOrig));
        indicesRef.push_back(boost::make_shared<Jibar>(6 * Months, ytsOrig));

        indices.push_back(boost::make_shared<Tibor>(6 * Months, ytsOrig));
        indicesRef.push_back(boost::make_shared<Tibor>(6 * Months, ytsOrig));

        // Configure the Bermudan swaptions for the reference XIBOR curves
        // the maturity of the bermudan swaption in years    
        Size length = 10;
                // instrument setup
        Real strike = 0.02; // near atm option
        Date effectiveDate = TARGET().advance(refDateOrig, 2 * Days);
        Date startDate = TARGET().advance(effectiveDate, 1 * Years);
        Date maturityDate = TARGET().advance(startDate, length * Years);

        Schedule fixedSchedule(startDate, maturityDate, 1 * Years, TARGET(),
                               ModifiedFollowing, ModifiedFollowing,
                               DateGeneration::Forward, false);
        Schedule floatingSchedule(startDate, maturityDate, 6 * Months, TARGET(),
                                  ModifiedFollowing, ModifiedFollowing,
                                  DateGeneration::Forward, false);


        std::vector<Date> exerciseDates;
        for (Size i = 0; i < length; ++i) {
            exerciseDates.push_back(
                TARGET().advance(fixedSchedule[i], -2 * Days));
        }

        boost::shared_ptr<Exercise> exercise =
            boost::make_shared<BermudanExercise>(exerciseDates, false);

        std::list<boost::shared_ptr<IborIndex>>::iterator it;
        std::list<boost::shared_ptr<VanillaSwap>> underlyings;
        std::list<boost::shared_ptr<Swaption>> swaptions; 
        std::list<boost::shared_ptr<NonstandardSwaption>> swaptions2; 
        
        
        for(it = indices.begin(); it != indices.end(); ++it){
            
            underlyings.push_back(boost::make_shared<VanillaSwap>(VanillaSwap(VanillaSwap::Payer, 1.0, fixedSchedule, strike, Thirty360(),
                floatingSchedule, *it, 0.0, Actual360())));

            swaptions.push_back(boost::make_shared<Swaption>(underlyings.back(), exercise));
            swaptions2.push_back(boost::make_shared<NonstandardSwaption>(*swaptions.back()));
        }

        std::list<boost::shared_ptr<VanillaSwap>> underlyingRefs;
        std::list<boost::shared_ptr<Swaption>> swaptionRefs; 
        std::list<boost::shared_ptr<NonstandardSwaption>> swaption2Refs; 
        
        for(it = indicesRef.begin(); it != indicesRef.end(); ++it){
             underlyingRefs.push_back(boost::make_shared<VanillaSwap>(VanillaSwap(
                VanillaSwap::Payer, 1.0, fixedSchedule, strike, Thirty360(),
                floatingSchedule, *it, 0.0, Actual360())));
             swaptionRefs.push_back(boost::make_shared<Swaption>(underlyingRefs.back(), exercise));
             swaption2Refs.push_back(boost::make_shared<NonstandardSwaption>(*swaptionRefs.back()));
        }

    
        // our instrument is a swaption, but the engine is for non standard
        // swaptions
        // so we just convert it

        //boost::shared_ptr<NonstandardSwaption> swaption2 =
        //    boost::make_shared<NonstandardSwaption>(*swaption);
        //boost::shared_ptr<NonstandardSwaption> swaptionRef2 =
        //    boost::make_shared<NonstandardSwaption>(*swaptionRef);

        // just take any model volatility and reversion, we do not calibrate
        // them here. Also they are flat, so no steps needed really.

        std::vector<Date> stepDates;
        std::vector<Real> sigmas(1, 0.0070);
        Real reversion = 0.0030;
        Period oneMonth(1,Months);

        // the gsr model in T-forward measure, T=50 chosen arbitrary here
        // the first model uses the fixed yts, used for the mc pricing
        // generating the proxy

        boost::shared_ptr<Gsr> gsrFixed = boost::make_shared<Gsr>(
            ytsOrig, stepDates, sigmas, reversion, 50.0);

        // the second model is used for the reference pricing, therefore
        // using the floating yts

        boost::shared_ptr<Gsr> gsrFloating =
            boost::make_shared<Gsr>(ytsRef, stepDates, sigmas, reversion, 50.0);

        // the integral engine for reference pricings

        boost::shared_ptr<PricingEngine> integralEngine =
            boost::make_shared<Gaussian1dNonstandardSwaptionEngine>(
                gsrFloating, 64, 7.0, true, false, Handle<Quote>(), ytsRef);

        // compute a reference price for the inital pricing

        timer.start();
        std::list<boost::shared_ptr<NonstandardSwaption>>::iterator it_;
        Real npvOrigIntegral; //npvRef;

        boost::shared_ptr<SimpleQuote> maturityRefQuote =
            boost::make_shared<SimpleQuote>();
        Handle<Quote> maturityRef(maturityRefQuote);
        Real cvaOrig = 0.0;
        for(it_ = swaption2Refs.begin(); it_ != swaption2Refs.end(); ++it_){
            (*it_)->setPricingEngine(integralEngine);

            Date currentDate = todaysDate;

            for (int i=0; i<1;i++)
            {
                currentDate += oneMonth;
                Settings::instance().evaluationDate() = currentDate;
                rateLevelRefQuote->setValue(0.02); // no change
                maturityRefQuote->setValue(10.5-float(i)/360.0);  // maturity of the underlying
                npvOrigIntegral = (*it_)->NPV();
                cvaOrig += (1.0-exp(-0.01*i)) *(1-recovery_rate)*npvOrigIntegral;//Analytic CVA estimation as in Eq 1.12 in the technical note
                //cvaOrig += hazardRateStructure->defaultProbability(currentDate)*(1-recovery_rate)*npvOrigIntegral;
                std::clog << "Integral engine npv = " << npvOrigIntegral << "\n"; \
            }
        }

        timer.stop();
        Real npvOrigIntegralTiming = timer.elapsed();

        // the mc engine, note that we use the fixed model here

 // reference maturity for the scenario rate

        

        timer.start();
        Real npvOrigMc;
        Real errorOrigMc;
        Real npvProxy; //npvProxyTiming;
        Real cvaProxy = 0.0;

        
        for(it_ = swaptions2.begin(); it_ != swaptions2.end(); ++it_){

            boost::shared_ptr<PricingEngine> mcEngine =
            MakeMcGaussian1dNonstandardSwaptionEngine<>(gsrFixed)
                .withSteps(1) // the gsr model allows for large steps
                .withSamples(10000)
                .withSeed(42)
                .withCalibrationSamples(10000)
                .withProxy(true);
            (*it_)->setPricingEngine(mcEngine);
            npvOrigMc = (*it_)->NPV();
            errorOrigMc = (*it_)->errorEstimate();

            boost::shared_ptr<PricingEngine> proxyEngine =
            boost::make_shared<ProxyNonstandardSwaptionEngine>(
                (*it_)->proxy(), rateLevelRef, maturityRef, 64, 7.0, false);

            Date currentDate = todaysDate;
            for (int i=0; i<1;i++)
            {
                currentDate += oneMonth;
                Settings::instance().evaluationDate() = currentDate;
                rateLevelRefQuote->setValue(0.02); // no change
                maturityRefQuote->setValue(10.5-float(i)/360.0);  // maturity of the underlying
                npvProxy = (*it_)->NPV();
                cvaProxy += (1.0-exp(-0.01*i))*(1-recovery_rate)*npvProxy; // CVA proxy estimation as in Eq 1.12 in the technical note
                //cvaProxy += hazardRateStructure->defaultProbability(currentDate)*(1-recovery_rate)*npvProxy;
                //errorProxyMc = (*it)->errorEstimate();
                std::clog << "Proxy engine npv = " << npvProxy << "\n"; \
            }
            

        }
        timer.stop();
        
        // compute the mc price
        Real npvOrigMcTiming = timer.elapsed();

        // output the results

        //std::clog << "Pricing results on the original reference date ("
        //          << refDateOrig << "):\n";
        std::clog << "Integral engine cva = " << cvaOrig
                  << " (timing: " << npvOrigIntegralTiming*1000000.0 << "mus)\n";
        std::clog << "MC       engine cva = " << cvaProxy
                  <<" (timing: " << npvOrigMcTiming*1000000.0 << "mus)\n";


        // Settings::instance().evaluationDate() = Date(10, January, 2020);
        // rateLevelRefQuote->setValue(0.005);
        // maturityRefQuote->setValue(6.0);
        // npvProxy = swaption2->NPV();
        // underlyingProxy = swaption2->result<Real>("exerciseValue");
        // std::clog << "\nExercise check (" << Settings::instance().evaluationDate() << "):\n";
        // std::clog << "otm option: exercise value=" << underlyingProxy << " npv=" << npvProxy << std::endl;

        // rateLevelRefQuote->setValue(0.04);
        // npvProxy = swaption2->NPV();
        // underlyingProxy = swaption2->result<Real>("exerciseValue");
        // std::clog << "itm option: exercise value=" << underlyingProxy << " npv=" << npvProxy << std::endl;

        return 0;

    } catch (QuantLib::Error e) {
        std::clog << "terminated with a ql exception: " << e.what()
                  << std::endl;
        return 1;
    } catch (std::exception e) {
        std::clog << "terminated with a general exception: " << e.what()
                  << std::endl;
        return 1;
    }
}
