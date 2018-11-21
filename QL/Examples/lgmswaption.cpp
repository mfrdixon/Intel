#include <ql/quantlib.hpp>

using namespace QuantLib;

int main() {

    try {

        // ------------------------------------------
        // read parameters
        // ------------------------------------------

        char *p_intpoints = getenv("intpoints");
        char *p_stddevs = getenv("stddevs");
        Size intpoints = p_intpoints != NULL ? atoi(p_intpoints) : 24;
        double stddevs = p_stddevs != NULL ? atof(p_stddevs) : 4.0;

        Real rateLevel = 0.02;
        Date evalDate(12, January, 2015);

        Settings::instance().evaluationDate() = evalDate;

        std::clog << "reference date " << evalDate << std::endl;

        Handle<YieldTermStructure> yts(boost::make_shared<FlatForward>(
            evalDate, rateLevel, Actual365Fixed()));

        boost::shared_ptr<IborIndex> euribor6m =
            boost::make_shared<Euribor>(6 * Months, yts);

        // -------------------------------------------
        // set up a bunch of swaptions
        // -------------------------------------------

        std::vector<boost::shared_ptr<Swaption> > swaptions;

        std::vector<Date> exerciseDates;

        Date effectiveDate = TARGET().advance(evalDate, 2 * Days);
        Date startDate = TARGET().advance(effectiveDate, 1 * Years);
        Date maturityDate =
            TARGET().advance(startDate, 10 * Years); // maturity here !!!

        for (Size i = 0; i < 10; ++i) {

            Real strike = 0.02 + static_cast<Real>(i) * 0.0010;

            Schedule fixedSchedule(startDate, maturityDate, 1 * Years, TARGET(),
                                   ModifiedFollowing, ModifiedFollowing,
                                   DateGeneration::Forward, false);
            Schedule floatingSchedule(startDate, maturityDate, 6 * Months,
                                      TARGET(), ModifiedFollowing,
                                      ModifiedFollowing,
                                      DateGeneration::Forward, false);

            // std::clog << "exercise dates:" << std::endl;
            exerciseDates.clear();
            for (Size ii = 0; ii < 9; ++ii) {
                exerciseDates.push_back(
                    TARGET().advance(fixedSchedule[ii], -2 * Days));
                // std::clog << exerciseDates.back() << "\n";
            }

            boost::shared_ptr<Exercise> exercise =
                boost::make_shared<BermudanExercise>(exerciseDates, false);

            boost::shared_ptr<VanillaSwap> underlying =
                boost::make_shared<VanillaSwap>(VanillaSwap(
                    VanillaSwap::Payer, 1.0, fixedSchedule, strike, Thirty360(),
                    floatingSchedule, euribor6m, 0.0, Actual360()));

            boost::shared_ptr<Swaption> swaption =
                boost::make_shared<Swaption>(underlying, exercise);

            swaptions.push_back(swaption);
        }

        // ------------------------------------
        // models and calibration swaptions
        // ------------------------------------

        std::vector<Date> stepDates(exerciseDates.begin(),
                                    exerciseDates.end() - 1);

        std::vector<Real> sigmas(stepDates.size() + 1, 0.0050);
        Real reversion = 0.01;

        boost::shared_ptr<Gsr> gsr =
            boost::make_shared<Gsr>(yts, stepDates, sigmas, reversion, 50.0);

        boost::shared_ptr<Lgm1> lgm =
            boost::make_shared<Lgm1>(yts, stepDates, sigmas, reversion);

        boost::shared_ptr<PricingEngine> swaptionEngineGsr =
            boost::make_shared<Gaussian1dSwaptionEngine>(gsr, intpoints,
                                                         stddevs, false, false);

        boost::shared_ptr<PricingEngine> swaptionEngineLgm =
            boost::make_shared<Gaussian1dSwaptionEngine>(lgm, intpoints,
                                                         stddevs, false, false);

        boost::shared_ptr<PricingEngine> swaptionEngineLgmAD =
            boost::make_shared<LgmSwaptionEngineAD>(lgm, intpoints, stddevs);


        std::vector<boost::shared_ptr<CalibrationHelper> > basket;
        std::vector<boost::shared_ptr<SimpleQuote> > marketQuotes;
        std::vector<Real> marketVegas;
        for (Size i = 0; i < exerciseDates.size(); ++i) {
            marketQuotes.push_back(boost::make_shared<SimpleQuote>(0.20));
            boost::shared_ptr<SwaptionHelper> tmp =
                boost::make_shared<SwaptionHelper>(
                    exerciseDates[i], maturityDate,
                    Handle<Quote>(marketQuotes[i]), euribor6m, 1 * Years,
                    Thirty360(), Actual360(), yts);
            tmp->setPricingEngine(swaptionEngineLgmAD);
            basket.push_back(tmp);
            // compute black vega of helpers
            boost::shared_ptr<PricingEngine> blackEngine =
                boost::make_shared<BlackSwaptionEngine>(
                    yts, Handle<Quote>(marketQuotes[i]), Actual365Fixed());
            boost::shared_ptr<Swaption> tmpSwaption = tmp->swaption();
            tmpSwaption->setPricingEngine(blackEngine);
            marketVegas.push_back(tmpSwaption->result<Real>("vega") / 100.0);
        }
        LevenbergMarquardt method;
        EndCriteria ec(1000, 500, 1E-8, 1E-8, 1E-8);
        // lgm->calibrateAlphasIterative(basket, method, ec);
        std::cout << "Calibration results:\n";
        Array lgm_alpha = lgm->alpha();
        for (Size i = 0; i < basket.size(); ++i) {
            std::cout << "#" << i << " model_vol " << lgm_alpha[i]
                      << " helper market " << basket[i]->marketValue()
                      << " model " << basket[i]->modelValue() << " vega "
                      << marketVegas[i] << std::endl;
        }

        // -------------------------------------------------------
        // pricing (usual pricing engine), vega (bump and revalue)
        // -------------------------------------------------------

        // swaptions[0]->setPricingEngine(swaptionEngineLgmAD);
        // Real npv0 = swaptions[0]->NPV();
        // std::cout << "i=0, npv (Lgm) = " << npv0 << std::endl;
        // for (Size i = 0; i < basket.size(); ++i) {
        //     marketQuotes[i]->setValue(0.201);
        //     lgm->calibrateAlphasIterative(basket, method, ec);
        //     Real npvi = swaptions[0]->NPV();
        //     std::cout << "   vega (" << i << ") = " << 10.0*(npvi - npv0) << std::endl;
        //     marketQuotes[i]->setValue(0.20);
        // }
        // lgm->calibrateAlphasIterative(basket, method, ec);

        // ----------------------------------
        // test fortran (ad) engine
        // ----------------------------------

        for (Size i = 0; i < 1; ++i) {
            swaptions[i]->setPricingEngine(swaptionEngineLgmAD);
            Real npvLgmAD = swaptions[i]->NPV();
            std::clog << "i=" << i << " npv (LgmAD) = " << npvLgmAD
                      << std::endl;
            std::vector<Real> times =
                swaptions[i]->result<std::vector<Real> >("sensitivityTimes");
            std::vector<Real> sensH =
                swaptions[i]->result<std::vector<Real> >("sensitivityH");
            std::vector<Real> sensZ =
                swaptions[i]->result<std::vector<Real> >("sensitivityZeta");
            std::vector<Real> sensD =
                swaptions[i]->result<std::vector<Real> >("sensitivityDiscount");
            for (Size j = 0; j < times.size(); ++j) {
                std::clog << i << ";" << j << ";" << times[j] << ";" << sensH[j]
                          << ";" << sensZ[j] << ";" << sensD[j] << std::endl;
            }

            // --------------------------------------
            // compute vega matrix
            // --------------------------------------

            // dNPV / dzeta
            Matrix A(1, basket.size(), 0.0);
            std::vector<Size> calIndices;
            Size ii = 0;
            for (Size i = 0; i < times.size() - 1; ++i) {
                if (!close_enough(std::fabs(sensZ[i + 1]), 0.0)) {
                    A[0][ii] = sensZ[i + 1];
                    std::cout << "found zeta sensi " << ii << " is "
                              << sensZ[i + 1] << std::endl;
                    QL_REQUIRE(ii < basket.size(),
                               "too many sensitivities to zeta !");
                    calIndices.push_back(i + 1);
                    ++ii;
                }
            }
            QL_REQUIRE(ii == basket.size(), "too few sensitivities to zeta !");

            // dImplVol / dzeta
            Matrix B(basket.size(), basket.size(), 0.0);
            for (Size i = 0; i < basket.size(); ++i) {
                boost::shared_ptr<Swaption> tmp =
                    boost::dynamic_pointer_cast<SwaptionHelper>(basket[i])
                        ->swaption();
                tmp->setPricingEngine(swaptionEngineLgmAD);
                std::vector<Real> tmpTimes =
                    tmp->result<std::vector<Real> >("sensitivityTimes");
                std::vector<Real> tmpSensH =
                    tmp->result<std::vector<Real> >("sensitivityH");
                std::vector<Real> tmpSensZ =
                    tmp->result<std::vector<Real> >("sensitivityZeta");
                std::vector<Real> tmpSensD =
                    tmp->result<std::vector<Real> >("sensitivityDiscount");
                Size jj = 0;
                while (jj < tmpTimes.size() - 1 &&
                       !close_enough(tmpTimes[jj + 1], times[calIndices[i]]))
                    ++jj;
                QL_REQUIRE(jj < tmpTimes.size() - 1,
                           "did not find calibration instruments #"
                               << i << " time (" << times[calIndices[i]] << ")")
                B[i][i] = tmpSensZ[jj + 1] / marketVegas[i];
            }
            std::cout << "B=" << B << std::endl;
            // final derivative
            Matrix C = A * inverse(B);
            for (Size i = 0; i < basket.size(); ++i) {
                std::cout << "AD vega (" << i << ") = " << C[0][i] << std::endl;
            }
        } // loop over bermudan swaptions
        return 0;
    } catch (QuantLib::Error e) {
        std::clog << e.what() << std::endl;
    } catch (std::exception e) {
        std::clog << e.what() << std::endl;
    }

} // main
