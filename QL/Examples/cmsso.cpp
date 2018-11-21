#include <ql/quantlib.hpp>

using namespace QuantLib;

int main() {

    double sigma = atof(getenv("SIGMA"));
    double shift = atof(getenv("SHIFT"));
    double rho = atof(getenv("RHO"));
    int normal = atoi(getenv("NORMAL"));
    double cap = Null<Real>(), floor = Null<Real>();
    if (getenv("CAP") != NULL)
        cap = atof(getenv("CAP"));
    if (getenv("FLOOR") != NULL)
        floor = atof(getenv("FLOOR"));

    Date refDate(30, October, 2015);
    Settings::instance().evaluationDate() = refDate;

    Handle<YieldTermStructure> yts1(
        boost::make_shared<FlatForward>(refDate, 0.03, Actual365Fixed()));
    Handle<YieldTermStructure> yts2(
        boost::make_shared<FlatForward>(refDate, 0.04, Actual365Fixed()));

    boost::shared_ptr<SwapIndex> cms2y =
        boost::make_shared<EuriborSwapIsdaFixA>(2 * Years, yts1, yts1);
    boost::shared_ptr<SwapIndex> cms10y =
        boost::make_shared<EuriborSwapIsdaFixA>(10 * Years, yts2, yts2);

    boost::shared_ptr<SwapSpreadIndex> cms10y_2y =
        boost::make_shared<SwapSpreadIndex>("dummy", cms10y, cms2y);

    RelinkableHandle<SwaptionVolatilityStructure> vol;

    if (normal == 0) {
        vol.linkTo(boost::make_shared<ConstantSwaptionVolatility>(
            refDate, TARGET(), ModifiedFollowing,
            Handle<Quote>(boost::make_shared<SimpleQuote>(sigma)),
            Actual365Fixed(), ShiftedLognormal, shift));
    } else {
        vol.linkTo(boost::make_shared<ConstantSwaptionVolatility>(
            refDate, TARGET(), ModifiedFollowing,
            Handle<Quote>(boost::make_shared<SimpleQuote>(sigma)),
            Actual365Fixed(), Normal));
    }

    Date payDate(30, October, 2025);
    Date fixDate(30, October, 2024);

    boost::shared_ptr<CmsCoupon> cms10yCoupon = boost::make_shared<CmsCoupon>(
        payDate, 1.0, fixDate, payDate, 0, cms10y, 1.0);
    boost::shared_ptr<CmsCoupon> cms2yCoupon = boost::make_shared<CmsCoupon>(
        payDate, 1.0, fixDate, payDate, 0, cms2y, 1.0);

    boost::shared_ptr<CappedFlooredCmsSpreadCoupon> spreadCoupon =
        boost::make_shared<CappedFlooredCmsSpreadCoupon>(
            payDate, 1.0, fixDate, payDate, 0, cms10y_2y, 1.0, 0.0, cap, floor);
    boost::shared_ptr<StrippedCappedFlooredCoupon> spreadFloor =
        boost::make_shared<StrippedCappedFlooredCoupon>(spreadCoupon);

    Handle<Quote> reversion(boost::make_shared<SimpleQuote>(0.01));

    boost::shared_ptr<CmsCouponPricer> cmsPricer =
        boost::make_shared<LinearTsrPricer>(
            vol, reversion, Handle<YieldTermStructure>(),
            LinearTsrPricer::Settings().withRateBound(
                normal != 0 ? -2.000 : 0.0001, 2.000));

    Handle<Quote> correlation(boost::make_shared<SimpleQuote>(rho));

    boost::shared_ptr<LognormalCmsSpreadPricer> cmsSpreadPricer =
        boost::make_shared<LognormalCmsSpreadPricer>(
            cmsPricer, correlation, Handle<YieldTermStructure>(), 16,
            boost::none, Null<Real>(), Null<Real>());

    cms2yCoupon->setPricer(cmsPricer);
    cms10yCoupon->setPricer(cmsPricer);

    std::cout << "2y rate " << cms2yCoupon->indexFixing() << " adjusted "
              << cms2yCoupon->adjustedFixing() << std::endl;
    std::cout << "10y rate " << cms10yCoupon->indexFixing() << " adjusted "
              << cms10yCoupon->adjustedFixing() << std::endl;

    std::cout << "sigma=" << sigma << " normal=" << normal << " shift=" << shift
              << " rho=" << rho << " cap=" << cap << " floor=" << floor
              << std::endl;

    spreadFloor->setPricer(cmsSpreadPricer);

    std::cout << "spread coupon rate " << spreadFloor->adjustedFixing()
              << std::endl;

    return 0;
}
