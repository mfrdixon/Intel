#include "vsl.hpp"
#include "svd.hpp"
#include "result.hpp"

#include <tbb/task_scheduler_init.h>
#include <pstl/execution>
#include <pstl/numeric>
#include <mkl.h>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/range/counting_range.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/reversed.hpp>

#include <cmath>
#include <iterator>
#include <vector>
#include <limits>
#include <iomanip>
#include <iostream>
#include <memory>
#include <chrono>

#include <sched.h>

template<typename T>
struct svd_input {
    using vec_t = result::vec_t<T>;

    using generator_type = vsl::gaussian_generator<T>;
    struct matrix_dims {
        using dim_type = int;
        dim_type m_, n_;

        dim_type min() const
        { return std::min(m_, n_); }

        dim_type max() const
        { return std::max(m_, n_); }

        size_t size() const { return m_ * n_; }
    };

    svd_input(vsl::partitioned_stream& strm, generator_type gen, matrix_dims dims)
        : strm_{ strm }
        , gen_{ gen }
        , dims_{ dims }
    { }

    using result_type = std::pair<matrix_dims, vec_t>;
    result_type operator()(int) const {
        vsl::generated_range<generator_type> gr{ strm_.get_next_stream(), gen_ };

        vec_t a;
        a.reserve(dims_.size());
        std::copy_n(std::begin(gr), dims_.size(), std::back_inserter<vec_t>(a));

        return std::make_pair(dims_, a);
    }

private:
    vsl::partitioned_stream& strm_;
    generator_type gen_;

    matrix_dims dims_;
};

template<typename T>
struct svd {
    using vec_t = result::vec_t<T>;
    using res_t = result::res_t<T>;

    using arg_t = typename svd_input<T>::result_type;
    res_t operator()(arg_t const& a) const {
        auto dims = a.first;
        auto aa = const_cast<double*>(a.second.data());

        vec_t s( dims.m_ );
        vec_t superb( dims.min() - 1 );
        vec_t u( dims.m_ * dims.m_ );
        vec_t vt( dims.n_ * dims.n_ );

        auto ec = lapacke::gesvd(LAPACK_COL_MAJOR, 'A', 'A', dims.m_, dims.n_,
                                    aa, dims.m_, s.data(), u.data(), dims.m_,
                                    vt.data(), dims.n_, superb.data());
        if (ec < 0) {
            static auto msg = boost::format("lapacke_gsvd: illegal parameter, ec= %1%");
            throw std::runtime_error(boost::str(msg % ec));
        }
        return res_t{ ec, std::move(s) };
    }
};

template<typename T>
struct merge_svd {
    using result_type = typename svd<T>::res_t;

    result_type operator()(result_type const& a, result_type const& b) {
        return merge(a, b);
    }
};

template<typename T>
double dtime(T t) {
    using dseconds = std::chrono::duration<double, std::ratio<1>>;
    return std::chrono::duration_cast<dseconds>(t).count();
}

using res_t = std::tuple<double, double, double>;
res_t do_it(int M, unsigned sample_ct) {
    auto const dt = sqrt(1.0/ M);

    auto const m = 400;
    auto const n = 250;

    using dsvd_input = svd_input<double>;
    using dsvd = svd<double>;
    using generator_type = vsl::gaussian_generator<double>;

    using clock = std::chrono::high_resolution_clock;
    std::vector<clock::duration> samples;
    for (auto i = 0; i < sample_ct; ++i) {
        vsl::partitioned_stream strm{
                vsl::stream{ vsl::stream::rng_kind::sfmt19937, 1 }
            };

        auto irange = boost::counting_range( 0, M )
                      | boost::adaptors::transformed(dsvd_input{ strm, generator_type{ 0, dt },
                                                                 dsvd_input::matrix_dims{ m, n } });

        std::vector<dsvd::arg_t> input;
        std::copy(std::begin(irange), std::end(irange), std::back_inserter(input));

        auto t0 = clock::now();
        auto res = std::transform_reduce(std::execution::par, std::begin(input), std::end(input),
                                            dsvd::res_t{}, merge_svd<double>{ },dsvd{});

        auto t1 = clock::now();

        samples.push_back(t1 - t0);
    }

    if (sample_ct > 1) {
        std::sort(std::begin(samples), std::end(samples));
        return std::make_tuple(dtime(samples.front()),
                               dtime(samples.at(sample_ct/2)),
                               dtime(samples.back()));
    }
    auto r = dtime(samples.front());
    return std::make_tuple(r, r, r);
}

int main(int argc, char const* argv[]) {
    // Knights Landing topology
    enum { core_ct = 68,
           threads_per_core = 4 };
    using tpc_t = std::array<int, threads_per_core>; 
    std::array<tpc_t, core_ct> core_assignment;
    for (auto i = 0; i < core_assignment.size(); ++i) {
        core_assignment[i] = { i, i+core_ct, i+2*core_ct, i+3*core_ct};
    }

    auto M = 10000;
    if (argc > 1) {
        M = boost::lexical_cast<decltype(M)>(argv[1]);
        std::cout << "Evaluating " << M << " parallel SVD calls" << std::endl;
    }

    auto cores = 68;
    auto threads = 4;
    if (argc == 4) {
        cores = boost::lexical_cast<decltype(cores)>(argv[2]);
        threads = boost::lexical_cast<decltype(threads)>(argv[3]);
    }

    auto sample_ct = 1u;

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    for (auto i = 0; i < cores; ++i) {
        for (auto j = 0; j < threads; ++j) {
            CPU_SET(core_assignment[i][j], &cpuset);
        }
    }
    auto ec = sched_setaffinity(0, sizeof(cpuset), &cpuset);
    if (ec < 0) {
        std::cerr << "sched_setaffinity returned ec=" << ec << std::endl;
        return ec;
    }
    
    auto const cpu_count = CPU_COUNT(&cpuset);
    tbb::task_scheduler_init init(cpu_count);
    std::cout << "sampling cores=" << cores << ", threads per core=" << threads << "....";
    std::cout.flush();
    double min, med, max;
    std::tie(min, med, max) = do_it(M, sample_ct);
    std::cout << "done (" << med << ")." << std::endl;

    std::cout << "NT, TMIN, TMED, TMAX\n" << std::setprecision(4)
              << threads << ", "
              << min << ", "
              << med << ", "
              << max << std::endl;
    return 0;
}
