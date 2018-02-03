#pragma once

#include <mkl_vsl.h>

#include <boost/assert.hpp>
#include <boost/optional.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include <atomic>
#include <mutex>
#include <memory>
#include <vector>
#include <deque>
#include <unordered_map>
#include <thread>
#include <stdexcept>
#include <sstream>

namespace vsl {
    struct generator_core_access;

    struct stream_base {
        stream_base(stream_base const& other)
            : strm_{ copy_stream(other.strm_) }
        { }

        stream_base& operator=(stream_base const& other) {
            strm_ = copy_stream(other.strm_);
            return *this;
        }

        bool valid() const { return strm_.get() != nullptr; }

    protected:
        friend generator_core_access;

        struct deleter {
            void operator()(VSLStreamStatePtr p) {
                vslDeleteStream(&p);
            }
        };

        using ptr_t = std::unique_ptr<void, deleter>;
        ptr_t::pointer get_stream() { return strm_.get(); }

        stream_base(ptr_t strm)
            : strm_{ std::move(strm) }
        { }

    private:
        ptr_t strm_;

        static ptr_t copy_stream(ptr_t const& other) {
            ptr_t::pointer res;
            auto ec = vslCopyStream(&res, other.get());
            if (ec != VSL_STATUS_OK) {
                std::ostringstream stm;
                stm << "Error copying VSL stream, ec=" << ec;
                throw std::runtime_error(stm.str());
            }
            return ptr_t{ res };
        }
    };

    struct stream : stream_base {
        enum class rng_kind : MKL_INT {
            sfmt19937 = VSL_BRNG_SFMT19937
        };

        stream(rng_kind brng, MKL_UINT seed)
            : stream_base{ make_stream(brng, seed) }
        { }

    private:
        static ptr_t make_stream(rng_kind brng, MKL_UINT seed) {
            ptr_t::pointer res;
            auto ec = vslNewStream(&res, static_cast<MKL_INT>(brng), seed);
            if (ec != VSL_STATUS_OK) {
                std::ostringstream stm;
                stm << "Error initializing VSL stream, ec=" << ec;
                throw std::runtime_error(stm.str());
            }
            return ptr_t{ res };
        }
    };

    struct skip_ahead_stream : stream_base {
        skip_ahead_stream(stream_base const& strm, long long skip)
            : stream_base{ strm }
        { make_skip_ahead(get_stream(), skip); }

    private:
        static void make_skip_ahead(ptr_t::pointer const& strm, long long skip) {
            auto ec = vslSkipAheadStream(strm, skip);
            if (ec != VSL_STATUS_OK) {
                std::ostringstream stm;
                stm << "Error initializing VSL stream, ec=" << ec;
                throw std::runtime_error(stm.str());
            }
        }
    };

    struct generator_core_access {
        generator_core_access(stream_base& strm)
            : strm_{ strm }
        { BOOST_ASSERT(strm_.valid()); }

        using pointer_t = stream::ptr_t::pointer;
        pointer_t get_stream() { return strm_.get_stream(); }

    private:
        stream_base& strm_;
    };

    struct partitioned_stream : stream_base {
        partitioned_stream(stream_base strm, unsigned num_shards = 64)
            : stream_base{ std::move(strm) }
            , pimpl_{ std::make_shared<shared_streams>(*this, num_shards) }
        { }

        static constexpr short max_stream_partitions() { return std::numeric_limits<short>::max(); }
        static constexpr size_t partition() { return (1ull << 48) - 1; }

        skip_ahead_stream get_next_stream() {
            BOOST_ASSERT(pimpl_);
            auto res = pimpl_->get_next_stream();
            return res;
        }

    private:
        std::atomic<short> next_;

        struct shard_t {
            stream_base strm_;
            std::mutex mtx_;

            shard_t(stream_base const& strm)
                : strm_{ strm }
            { }
        };

        struct shared_streams {
            std::atomic<short> next_;
            std::deque<shard_t> shards_;

            using lock_type = std::unique_lock<decltype(shard_t::mtx_)>;

            shared_streams(stream_base const& strm, unsigned num_shards)
                : next_{ 0 }
            {
                while (--num_shards)
                    shards_.emplace_back(strm);
            }

            skip_ahead_stream get_next_stream() {
                static std::hash<std::thread::id> hash;
                auto tid = std::this_thread::get_id();
                auto& shard = shards_[hash(tid) % shards_.size()];
                lock_type l(shard.mtx_);
                auto skip = next_++;
                if (next_ < 0)
                    throw std::runtime_error("to many partitions");
                return skip_ahead_stream{ shard.strm_, static_cast<long long>(skip) << 48 };
            }
        };
        std::shared_ptr<shared_streams> pimpl_;
    };

    namespace detail {
        template<typename Q>
        int rng_gaussian(int, void*, size_t, Q*, Q, Q);

        template<>
        int rng_gaussian<double>(int mthd, void* strm, size_t N,
                                    double* res, double a, double sigma) {
            BOOST_ASSERT(strm != nullptr);
            BOOST_ASSERT(res != nullptr);
            return vdRngGaussian(mthd, strm, N, res, a, sigma);
        }

        template<>
        int rng_gaussian<float>(int mthd, void* strm, size_t N,
                                            float* res, float a, float sigma) {
            BOOST_ASSERT(strm != nullptr);
            BOOST_ASSERT(res != nullptr);
            return vsRngGaussian(mthd, strm, N, res, a, sigma);
        }
    } // namespace detail 

    template<typename T>
    struct gaussian_generator {
        using value_type = T;

        enum class method : int {
            boxmuller2 = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2
        };

        gaussian_generator(value_type a, value_type sigma,
                           method mthd = method::boxmuller2)
            : a_{ a }
            , sigma_{ sigma }
            , mthd_{ mthd }
        { }

        template<typename Stream, typename ContiguousIterator>
        void operator()(Stream& strm, ContiguousIterator begin, ContiguousIterator end) const {
            BOOST_ASSERT(strm.valid());
            generator_core_access s{ strm };
            auto ss = s.get_stream();
            BOOST_ASSERT(ss != nullptr);

            auto ec = detail::rng_gaussian(static_cast<int>(mthd_), s.get_stream(),
                                                std::distance(begin, end), &(*begin), a_, sigma_);
            if (ec != VSL_STATUS_OK) {
                std::ostringstream stm;
                stm << "Error generating VSL data, ec=" << ec;
                throw std::runtime_error(stm.str());
            }
        }

    private:
        value_type a_, sigma_;
        method mthd_;
    };

    template<typename Generator>
    struct generated_range {
        using value_type = typename Generator::value_type;

        generated_range(skip_ahead_stream strm, Generator gen, unsigned buf_max = 8192)
            : state_{ std::make_shared<state>(std::move(strm), std::move(gen), buf_max) }
        { }

    private:
        using vec_t = std::vector<value_type>;

        struct state {
            skip_ahead_stream strm_;
            Generator gen_;
            vec_t buf_;
            unsigned last_;

            state(skip_ahead_stream strm, Generator gen, unsigned buf_max)
                : strm_{ std::move(strm) }
                , gen_{ gen }
                , buf_( buf_max )
                , last_( buf_max ) // <- Unicorn Initialization
            { }

            value_type next(long long& ix) {
                auto cur = ++ix % buf_.size();
                if (cur < last_)
                    gen_(strm_, std::begin(buf_), std::end(buf_));
                last_ = cur;
                return buf_[cur];
            }
        };
        using ptr_t = std::shared_ptr<state>;
        ptr_t state_;

        struct const_iterator
            : boost::iterator_facade<
                  const_iterator
                , value_type const
                , boost::forward_traversal_tag
              > {

            const_iterator(ptr_t s, long long ix)
                : s_{ std::move(s) }
                , ix_{ ix }
            { }

        private:
            friend boost::iterator_core_access;
            ptr_t s_;
            long long ix_;
            value_type val_;

            void increment() { val_ = s_->next(ix_); }

            bool equal(const_iterator const& other) const
            { return s_ == other.s_ && ix_ == other.ix_; }

            value_type const& dereference() const { val_; }
        };

    public:
        const_iterator begin() const {
            return const_iterator{ state_, 0 };
        }

        const_iterator end() const {
            return const_iterator{ state_, partitioned_stream::partition() };
        }
    };
}
