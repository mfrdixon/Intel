#pragma once

#include <boost/optional.hpp>

#include <vector>
#include <memory>
#include <tuple>

namespace result {
    template<typename T>
    using vec_t = std::vector<T>;

    template<typename T>
    struct res_t {
        res_t() = default;

        res_t(int ec, vec_t<T> v)
            : pimpl_{ std::make_shared<impl>(ec, std::move(v)) }
        { }

        void swap(res_t& other) { std::swap(pimpl_, other.pimpl_); }

        size_t size() const {
            if (!pimpl_)
                return 0;
            return pimpl_->size();
        }

        friend res_t merge(res_t lhs, res_t rhs) {
            if (lhs.size() == 0)
                return std::move(rhs);
            if (rhs.size()) {
                // always merge into larger of the two
                if (lhs.size() < rhs.size())
                    std::swap(lhs, rhs);

                lhs.pimpl_->merge(*rhs.pimpl_);
            }
            return std::move(lhs);
        }

    private:
        // this is hideous, but Thrust insists on copies, we want to
        // make them cheap
        using el_t = std::tuple<int, vec_t<T>>;
        using opt_t = boost::optional<el_t>;
        struct impl {
            opt_t v_;
            std::vector<el_t> vs_;

            impl(int ec, vec_t<T> v)
                : v_{ std::make_pair(ec, std::move(v)) }
            { }

            size_t size() const {
                if (v_)
                    return 1;
                return vs_.size();
            }

            void emplace_back(opt_t& v) {
                BOOST_ASSERT(static_cast<bool>(v));
                vs_.emplace_back(std::make_pair(std::get<0>(*v),
                                 std::move(std::get<1>(*v))));
                v = boost::none;
            }

            void merge(impl& other) {
                if (!size())
                    return; // empty

                BOOST_ASSERT(size() >= other.size());
                vs_.reserve(size() + other.size());

                if (vs_.empty()) {
                    emplace_back(v_);
                }

                if (other.v_) {
                    emplace_back(other.v_);
                } else {
                    for (auto&& el : other.vs_)
                        vs_.emplace_back(std::move(el));
                }
            }
        };
        std::shared_ptr<impl> pimpl_;
    };
}
