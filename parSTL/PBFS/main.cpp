#include <boost/assert.hpp>
#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <vector>
#include <memory>
#include <functional>
#include <string>
#include <atomic>

namespace detail {
struct pennant {
    pennant() { reset(); }
    virtual void reset() {
        left_ = nullptr;
        right_ = nullptr;
    }

    static pennant& merge(pennant& x, pennant& y) {
        y.right_ = x.left_;
        x.left_ = &y;
        return x;
    }

    static pennant& split(pennant& x) {
        // must contain at least two elements
        BOOST_ASSERT(x.left_ != nullptr);
        BOOST_ASSERT(x.right_ != nullptr);

        auto y = x.left_;
        x.left_ = y->right_;
        y->right_ = nullptr;
        return *y;
    }

    size_t size() const {
        auto res = left_ ? left_->size() + 1 : 1;
        return res + right_ ? right_->size() : 0;
    }

    template<typename Function>
    void visit(Function&& f) {
        f(*this);
        if (left_)
            left_->visit(std::forward<Function>(f));
        if (right_)
            right_->visit(std::forward<Function>(f));
    }

private:
    pennant* left_;
    pennant* right_;
};

struct bag {
    bag(long r = 1024)
        : spine_(r, nullptr)
    { }

    void swap(bag& rhs) {
        std::swap(spine_, rhs.spine_);
    }

    void insert(pennant& x) {
        auto k = 0;
        while (spine_[k]) {
            x = pennant::merge(*spine_[k], x);
            spine_[k++] = nullptr;
        }
        spine_[k] = &x;
    }

    bool empty() const {
        return std::all_of(std::begin(spine_), std::end(spine_),
                           [](auto const& x) { return x == nullptr; });
    }

    size_t size() const { return spine_.size(); }
    pennant* operator[](size_t k) { return spine_[k]; }
    pennant const* operator[](size_t k) const { return spine_[k]; }

    static bag& merge(bag& l, bag& r) {
        // unfortunately these are necessary without transform_reduce
        if (l.spine_.size() < r.spine_.size())
            l.spine_.resize(r.spine_.size(), nullptr);
        else if (r.spine_.size() < l.spine_.size())
            r.spine_.resize(l.spine_.size(), nullptr);

        pennant* y = nullptr;
        for (auto k = 0; k < l.spine_.size(); ++k) {
            pennant* sk;
            std::tie(sk, y) = fa(l.spine_[k], r.spine_[k], y);
            l.spine_[k] = sk;
        }
        return l;
    }

private:
    std::vector<pennant*> spine_;

    using res_t = std::pair<pennant*, pennant*>;
    static res_t fa(pennant* x, pennant* y, pennant* z) {
        if (!(x || y || z))
            return std::make_pair(nullptr, nullptr);
        if (x && !(y || z))
            return std::make_pair(x, nullptr);
        if (y && !(x || z))
            return std::make_pair(y, nullptr);
        if (z && !(x || y))
            return std::make_pair(z, nullptr);
        if (x && z && !y)
            return std::make_pair(nullptr, &pennant::merge(*x, *z));
        if (y && z && !x)
            return std::make_pair(nullptr, &pennant::merge(*y, *z));
        return std::make_pair(x, &pennant::merge(*y, *z));
    }
};
} // namespace detail

template<typename T, typename E>
struct vertex : detail::pennant {
    using this_type = vertex<T, E>;
    using value_type = T;
    using name_type = E;

    this_type* up_ = nullptr;
    std::atomic<long> distance_;

    using optional_t = boost::optional<value_type>;
    optional_t val_;

    struct identity_fn {
        optional_t operator()(T const& t) const { return t; }
    };

    struct filter_fn {
        using type = std::function<value_type(value_type const&)>;
        type f_;
        filter_fn(type&& f) : f_{ std::forward<type>(f) } { }

        optional_t operator()(T const& t) const { return f_(t); } 
    };

    struct sink_fn {
        using type = std::function<void(value_type const&)>;

        type f_;
        sink_fn(type&& f) : f_{ std::forward<type>(f) } { }

        optional_t operator()(T const& t) const {
            f_(t);
            return boost::none;
        }
    };

    boost::variant<identity_fn, filter_fn, sink_fn> filter_;
    using pipes_t = std::vector<name_type>;
    pipes_t pipes;

    vertex()
        : filter_{ identity_fn{ } }
        , distance_{ -1 }
    { }

    vertex(vertex const& other)
        : filter_{ other.filter_ }
        , distance_{ other.distance_.load() }
    { }

    vertex& operator=(vertex const& other) {
        filter_ = other.filter_;
        distance_ = other.distance_;
        return *this;
    }

    vertex(filter_fn f)
        : filter_{ std::move(f) }
        , distance_{ -1 }
    { }

    vertex(sink_fn f)
        : filter_{ std::move(f) }
        , distance_{ -1 }
    { } 

    optional_t apply() const {
        if (!val_)
            return boost::none;

        apply_fn_visitor v{ val_ };
        return boost::apply_visitor(v, filter_);
    }

    bool reached() const { return distance_ != -1; }
    void reset() override {
        pennant::reset();
        distance_ = -1;
        val_ = boost::none;
    }

private:
    struct apply_fn_visitor : boost::static_visitor<optional_t> {
        optional_t opt;
        apply_fn_visitor(optional_t o)
            : opt{ o }
        { }

        template<typename F>
        optional_t operator()(F const& f) { return f(opt); }
    };
};

template<typename T, typename E>
struct graph {
    using vertex_type = vertex<T, E>;
    using value_type = typename vertex_type::value_type;
    using name_type = E;

    bool add_source(name_type&& name) {
        auto [_, b] = vertices_.emplace(std::make_pair(std::forward<name_type>(name),
                                                       vertex_type{ }));
        return b;
    }

    bool add_filter(name_type&& name, std::function<value_type(value_type)> f) {
        using fn_t = typename vertex_type::filter_fn;
        auto [_, b] = vertices_.emplace(std::make_pair(std::forward<name_type>(name),
                                                       vertex_type{ fn_t{ std::move(f) } }));
        return b;
    }

    bool add_sink(name_type&& name, std::function<void(value_type)> f) {
        using fn_t = typename vertex_type::sink_fn;
        auto [_, b] = vertices_.emplace(std::make_pair(std::forward<name_type>(name),
                                                       vertex_type{ fn_t{ std::move(f) } }));
        return b;
    }

    bool add_pipe(name_type const& from, name_type const& to) {
        auto itf = vertices_.find(from);
        if (itf == std::end(vertices_))
            return false;

        auto itt = vertices_.find(to);
        if (itt == std::end(vertices_))
            return false;

        if (up_find(&itf->second) == up_find(&itt->second))
            return false; // edge would violate the DAG

        itt->second.up_ = &itf->second;
        itf->second.pipes.emplace_back(to);
        return true;
    }

    vertex_type& operator[](name_type const& name)
    { return vertices_[name]; }

    vertex_type const& operator[](name_type const& name) const
    { return vertices_[name]; }

    void reset() {
        std::for_each(std::begin(vertices_), std::end(vertices_),
                      [](auto& v) { v.reset(); });
    }

private:
    using map_type = std::unordered_map<name_type, vertex_type>;
    map_type vertices_;

    vertex_type* up_find(vertex_type* v) {
        auto last = v;
        auto up = v->up_;
        while (up) {
            v->up_ = up->up_; // compress path
            last = v;
            v = up;
            up = v->up_;
        }
        return last;
    }

public:
    using iterator = typename map_type::iterator;
    using const_iterator = typename map_type::const_iterator;

    iterator begin() { return vertices_.begin(); }
    iterator end() { return vertices_.end(); }
    const_iterator begin() const { return vertices_.begin(); }
    const_iterator end() const { return vertices_.end(); }
};

namespace detail {
template<typename ExecutionPolicy,
         typename Graph>
bag process_pennant(ExecutionPolicy&& exec, Graph const& g, pennant* in, long d, size_t grain_size = 128) {
    bag res;
    using vertex_type = typename Graph::vertex_type;
    // if (in->size() < grain_size) {
        in->visit([&](auto& u) {
                auto v = reinterpret_cast<vertex_type&>(u);
                auto val = v.apply();
                std::for_each(std::forward<ExecutionPolicy>(exec), std::begin(v.pipes), std::end(v.pipes),
                              [&](auto const& n) {
                                    auto& adjv = g[n];
                                    long expected = -1;
                                    if (std::atomic_compare_exchange_strong(adjv.distance_, expected, d + 1)) {
                                        adjv.val_ = val;
                                        res.insert(adjv);
                                    }
                                });
            });
        return res;
    // }
    // TODO implement task spawning to split large inputs
}

template<typename ExecutionPolicy,
         typename Graph>
bag process_layer(ExecutionPolicy exec, Graph const& g, bag& in, long d) {
    return std::transform_reduce(std::forward<ExecutionPolicy>(exec), 0, in.size(), bag{}
                                 bag::merge
                                 [&](auto k) {
                                        if (in[k] != nullptr) {
                                            auto out = process_pennant(std::forward<ExecutionPolicy>(exec), g, in[k], d);
                                            bag::merge(res, out);
                                        }
                                 });
}
} // detail

template<typename ExecutionPolicy,
         typename Graph>
void post(ExecutionPolicy&& exec, Graph const& g, typename Graph::name_type const& v0, typename Graph::value_type val) {
    std::for_each(std::forward<ExecutionPolicy>(exec), std::begin(g), std::end(g),
                  [](auto& v) { v.reset(); });

    auto d = 0;
    auto& v = g[v0];
    v.distance_ = d;
    v.val_ = val;

    detail::bag bag{};
    bag.insert(v);
    while (!bag.empty()) {
        auto out = detail::process_layer(std::forward<ExecutionPolicy>(exec), bag, d);
        std::swap(bag, out);
        ++d;
    }
}

int main(int, const char* argv[]) {
    graph<int, int> g;
    g.add_source(0);
    g.add_filter(1, [](int x) { return x * 2; });
    g.add_filter(2, [](int x) { return x * 3; });
    g.add_sink(3, [](int x) { (void) x; });

    g.add_pipe(0, 1);
    g.add_pipe(0, 2);
    g.add_pipe(1, 3);
    g.add_pipe(2, 3);

    post(std::execution::seq, g, 0, 42);
    return 0;
}
