/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Klaus Spanderen

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the QuantLib license for more details.
*/

#ifndef quantlib_test_mpi_parallel_test_runner_hpp
#define quantlib_test_mpi_parallel_test_runner_hpp

#include <ql/types.hpp>

#include <boost/timer.hpp>
#include <boost/thread/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/serialization/utility.hpp>

#define BOOST_TEST_NO_MAIN 1
#include <boost/test/results_collector.hpp>

#include <boost/test/included/unit_test.hpp>
#include <boost/algorithm/string.hpp>

#include <boost/mpi.hpp>

#include <map>
#include <list>
#include <sstream>
#include <utility>
#include <fstream>
#include <iostream>

#include <string>
#include <cstring>
#include <cstdlib>

#ifdef BOOST_MSVC
#  error parallel test suite runner is not available on Windows
#endif

using boost::unit_test::test_results;
using namespace boost::unit_test_framework;


namespace {
    class TestCaseCollector : public test_tree_visitor {
      public:
        typedef std::map<test_unit_id, std::list<test_unit_id> > id_map_t;

        const id_map_t& map() const { return idMap_; }
        test_unit_id testSuiteId() const { return testSuiteId_; }

        bool visit(test_unit const& tu) {
            if (tu.p_parent_id == framework::master_test_suite().p_id) {
                QL_REQUIRE(!tu.p_name.get().compare("QuantLib test suite"),
                     "could not find QuantLib test suite");
                testSuiteId_ = tu.p_id;
            }
            return test_tree_visitor::visit(tu);
        }

        void visit(test_case const& tc) {
            idMap_[tc.p_parent_id].push_back(tc.p_id);
        }

        std::list<test_unit_id>::size_type numberOfTests() {
            std::list<test_unit_id>::size_type n=0;
            for (id_map_t::const_iterator p_it = idMap_.begin();
                p_it != idMap_.end(); ++p_it) n+=p_it->second.size();

            return n;
        }
      private:
        id_map_t idMap_;
        test_unit_id testSuiteId_;
    };

    struct TestCaseId {
        friend class boost::serialization::access;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version) {
            ar & id;
            ar & terminate;
        }

        test_unit_id id;
        bool terminate;
    };

    const char* const namesLogMutexName = "named_log_mutex";

    std::ostream& log_message(std::ostream& out, const std::string& msg) {
        std::vector<std::string> tok;
        boost::split(tok, msg, boost::is_any_of("\n"));

        for (std::vector<std::string>::const_iterator
            iter = tok.begin(); iter != tok.end(); ++iter) {
            if (iter->length() && iter->compare("Running 1 test case...")) {
                out << *iter << std::endl;
            }
        }

        return out;
    }

    void output_logstream(
        std::ostream& out, std::streambuf* outBuf, std::stringstream& s) {

        out.flush();
        out.rdbuf(outBuf);

        log_message(out, s.str());

        s.str(std::string());
        out.rdbuf(s.rdbuf());
    }
}

test_suite* init_unit_test_suite(int, char* []);

int main( int argc, char* argv[] )
{
    typedef QuantLib::Time Time;
    const char* const profileFileName = ".unit_test_profile.txt";

    typedef std::vector<std::pair<std::string, QuantLib::Time> >
        run_time_list_type;

    boost::mpi::environment env;
    boost::mpi::communicator world;

    if (!world.rank()) {
        framework::init(init_unit_test_suite, argc, argv );
        framework::finalize_setup_phase();

        std::map<std::string, Time> runTimeLog;

        std::ifstream in(profileFileName);
        if (in.good()) {
            for (std::string line; std::getline(in, line);) {
                std::vector<std::string> tok;
                boost::split(tok, line, boost::is_any_of(" "));

                QL_REQUIRE(tok.size() == 2,
                    "every line should consists of two entries");
                runTimeLog[tok[0]] = std::atof(tok[1].c_str());
            }
        }
        in.close();

        TestCaseCollector tcc;
        traverse_test_tree(framework::master_test_suite(), tcc , true);

        s_log_impl().stream() << "Total number of test cases: "
            << tcc.numberOfTests() << std::endl;

        // run root test cases in master process
        const std::list<test_unit_id>& qlRoot
            = tcc.map().find(tcc.testSuiteId())->second;

        std::multimap<Time, test_unit_id> testsSortedByRunTime;

        for (TestCaseCollector::id_map_t::const_iterator
            p_it = tcc.map().begin();
            p_it != tcc.map().end(); ++p_it) {

            if (p_it->first != tcc.testSuiteId()) {
                for (std::list<test_unit_id>::const_iterator
                    it =  p_it->second.begin();
                    it != p_it->second.end(); ++it) {

                    const std::string& name
                        = framework::get(*it, TUT_ANY).p_name;

                    if (runTimeLog.count(name)) {
                        testsSortedByRunTime.insert(
                            std::make_pair(runTimeLog[name], *it));
                    }
                    else {
                        testsSortedByRunTime.insert(
                            std::make_pair(
                                std::numeric_limits<Time>::max(), *it));
                    }
                }
            }
        }

        std::vector<test_unit_id> ids;
        for (std::multimap<Time, test_unit_id>::const_iterator
            iter = testsSortedByRunTime.begin();
            iter != testsSortedByRunTime.end(); ++iter) {
            ids.push_back(iter->second);
        }
        QL_REQUIRE(ids.size() + qlRoot.size() == tcc.numberOfTests(),
            "missing test case in distrubtion list");

        testsSortedByRunTime.clear();

        test_results results;

        std::stringstream logBuf;
        std::streambuf* const oldBuf = s_log_impl().stream().rdbuf();
        s_log_impl().stream().rdbuf(logBuf.rdbuf());

        for (std::list<test_unit_id>::const_iterator iter = qlRoot.begin();
            std::distance(qlRoot.begin(), iter) < int(qlRoot.size())-1;
            ++iter) {
            framework::run(*iter);
            results += boost::unit_test::results_collector.results(*iter);
        }
        output_logstream(s_log_impl().stream(), oldBuf, logBuf);
        s_log_impl().stream().rdbuf(oldBuf);

        const unsigned nSlaves =
            std::min(size_t(world.size()-1), ids.size());

        for (unsigned i=nSlaves+1; i < world.size(); ++i) {
            const TestCaseId id = { 0, true };
            world.send(i, 0, id);
        }

        for (unsigned i=0; i < ids.size() + nSlaves; ++i) {
            if (i < nSlaves) {
                const TestCaseId id = { ids[ids.size()-i-1], false };
                world.send(i+1, 0, id);
            }
            else {
                std::string msg;
//                while(!world.iprobe(boost::mpi::any_source, 2))
//                    boost::this_thread::sleep(boost::posix_time::microsec(1000));

                const boost::mpi::status s =
                        world.recv(boost::mpi::any_source, 2, msg);

                if (i < ids.size()) {
                    const TestCaseId id = { ids[ids.size()-i-1], false };
                    world.send(s.source(), 0, id);
                }
                else {
                    const TestCaseId id = { 0, true };
                    world.send(s.source(), 0, id);
                }

                log_message(std::cout, msg);
            }
        }

        run_time_list_type runTimeLogs;

        for (unsigned i=0; i < world.size()-1; ++i) {
            run_time_list_type log;
            world.recv(boost::mpi::any_source, 3, log);

            for (run_time_list_type::const_iterator
                iter = log.begin(); iter != log.end(); ++iter) {
                runTimeLog[iter->first] = iter->second;
            }

            test_results r;
            world.recv(boost::mpi::any_source, 1,
                (char*) &r, sizeof(test_results));
            results+=r;
        }

        if (!qlRoot.empty())
            framework::run(qlRoot.back());

        results += boost::unit_test::results_collector.results(
            qlRoot.back());

        s_log_impl().stream() << "Test module \""
            << framework::master_test_suite().p_name
            <<"\" has passed with:"
            << std::endl
            << " " << results.p_assertions_failed << " assertions failed"
            << std::endl
            << " " << results.p_assertions_passed << " assertions passed"
            << std::endl;

        std::ofstream out(
            profileFileName, std::ios::out | std::ios::trunc);
        out << std::setprecision(6);
        for (std::map<std::string, QuantLib::Time>::const_iterator
            iter = runTimeLog.begin(); iter != runTimeLog.end(); ++iter) {
            out << iter->first << " " << iter->second << std::endl;
        }
        out.close();
    }
    else {
        std::stringstream logBuf;
        s_log_impl().stream().rdbuf(logBuf.rdbuf());

        framework::init(init_unit_test_suite, argc, argv );
        framework::finalize_setup_phase();

        logBuf.str(std::string());

        test_results r;
        run_time_list_type runTimeLogs;

        TestCaseId id;
        world.recv(0, 0, id);

        while (!id.terminate) {

            boost::timer t;

            framework::run(id.id);
            r += boost::unit_test::results_collector.results(id.id);

            runTimeLogs.push_back(std::make_pair(
                framework::get(id.id, TUT_ANY).p_name, t.elapsed()));

            s_log_impl().stream().flush();
            const std::string log = logBuf.str();
            logBuf.str(std::string());

            world.send(0, 2, log);
            world.recv(0, 0, id);
        }

        world.send(0, 3, runTimeLogs);
//        boost::mpi::request t =
                world.send(0, 1, (char*) &r, sizeof(test_results));
//        while(!boost::mpi::test_all(&t, &t+1))
//            boost::this_thread::sleep(boost::posix_time::microsec(1000));
    }

    framework::shutdown();
}

#endif
