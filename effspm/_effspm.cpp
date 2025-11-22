// _effspm.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>
#include <fstream>
#include <cstdio>   // std::remove
#include <vector>
#include <string>
#include <ctime>
#include <cmath>

namespace py = pybind11;

// PrefixProjection headers (global namespace)
#include "freq_miner.hpp"
#include "load_inst.hpp"
#include "utility.hpp"

// BTMiner (namespaced)
#include "btminer/src/freq_miner.hpp"
#include "btminer/src/load_inst.hpp"
#include "btminer/src/utility.hpp"
#include "btminer/src/build_mdd.hpp"

// HTMiner (namespaced)
#include "htminer/src/build_mdd.hpp"
#include "htminer/src/freq_miner.hpp"
#include "htminer/src/load_inst.hpp"
#include "htminer/src/utility.hpp"

// LargePrefixProjection
#include "largepp/src/freq_miner.hpp"
#include "largepp/src/load_inst.hpp"
#include "largepp/src/utility.hpp"

// LargeBTMiner
#include "largebm/src/freq_miner.hpp"
#include "largebm/src/load_inst.hpp"
#include "largebm/src/utility.hpp"
#include "largebm/src/build_mdd.hpp"

// LargeHTMiner
#include "largehm/src/freq_miner.hpp"
#include "largehm/src/load_inst.hpp"
#include "largehm/src/utility.hpp"
#include "largehm/src/build_mdd.hpp"

namespace {

// RAII helper for temp file
struct TempFile {
    std::string path;
    ~TempFile() {
        if (!path.empty()) {
            std::remove(path.c_str());
        }
    }
};

// Write Python list[list[int]] to a temp file in professor’s format:
// one sequence per line, items separated by spaces.
std::string write_temp_seq_file(const std::vector<std::vector<int>>& seqs) {
    char tmp_name[L_tmpnam];
    if (!std::tmpnam(tmp_name)) {
        throw std::runtime_error("Failed to create temporary file name");
    }
    std::string path = std::string(tmp_name) + ".txt";

    std::ofstream ofs(path);
    if (!ofs) {
        throw std::runtime_error("Failed to open temporary file for writing: " + path);
    }

    for (const auto& seq : seqs) {
        for (size_t i = 0; i < seq.size(); ++i) {
            if (i) ofs << ' ';
            ofs << seq[i];
        }
        ofs << '\n';
    }

    ofs.close();
    return path;
}

} // anonymous namespace


PYBIND11_MODULE(_effspm, m) {
    m.doc() = "Unified SPM library: PrefixProjection, BTMiner, HTMiner, Large* variants";

    // ─────────────────────────────────────────────────────────────
    // PrefixProjection  (works directly on Python lists or files)
    // ─────────────────────────────────────────────────────────────
    m.def("PrefixProjection",
        [](py::object data,
           double minsup,
           unsigned int time_limit,
           bool preproc,
           bool use_dic,
           bool verbose,
           const std::string &out_file)
        {
            ::time_limit = time_limit;
            ::pre_pro    = preproc;
            ::use_dic    = use_dic;
            ::use_list   = false;
            ::b_disp     = verbose;                  // controls prints in original code
            ::b_write    = !out_file.empty();
            ::out_file   = out_file;

            ClearCollected();
            start_time = std::clock();

            if (py::isinstance<py::str>(data)) {
                std::string path = data.cast<std::string>();
                if (!Load_instance(path, minsup))
                    throw std::runtime_error("PrefixProjection: failed to load file: " + path);
            } else {
                auto seqs = data.cast<std::vector<std::vector<int>>>();
                items = std::move(seqs);
                N = items.size();

                int max_id = 0;
                for (auto &seq : items)
                    for (int x : seq)
                        max_id = std::max(max_id, std::abs(x));
                L = max_id;

                theta = (minsup < 1.0) ? std::ceil(minsup * N) : minsup;

                DFS.clear();
                DFS.reserve(L);
                for (unsigned int i = 0; i < L; ++i)
                    DFS.emplace_back(-static_cast<int>(i) - 1);

                M = 0;
                E = 0;
                for (auto &seq : items) {
                    M = std::max<unsigned int>(M, seq.size());
                    E += seq.size();
                }
            }

            Freq_miner();

            py::dict out;
            out["patterns"] = GetCollected();
            out["time"]     = give_time(std::clock() - start_time);
            return out;
        },
        py::arg("data"),
        py::arg("minsup")     = 0.01,
        py::arg("time_limit") = 36000,
        py::arg("preproc")    = false,
        py::arg("use_dic")    = false,
        py::arg("verbose")    = false,
        py::arg("out_file")   = ""
    );

    // ─────────────────────────────────────────────────────────────
    // BTMiner  (always uses professor's Load_instance)
    // ─────────────────────────────────────────────────────────────
   // ─────────────────────────────────────────────────────────────
// BTMiner  (always uses professor's Load_instance)
// ─────────────────────────────────────────────────────────────
/*m.def("BTMiner",
    [](py::object data,
       double minsup,
       unsigned int time_limit,
       bool preproc,
       bool use_dic,
       bool verbose,
       const std::string &out_file)
    {
        // 1) Configure professor globals
        btminer::time_limit = static_cast<int>(time_limit);
        btminer::pre_pro    = preproc;
        btminer::use_dic    = use_dic;
        btminer::b_disp     = verbose;
        btminer::b_write    = !out_file.empty();
        btminer::out_file   = out_file;
        btminer::N_mult     = 1;
        btminer::M_mult     = 1;
        btminer::just_build = false;

        // 2) HARD RESET of *known* global state for BTMiner
        // (Only touch what we know exists in btminer namespace)
        btminer::ClearCollected();   // clear collected patterns
        btminer::Tree.clear();       // clear MDD tree
        btminer::DFS.clear();        // clear DFS patterns

        btminer::M      = 0;
        btminer::L      = 0;
        btminer::N      = 0;
        btminer::theta  = 0;
        btminer::E      = 0;
        btminer::num_patt = 0;       // reset pattern counter if defined

        // NOTE: we do NOT reinsert root here; btminer::Load_instance()
        // is responsible for calling Tree.emplace_back(0,0,0) as needed.

        btminer::start_time = std::clock();

        // 3) Handle input (path or list-of-lists)
        TempFile tmp;
        std::string path;

        if (py::isinstance<py::str>(data)) {
            // File path: use directly
            path = data.cast<std::string>();
        } else {
            // Python list → write to a temp file in professor’s format
            auto seqs = data.cast<std::vector<std::vector<int>>>();
            tmp.path  = write_temp_seq_file(seqs);
            path      = tmp.path;
        }

        if (verbose) {
            std::cerr << "[BTMiner] path=" << path
                      << " minsup=" << minsup
                      << " preproc=" << preproc
                      << " use_dic=" << use_dic
                      << std::endl;
        }

        // 4) Build MDD + run miner
        if (!btminer::Load_instance(path, minsup)) {
            throw std::runtime_error("BTMiner: failed to load instance from: " + path);
        }

        btminer::Freq_miner();

        // 5) Return results
        py::dict out;
        out["patterns"]     = btminer::GetCollected();
        out["num_patterns"] = btminer::num_patt;
        out["time"]         = btminer::give_time(std::clock() - btminer::start_time);
        out["N"]            = btminer::N;
        out["L"]            = btminer::L;
        out["theta"]        = btminer::theta;
        return out;
    },
    py::arg("data"),
    py::arg("minsup")     = 0.01,
    py::arg("time_limit") = 36000,
    py::arg("preproc")    = false,
    py::arg("use_dic")    = false,
    py::arg("verbose")    = false,
    py::arg("out_file")   = ""
); */
m.def("BTMiner",
    [](py::object data,
       double minsup,
       unsigned int time_limit,
       bool preproc,
       bool use_dic,
       bool verbose,
       const std::string &out_file)
    {
        // 1) Configure globals
        btminer::time_limit = static_cast<int>(time_limit);
        btminer::pre_pro    = preproc;
        btminer::use_dic    = use_dic;
        btminer::b_disp     = verbose;
        btminer::b_write    = !out_file.empty();
        btminer::out_file   = out_file;
        btminer::N_mult     = 1;
        btminer::M_mult     = 1;
        btminer::just_build = false;

        // 2) CRASH FIX: Force Memory Release
        // .clear() keeps capacity reserved (leading to OOM). .swap() releases it.
        std::vector<btminer::Arc>().swap(btminer::Tree);
        std::vector<btminer::Pattern>().swap(btminer::DFS);
        std::vector<std::vector<int>>().swap(btminer::collectedPatterns);
        std::vector<int>().swap(btminer::freq);
        std::vector<int>().swap(btminer::item_dic);
        btminer::item_map.clear();
        btminer::item_map_rev.clear();

        // reset scalar globals
        btminer::M         = 0;
        btminer::L         = 0;
        btminer::N         = 0;
        btminer::theta     = 0;
        btminer::E         = 0;
        btminer::num_patt  = 0;
        btminer::num_nodes = 0;
        btminer::cur_node  = 0;

        btminer::start_time = std::clock();

        // 3) Handle input (path or list-of-lists)
        TempFile tmp;
        std::string path;

        if (py::isinstance<py::str>(data)) {
            path = data.cast<std::string>();
        } else {
            auto seqs = data.cast<std::vector<std::vector<int>>>();
            tmp.path  = write_temp_seq_file(seqs);
            path      = tmp.path;
        }

        if (verbose) {
            std::cerr << "[BTMiner] path=" << path
                      << " minsup=" << minsup
                      << " preproc=" << preproc
                      << " use_dic=" << use_dic
                      << std::endl;
        }

        // 4) Build MDD + run miner
        if (!btminer::Load_instance(path, minsup)) {
            throw std::runtime_error("BTMiner: failed to load instance from: " + path);
        }

        btminer::Freq_miner();

        // 5) Return results
        py::dict out;
        out["patterns"]     = btminer::GetCollected();
        out["num_patterns"] = btminer::num_patt;
        out["time"]         = btminer::give_time(std::clock() - btminer::start_time);
        out["N"]            = btminer::N;
        out["L"]            = btminer::L;
        out["theta"]        = btminer::theta;
        return out;
    },
    py::arg("data"),
    py::arg("minsup")     = 0.01,
    py::arg("time_limit") = 36000,
    py::arg("preproc")    = false,
    py::arg("use_dic")    = false,
    py::arg("verbose")    = false,
    py::arg("out_file")   = ""
); 
m.def("HTMiner",
    [](py::object data,
       double minsup,
       unsigned int time_limit,
       bool /*preproc*/,   // Python arg is ignored internally
       bool use_dic,
       bool verbose,
       const std::string &out_file)
    {
        using namespace htminer;

        // ───────── Global parameter setup ─────────
        htminer::time_limit = time_limit;

        // IMPORTANT: always run with preprocessing ON,
        // regardless of the Python `preproc` flag.
        htminer::pre_pro    = true;
        htminer::use_dic    = use_dic;
        htminer::just_build = false;
        htminer::b_disp     = verbose;
        htminer::b_write    = !out_file.empty();
        htminer::out_file   = out_file;

        // ───────── HARD RESET of HTMiner globals ─────────
        htminer::ClearCollected();
        htminer::Tree.clear();
        htminer::VTree.clear();
        htminer::CTree.clear();
        htminer::DFS.clear();
        htminer::VDFS.clear();
        htminer::item_dic.clear();

        htminer::M             = 0;
        htminer::N             = 0;
        htminer::L             = 0;
        htminer::E             = 0;
        htminer::theta         = 0;
        htminer::mlim          = 0;
        htminer::itmset_exists = false;

        // NOTE: do NOT add a root arc here;
        // htminer::Load_instance() already does Tree.emplace_back(0,0,0)
        htminer::start_time = std::clock();

        // ───────── Handle input (path or in-memory sequences) ─────────
        TempFile    tmp;
        std::string path;

        if (py::isinstance<py::str>(data)) {
            // data is a file path
            path = data.cast<std::string>();
        } else {
            // data is a list[list[int]] → write a temp file in the same text format
            auto seqs = data.cast<std::vector<std::vector<int>>>();
            tmp.path  = write_temp_seq_file(seqs);
            path      = tmp.path;
        }

        if (verbose) {
            std::cerr << "[HTMiner] path=" << path
                      << " minsup=" << minsup
                      << " preproc(always)=true"
                      << " use_dic=" << use_dic
                      << std::endl;
        }

        // ───────── Build MDD via professor's loader ─────────
        if (!htminer::Load_instance(path, minsup)) {
            throw std::runtime_error("HTMiner: failed to load instance from: " + path);
        }

        // ───────── Run miner ─────────
        htminer::Freq_miner();

        // ───────── Return results ─────────
        py::dict out;
        out["patterns"] = htminer::GetCollected();
        out["time"]     = htminer::give_time(std::clock() - htminer::start_time);
        return out;
    },
    py::arg("data"),
    py::arg("minsup")     = 0.01,
    py::arg("time_limit") = 36000,
    py::arg("preproc")    = false,   // kept for API symmetry, but IGNORED
    py::arg("use_dic")    = false,
    py::arg("verbose")    = false,
    py::arg("out_file")   = ""
);
m.def("LargePrefixProjection",
        [](py::object data,
           double minsup,
           unsigned int time_limit,
           bool preproc,
           bool use_dic,
           bool verbose,
           const std::string &out_file)
        {
            largepp::time_limit = time_limit;
            largepp::pre_pro    = preproc;
            largepp::use_dic    = use_dic;
            largepp::use_list   = true;
            largepp::b_disp     = verbose;
            largepp::b_write    = !out_file.empty();
            largepp::out_file   = out_file;
            largepp::just_build = false;

            // MEMORY FIX
            largepp::ClearCollected();
            std::vector<std::vector<int>>().swap(largepp::items);
            std::vector<int>().swap(largepp::item_dic);
            std::vector<largepp::Pattern>().swap(largepp::DFS);

            largepp::M = 0; largepp::L = 0; largepp::N = 0; 
            largepp::theta = 0; largepp::E = 0; largepp::num_patt = 0;
            largepp::start_time = std::clock();

            if (py::isinstance<py::str>(data)) {
                std::string fname = data.cast<std::string>();
                if (!largepp::Load_instance(fname, minsup)) {
                     throw std::runtime_error("LargePrefixProjection: failed to load");
                }
            } else {
                largepp::Load_py(data, minsup);
            }

            largepp::Freq_miner();

            py::dict out;
            out["patterns"] = largepp::GetCollected();
            out["time"]     = largepp::give_time(std::clock() - largepp::start_time);
            return out;
        },
        py::arg("data"), py::arg("minsup") = 0.01, py::arg("time_limit") = 36000,
        py::arg("preproc") = false, py::arg("use_dic") = false,
        py::arg("verbose") = false, py::arg("out_file") = ""
    );


    // ─────────────────────────────────────────────────────────────
    // LargeBTMiner  (always uses professor's largebm::Load_instance)
    // ─────────────────────────────────────────────────────────────
    m.def("LargeBTMiner",
        [](py::object data,
           double minsup,
           unsigned int time_limit,
           bool preproc,
           bool use_dic,
           bool verbose,
           const std::string &out_file)
        {
            using namespace largebm;

            largebm::time_limit = time_limit;
            largebm::pre_pro    = preproc;
            largebm::use_dic    = use_dic;
            largebm::use_list   = false;   // MDD-based
            largebm::b_disp     = verbose;
            largebm::b_write    = !out_file.empty();
            largebm::out_file   = out_file;
            largebm::just_build = false;

            largebm::ClearCollected();
            largebm::items.clear();
            largebm::item_dic.clear();
            largebm::inv_item_dic.clear();
            largebm::Tree.clear();
            largebm::DFS.clear();

            largebm::start_time = std::clock();

            TempFile tmp;
            std::string path;

            if (py::isinstance<py::str>(data)) {
                path = data.cast<std::string>();
            } else {
                auto seqs = data.cast<std::vector<std::vector<int>>>();
                tmp.path  = write_temp_seq_file(seqs);
                path      = tmp.path;
            }

            if (verbose) {
                std::cerr << "[LargeBTMiner] path=" << path
                          << " minsup=" << minsup
                          << " preproc=" << preproc
                          << " use_dic=" << use_dic
                          << std::endl;
            }

            if (!largebm::Load_instance(path, minsup)) {
                throw std::runtime_error("LargeBTMiner: failed to load instance from: " + path);
            }

            largebm::Freq_miner();

            py::dict out;
            out["patterns"] = largebm::GetCollected();
            out["time"]     = largebm::give_time(std::clock() - largebm::start_time);
            return out;
        },
        py::arg("data"),
        py::arg("minsup")     = 0.01,
        py::arg("time_limit") = 36000,
        py::arg("preproc")    = false,
        py::arg("use_dic")    = false,
        py::arg("verbose")    = false,
        py::arg("out_file")   = ""
    );

    // ─────────────────────────────────────────────────────────────
    // LargeHTMiner  (always uses professor's largehm::Load_instance; pre_pro forced ON)
    // ─────────────────────────────────────────────────────────────
      // ─────────────────────────────────────────────────────────────
    // LargeHTMiner  (professor's Large HTMiner, namespaced as largehm)
    // ─────────────────────────────────────────────────────────────
    m.def("LargeHTMiner",
        [](py::object data,
           double minsup,
           unsigned int time_limit,
           bool /*preproc*/,   // kept for API symmetry; ignored
           bool use_dic,
           bool verbose,
           const std::string &out_file)
        {
            using namespace largehm;

            // 1) Global configuration (mirror professor's style)
            largehm::time_limit = time_limit;
            largehm::pre_pro    = true;       // always preprocess
            largehm::use_dic    = use_dic;
            largehm::just_build = false;
            largehm::b_disp     = verbose;
            largehm::b_write    = !out_file.empty();
            largehm::out_file   = out_file;

            // 2) HARD RESET of all global state for a fresh run
            largehm::ClearCollected();   // our helper in largehm::utility.cpp

            largehm::M             = 0;
            largehm::L             = 0;
            largehm::mlim          = 0;
            largehm::N             = 0;
            largehm::theta         = 0;
            largehm::E             = 0;
            largehm::itmset_exists = false;

            // containers
            // (item_dic reset is optional and not strictly needed here)
            largehm::DFS.clear();
            largehm::VDFS.clear();
            largehm::Tree.clear();
            largehm::VTree.clear();
            largehm::CTree.clear();

            largehm::start_time = std::clock();

            // 3) Handle input (file path or Python list)
            TempFile tmp;
            std::string path;

            if (py::isinstance<py::str>(data)) {
                path = data.cast<std::string>();
            } else {
                auto seqs = data.cast<std::vector<std::vector<int>>>();
                tmp.path  = write_temp_seq_file(seqs);
                path      = tmp.path;
            }

            if (verbose) {
                std::cerr << "[LargeHTMiner] path=" << path
                          << " minsup=" << minsup
                          << " preproc(always)=true"
                          << " use_dic=" << use_dic
                          << std::endl;
            }

            // 4) Build MDD / load instance.
            //    NOTE: Load_instance() itself does Tree.emplace_back(0,0,0),
            //    so we DO NOT create a root node here.
            if (!largehm::Load_instance(path, minsup)) {
                throw std::runtime_error("LargeHTMiner: failed to load instance from: " + path);
            }

            // 5) Run miner (same timing logic as original main)
            if (!largehm::just_build &&
                largehm::give_time(std::clock() - largehm::start_time) < largehm::time_limit)
            {
                largehm::Freq_miner();
                if (largehm::give_time(std::clock() - largehm::start_time) >= largehm::time_limit) {
                    std::cout << "TIME LIMIT REACHED\n";
                }
            }

            // 6) Return collected patterns + runtime
            py::dict out;
            out["patterns"] = largehm::GetCollected();
            out["time"]     = largehm::give_time(std::clock() - largehm::start_time);
            return out;
        },
        py::arg("data"),
        py::arg("minsup")     = 0.01,
        py::arg("time_limit") = 36000,
        py::arg("preproc")    = false,   // kept for API symmetry
        py::arg("use_dic")    = false,
        py::arg("verbose")    = false,
        py::arg("out_file")   = ""
    );


}


/*#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
#include <iostream>


// PrefixProjection headers
#include "freq_miner.hpp"
#include "load_inst.hpp"
#include "utility.hpp"

// BTMiner (wrapped in its own namespace in source files)
#include "btminer/src/freq_miner.hpp"
#include "btminer/src/load_inst.hpp"
#include "btminer/src/utility.hpp"
#include "btminer/src/build_mdd.hpp"

// HTMiner (wrapped in its own namespace in source files)
#include "htminer/src/build_mdd.hpp"    // ← ensure HTMiner MDD builder is available
#include "htminer/src/freq_miner.hpp"
#include "htminer/src/load_inst.hpp"
#include "htminer/src/utility.hpp"


#include "largepp/src/freq_miner.hpp"
#include "largepp/src/load_inst.hpp"
#include "largepp/src/utility.hpp"


#include "largebm/src/freq_miner.hpp"
#include "largebm/src/load_inst.hpp"
#include "largebm/src/utility.hpp"
#include "largebm/src/build_mdd.hpp"

#include "largehm/src/freq_miner.hpp"
#include "largehm/src/load_inst.hpp"
#include "largehm/src/utility.hpp"
#include "largehm/src/build_mdd.hpp"



PYBIND11_MODULE(_effspm, m) {
    m.doc() = "Unified SPM library: PrefixProjection, BTMiner, HTMiner";

    // ─────────────────────────────────────────────────────────────
    // PrefixProjection
    // ─────────────────────────────────────────────────────────────
    m.def("PrefixProjection",
        [](py::object data,
           double minsup,
           unsigned int time_limit,
           bool preproc,
           bool use_dic,
           bool verbose,
           const std::string &out_file)
        {
            ::time_limit = time_limit;
            ::pre_pro    = preproc;
            ::use_dic    = use_dic;
            ::use_list   = false;
            ::b_disp     = verbose;
            ::b_write    = !out_file.empty();
            ::out_file   = out_file;

            ClearCollected();
            start_time = std::clock();

            if (py::isinstance<py::str>(data)) {
                std::string path = data.cast<std::string>();
                if (!Load_instance(path, minsup))
                    throw std::runtime_error("Failed to load file: " + path);
            } else {
                auto seqs = data.cast<std::vector<std::vector<int>>>();
                items = std::move(seqs);
                N = items.size();

                int max_id = 0;
                for (auto &seq : items)
                    for (int x : seq)
                        max_id = std::max(max_id, std::abs(x));
                L = max_id;

                theta = (minsup < 1.0) ? std::ceil(minsup * N) : minsup;

                DFS.clear();
                DFS.reserve(L);
                for (unsigned int i = 0; i < L; ++i)
                    DFS.emplace_back(-static_cast<int>(i) - 1);

                M = 0;
                E = 0;
                for (auto &seq : items) {
                    M = std::max<unsigned int>(M, seq.size());
                    E += seq.size();
                }
            }

            Freq_miner();

            py::dict out;
            out["patterns"] = GetCollected();
            out["time"]     = give_time(std::clock() - start_time);
            return out;
        },
        py::arg("data"),
        py::arg("minsup")     = 0.01,
        py::arg("time_limit") = 36000,
        py::arg("preproc")    = false,
        py::arg("use_dic")    = false,
        py::arg("verbose")    = false,
        py::arg("out_file")   = ""
    );
m.def("BTMiner",
    [](py::object data,
       double minsup,
       unsigned int time_limit,
       bool preproc,
       bool use_dic,
       bool verbose,
       const std::string &out_file)
    {
        // We are calling the *professor* BTMiner, now namespaced as btminer::.
        // So we only set the globals the professor code actually has.

        // 1) configure professor globals
        btminer::time_limit = static_cast<int>(time_limit);
        btminer::pre_pro    = preproc;
        btminer::use_dic    = use_dic;
        btminer::b_disp     = verbose;
        btminer::b_write    = !out_file.empty();
        btminer::out_file   = out_file;
        btminer::N_mult     = 1;         // professor uses these too
        btminer::M_mult     = 1;
        btminer::just_build = false;     // we want full mining

        btminer::start_time = std::clock();

        // 2) load data
        //
        // Professor’s code is primarily file-based (Load_instance(const string&, double)).
        // So: if user passes a file path → use the professor loader directly.
        // If user passes a Python list-of-lists → we will build the MDD the same
        // way professor’s loader does, but without changing his logic.
        if (py::isinstance<py::str>(data)) {
            // ----- FILE MODE -----
            std::string path = data.cast<std::string>();

            if (verbose) {
                std::cerr << "[BT][binding] file=" << path
                          << " minsup=" << minsup
                          << " preproc=" << preproc << std::endl;
            }

            if (!btminer::Load_instance(path, minsup)) {
                throw std::runtime_error("BTMiner: failed to load file: " + path);
            }
        } else {
            // ----- PYTHON LIST MODE -----
            //
            // We mimic professor’s loader:
            //   - create root in Tree
            //   - compute N, M, L
            //   - compute theta from minsup
            //   - seed DFS (one Pattern per item, as in Preprocess branch)
            //   - call Build_MDD(...) for each sequence
            //
            // This DOES NOT change his mining logic; it just drives it from memory.

            auto seqs = data.cast<std::vector<std::vector<int>>>();

            // clear MDD and globals to a known state
            btminer::Tree.clear();
            btminer::Tree.emplace_back(0, 0, 0);      // root (exactly like professor)

            // compute basic stats
            int max_id = 0;
            int max_len = 0;
            int seq_count = 0;
            long long entries = 0;

            for (const auto &s : seqs) {
                if (s.empty()) continue;
                ++seq_count;
                max_len = std::max<int>(max_len, static_cast<int>(s.size()));
                for (int x : s) {
                    max_id = std::max(max_id, std::abs(x));
                    ++entries;
                }
            }

            btminer::N = seq_count;
            btminer::M = max_len;
            btminer::L = max_id;
            btminer::E = static_cast<int>(entries);

            // theta = abs support
            if (minsup < 1.0)
                btminer::theta = static_cast<int>(std::ceil(minsup * btminer::N * btminer::N_mult));
            else
                btminer::theta = static_cast<int>(minsup);

            // seed DFS exactly like professor does in the preprocessed branch:
            btminer::DFS.clear();
            btminer::DFS.reserve(btminer::L);
            for (int i = 0; i < btminer::L; ++i)
                btminer::DFS.emplace_back(-i - 1);

            // now build the MDD, sequence by sequence
            for (const auto &s : seqs) {
                if (s.empty()) continue;
                // professor’s Build_MDD takes a vector<int> by non-const ref
                std::vector<int> tmp = s;
                btminer::Build_MDD(tmp);
            }

            if (verbose) {
                std::cerr << "[BT][binding] PY mode: N=" << btminer::N
                          << " L=" << btminer::L
                          << " M=" << btminer::M
                          << " E=" << btminer::E
                          << " theta=" << btminer::theta
                          << " Tree.size()=" << btminer::Tree.size()
                          << std::endl;
            }
        }

        // 3) run professor’s miner
        btminer::Freq_miner();

        // 4) build python result
        // 4) build python result
            py::dict out;
            out["patterns"]     = btminer::GetCollected();    // ← NEW
            out["num_patterns"] = btminer::num_patt;
            out["time"]         = btminer::give_time(std::clock() - btminer::start_time);
            out["N"]            = btminer::N;
            out["L"]            = btminer::L;
            out["theta"]        = btminer::theta;
            return out;

    },
    py::arg("data"),
    py::arg("minsup")     = 0.01,
    py::arg("time_limit") = 36000,
    py::arg("preproc")    = false,
    py::arg("use_dic")    = false,
    py::arg("verbose")    = false,
    py::arg("out_file")   = ""
);



       
// HTMiner
// ─────────────────────────────────────────────────────────────
// HTMiner
m.def("HTMiner",
    [](py::object data,
       double minsup, unsigned int time_limit,
       bool preproc, bool use_dic,
       bool verbose, const std::string &out_file)
    {
        htminer::time_limit = time_limit;
        htminer::pre_pro    = preproc;
        htminer::use_dic    = use_dic;
        htminer::just_build = false;
        htminer::use_list   = false;
        htminer::b_disp     = verbose;
        htminer::b_write    = !out_file.empty();
        htminer::out_file   = out_file;
        htminer::ClearCollected();
        htminer::start_time = std::clock();

        if (py::isinstance<py::str>(data)) {
            std::string path = data.cast<std::string>();
            if (!htminer::Load_instance(path, minsup))
                throw std::runtime_error("Failed to load file: " + path);
        } else {
            auto seqs = data.cast<std::vector<std::vector<int>>>();
            htminer::items = std::move(seqs);
            htminer::N = htminer::items.size();

            int max_id = 0;
            htminer::M = 0;
            htminer::E = 0;
            for (auto &seq : htminer::items) {
                htminer::M = std::max<unsigned int>(htminer::M, seq.size());
                for (int x : seq)
                    max_id = std::max(max_id, std::abs(x));
                htminer::E += seq.size();
            }
            htminer::L = max_id;
            htminer::theta = (minsup < 1.0)
                              ? static_cast<unsigned long long>(std::ceil(minsup * htminer::N))
                              : static_cast<unsigned long long>(minsup);

            htminer::DFS.clear();
            htminer::DFS.reserve(htminer::L);
            for (unsigned int i = 0; i < static_cast<unsigned int>(htminer::L); ++i)
                htminer::DFS.emplace_back(-static_cast<int>(i) - 1);

            htminer::VDFS.clear();
            htminer::VDFS.resize(htminer::L);
        }

        htminer::Freq_miner();

        // 👇 now really respects verbose
        if (verbose) {
            std::cout << " total patterns = "
                      << htminer::collectedPatterns.size() << "\n";
        }

        py::dict out;
        out["patterns"] = htminer::GetCollected();
        out["time"]     = htminer::give_time(std::clock() - htminer::start_time);
        return out;
    },
    py::arg("data"),
    py::arg("minsup")     = 0.01,
    py::arg("time_limit") = 36000,
    py::arg("preproc")    = false,
    py::arg("use_dic")    = false,
    py::arg("verbose")    = false,
    py::arg("out_file")   = ""
);


 m.def("LargePrefixProjection",
    [](py::object data,
       double minsup,
       unsigned int time_limit,
       bool preproc,
       bool use_dic,
       bool verbose,
       const std::string &out_file)
    {
        largepp::time_limit = time_limit;
        largepp::pre_pro    = preproc;
        largepp::use_dic    = use_dic;
        largepp::use_list   = true;
        largepp::b_disp     = verbose;
        largepp::b_write    = !out_file.empty();
        largepp::out_file   = out_file;
        largepp::just_build = false;

        largepp::ClearCollected();
        largepp::start_time = std::clock();

        // 👇 this was the last noisy one
        if (verbose) {
            std::cerr << " minsup=" << minsup
                      << " preproc=" << preproc
                      << " verbose=" << verbose
                      << " out_file=" << (out_file.empty() ? "(none)" : out_file)
                      << " use_dic=" << use_dic << "\n";
        }

        if (py::isinstance<py::str>(data)) {
            std::string fname = data.cast<std::string>();
            largepp::Load_instance(fname, minsup);
        } else {
            largepp::Load_py(data, minsup);
        }

        largepp::Freq_miner();

        py::dict out;
        out["patterns"] = largepp::GetCollected();
        out["time"]     = largepp::give_time(std::clock() - largepp::start_time);
        return out;
    },
    py::arg("data"),
    py::arg("minsup")     = 0.01,
    py::arg("time_limit") = 36000,
    py::arg("preproc")    = false,
    py::arg("use_dic")    = false,
    py::arg("verbose")    = false,
    py::arg("out_file")   = ""
);



// ─────────────────────────────────────────────────────────────
// LargeBTMiner  -- Python wrapper for the largebm implementation
// ─────────────────────────────────────────────────────────────
// m.def(
//     "LargeBTMiner",
//     [](py::object        data,
//        double            minsup ,
//        unsigned int      time_limit,
//        bool              preproc ,
//        bool              use_dic,
//        bool              verbose,
//        const std::string &out_file ) 
//     {
//         /* 1) Global flags */
//         largebm::time_limit = time_limit;
//         largebm::pre_pro    = preproc;
//         largebm::use_dic    = use_dic;
//         largebm::use_list   = false;          // large-mode → always MDD
//         largebm::just_build = false;
//         largebm::b_disp     = verbose;
//         largebm::b_write    = !out_file.empty();
//         largebm::out_file   = out_file;

//         /* 2) Reset per-run state */
//         largebm::ClearCollected();
//         largebm::start_time = std::clock();

//         /* 3) Load the DB (file path or in-memory list<list<int>>) */
//         if (py::isinstance<py::str>(data)) {
//             std::string path = data.cast<std::string>();
//             if (!largebm::Load_instance(path, minsup))
//                 throw std::runtime_error("Failed to load file: " + path);
//         } else {
//             // In-memory sequences
//             largebm::items = std::move(data.cast<std::vector<std::vector<int>>>());
//             largebm::N     = static_cast<unsigned int>(largebm::items.size());

//             /* -- basic stats -- */
//             int max_id = 0;
//             largebm::M = 0;
//             largebm::E = 0;
//             for ( auto &seq : largebm::items) {
//                 largebm::M = std::max<unsigned int>(largebm::M,
//                                                     static_cast<unsigned int>(seq.size()));
//                 largebm::E += static_cast<unsigned long long>(seq.size());
//                 for (int x : seq) max_id = std::max(max_id, std::abs(x));
//             }
//             largebm::L     = static_cast<unsigned int>(max_id);
//             largebm::theta = (minsup < 1.0)
//                                ? static_cast<unsigned long long>(std::ceil(minsup * largebm::N))
//                                : static_cast<unsigned long long>(minsup);

//             /* -- DFS buffer (size = L) -- */
//             largebm::DFS.clear();
//             largebm::DFS.reserve(largebm::L);
//             for (unsigned int i = 0; i < largebm::L; ++i)
//                 largebm::DFS.emplace_back(-static_cast<int>(i) - 1);

//             /* -- Build the MDD -- */
//             largebm::Tree.clear();
//             largebm::Tree.emplace_back(0, 0, 0);              // dummy root
//             for ( auto &seq : largebm::items)
//                 largebm::Build_MDD(seq);
//         }

//         /* 4) Mine and return results */
//         largebm::Freq_miner();

//         py::dict out;
//         out["patterns"] = largebm::GetCollected();
//         out["time"]     = largebm::give_time(std::clock() - largebm::start_time);
//         return out;
//     },
//     py::arg("data"),
//     py::arg("minsup")     = 0.01,
//     py::arg("time_limit") = 36000,
//     py::arg("preproc")    = false,
//     py::arg("use_dic")    = false,
//     py::arg("verbose")    = false,
//     py::arg("out_file")   = ""
// );



    // ─────────────────────────────────────────────────────────────────────────
    // LargeBTMiner (MDD-based)
    // ─────────────────────────────────────────────────────────────────────────
 /*m.def("LargeBTMiner",
    [](py::object data,
       double minsup,
       unsigned int time_limit,
       bool preproc,
       bool use_dic,
       bool verbose,
       const std::string &out_file)
    {
        using namespace largebm;

        // 0) Set global flags and timers
        largebm::time_limit = time_limit;
        largebm::pre_pro    = preproc;
        largebm::use_dic    = use_dic;
        largebm::use_list   = false;         // large-mode → always MDD
        largebm::b_disp     = verbose;
        largebm::b_write    = !out_file.empty();
        largebm::out_file   = out_file;
        largebm::just_build = false;

        // 0.1) Clear any leftover data/state from previous runs
        largebm::items.clear();
        largebm::item_dic.clear();
        largebm::inv_item_dic.clear();
        largebm::Tree.clear();
        largebm::DFS.clear();
        largebm::ClearCollected();

        // 1) Load sequences (either from filename or from Python list)
        if (py::isinstance<py::str>(data)) {
            // ─────────── FILE-BASED MODE ───────────
            std::string path = data.cast<std::string>();
            if (!largebm::Load_instance(path, minsup))
                throw std::runtime_error("Failed to load file: " + path);

        } else {
            // ────────── IN-MEMORY MODE ──────────
            auto seqs = data.cast<std::vector<std::vector<int>>>();
            largebm::items = std::move(seqs);
            largebm::N     = largebm::items.size();

            // 1.1) Compute basic DB statistics (M, E, L) and absolute support θ
            int max_id = 0;
            largebm::M = 0;
            largebm::E = 0;
            for (auto &seq : largebm::items) {
                largebm::M = std::max<unsigned int>(largebm::M, static_cast<unsigned int>(seq.size()));
                largebm::E += static_cast<unsigned long long>(seq.size());
                for (int x : seq) max_id = std::max(max_id, std::abs(x));
            }
            largebm::L = static_cast<unsigned int>(max_id);
            largebm::theta = (minsup < 1.0)
                               ? static_cast<unsigned long long>(std::ceil(minsup * largebm::N))
                               : static_cast<unsigned long long>(minsup);

            // 1.2) Initialize DFS buffer (size = L)
            largebm::DFS.reserve(largebm::L);
            for (unsigned int i = 0; i < largebm::L; ++i)
                largebm::DFS.emplace_back(-static_cast<int>(i) - 1);

            // 1.3) Build the MDD “Tree”
            // Insert one dummy root node (item=0, freq=0, anct=0)
            largebm::Tree.emplace_back(0, 0, 0);
            for (auto &seq : largebm::items)
                largebm::Build_MDD(const_cast<std::vector<int>&>(seq));
        }

        // 2) Rebuild inverse-dictionary from fresh item_dic
        {
            std::vector<int> inv(largebm::item_dic.size() + 1);
            for (int old = 1; old <= static_cast<int>(largebm::item_dic.size()); ++old) {
                int cid = largebm::item_dic[old - 1];
                if (cid > 0) inv[cid] = old;
            }
            largebm::inv_item_dic = std::move(inv);
        }

        // 3) Start timing and run the miner
        largebm::start_time = std::clock();
        largebm::Freq_miner();

        // 4) Collect results and elapsed time
        const auto& pats = largebm::GetCollected();

        py::dict out;
        out["patterns"] = pats;
        out["time"]     = largebm::give_time(std::clock() - largebm::start_time);
        return out;
    },
    py::arg("data"),
    py::arg("minsup")     = 0.01,
    py::arg("time_limit") = 36000,
    py::arg("preproc")    = false,
    py::arg("use_dic")    = false,
    py::arg("verbose")    = false,
    py::arg("out_file")   = ""
);


   
m.def("LargeHTMiner",
    [](py::object data,
       double minsup,
       unsigned int time_limit,
       bool preproc,
       bool use_dic,
       bool verbose,
       const std::string &out_file)
    {
        // 0) Set global flags and timers:
        largehm::time_limit = time_limit;
        largehm::pre_pro    = preproc;
        largehm::use_dic    = use_dic;
        largehm::use_list   = true;     // force in‐memory mode
        largehm::b_disp     = verbose;
        largehm::b_write    = !out_file.empty();
        largehm::out_file   = out_file;
        largehm::just_build = false;

        largehm::ClearCollected();
        largehm::start_time = std::clock();

        if (py::isinstance<py::str>(data)) {
            // ───────────── FILE‐BASED MODE ─────────────
            // Force mlim so that every item lands in temp_vec (never temp_lim):
            largehm::mlim = UINT_MAX;

            std::string path = data.cast<std::string>();
            if (! largehm::Load_instance(path, minsup))
                throw std::runtime_error("Failed to load file: " + path);
        }
        else {
            // ───────────── IN‐MEMORY MODE ─────────────
            auto seqs = data.cast<std::vector<std::vector<int>>>();
            largehm::items = std::move(seqs);
            largehm::N     = largehm::items.size();

            // 1) Compute L = maximum absolute item ID
            int max_id = 0;
            for (auto &seq : largehm::items)
                for (int x : seq)
                    max_id = std::max(max_id, std::abs(x));
            largehm::L = static_cast<unsigned int>(max_id);

            // 2) Compute theta as absolute support threshold
            largehm::theta = (minsup < 1.0)
                               ? static_cast<unsigned long long>(std::ceil(minsup * largehm::N))
                               : static_cast<unsigned long long>(minsup);

            // 3) Initialize DFS (size = L)
            largehm::DFS.clear();
            largehm::DFS.reserve(largehm::L);
            for (unsigned int i = 0; i < largehm::L; ++i)
                largehm::DFS.emplace_back(-static_cast<int>(i) - 1);

            // 4) Compute M (max sequence length) and E (total entries)
            largehm::M = 0;
            largehm::E = 0;
            for (auto &seq : largehm::items) {
                largehm::M = std::max<unsigned int>(
                    largehm::M, static_cast<unsigned int>(seq.size()));
                largehm::E += seq.size();
            }

            // 5) ─── Build the MDD “manually” ───
            largehm::Tree.clear();
            largehm::VTree.clear();
            largehm::CTree.clear();

            // Insert exactly one dummy root node (chld=0, sibl=0, freq=0):
            largehm::Tree.emplace_back(0,0,0);

            // For each sequence “seq”, insert into MDD by placing a single −1 sentinel:
            for (auto &seq : largehm::items) {
                // Copy the item IDs:
                std::vector<int> temp_vec = seq;
                // Only a single “−1” is needed to force the suffix insertion:
                std::vector<int> temp_lim(1, -1);

                largehm::Build_MDD(temp_vec, temp_lim);
            }

            
        }

        // 6) Run the frequency miner (Tree is now properly built):
        largehm::Freq_miner();

        // 7) Return results to Python:
        py::dict out;
        out["patterns"] = largehm::GetCollected();
        out["time"]     = largehm::give_time(std::clock() - largehm::start_time);
        return out;
    },
    py::arg("data"),
    py::arg("minsup")     = 0.01,
    py::arg("time_limit") = 36000,
    py::arg("preproc")    = false,
    py::arg("use_dic")    = false,
    py::arg("verbose")    = false,
    py::arg("out_file")   = ""
);




} */