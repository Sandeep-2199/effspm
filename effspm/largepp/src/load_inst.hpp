#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <pybind11/pybind11.h>

#include "pattern.hpp"

namespace largepp {
using namespace std;

// public entry points
bool Load_instance(std::string& items_file, double thresh);
void Load_py(const pybind11::object& py_data, double thresh);

// shared state
extern std::vector<std::vector<int>> items;
extern std::string                    out_file;

// FIX: Expose item_dic so _effspm.cpp can clean it up
extern std::vector<int>               item_dic; 

extern bool b_disp, b_write, use_dic, just_build, ovr_count, pre_pro;
extern bool use_list;

extern unsigned int        M, L, time_limit;
extern unsigned long long  N;
extern double              theta;
extern unsigned long long  E;
extern std::clock_t        start_time;

// DFS queue of potential patterns to extend
extern std::vector<largepp::Pattern> DFS;

} // namespace largepp