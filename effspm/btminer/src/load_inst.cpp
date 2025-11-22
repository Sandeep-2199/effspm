// effspm/btminer/src/load_inst.cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <time.h>

#include "load_inst.hpp"
#include "utility.hpp"
#include "build_mdd.hpp"
#include "freq_miner.hpp"

namespace btminer {

using namespace std;

// ---------------------------------------------------------------------
// global definitions
// ---------------------------------------------------------------------
int  M = 0;
int  N = 0;
int  L = 0;
unsigned long long E = 0ULL;
int  num_nodes = 0;
int  theta = 0;
int  cur_node = 0;

map<string, int>  item_map;
map<int, string>  item_map_rev;

std::vector<int>  freq;
std::vector<int>  item_dic;
std::vector<std::vector<int>> items; 
std::vector<Pattern> DFS;

string out_file, folder;
bool   b_disp    = 0;
bool   b_write   = 0;
bool   use_dic   = 0;
bool   just_build= 0;
bool   pre_pro   = 1;

int    N_mult    = 1;
int    M_mult    = 1;
int    time_limit= 30 * 3600;

clock_t start_time;

void Load_items_pre(string &inst_name);
bool Load_items(string &inst_name);
bool Preprocess(string &inst, double thresh);

// ---------------------------------------------------------------------
// main loader
// ---------------------------------------------------------------------
bool Load_instance(string &items_file, double thresh) {
    clock_t kk = clock();

    Tree.emplace_back(0, 0, 0);

    if (pre_pro) {
        if (!Preprocess(items_file, thresh))
            return false;
        if (b_disp)
           cout << "\nPreprocess done in " << give_time(clock() - kk) << " seconds\n\n";

        DFS.clear();
        DFS.reserve(L);
        for (int i = 0; i < L; ++i)
            DFS.emplace_back(-i - 1);

        kk = clock();
        Load_items_pre(items_file);
    }
    else if (!Load_items(items_file)) {
        return false;
    }
    else {
        if (thresh < 1)
            theta = static_cast<int>(ceil(thresh * N * N_mult));
        else
            theta = static_cast<int>(thresh);
    }
    if (b_disp)
         cout << "\nMDD Database built in " << give_time(clock() - kk) << " seconds\n\n";
    if (b_disp)
        cout << "Found " << N * N_mult
             << " sequence, with max line len " << M
             << ", and " << L << " items, and " << E << " enteries\n";

    return true;
}

// ---------------------------------------------------------------------
// Preprocess (UNCHANGED as per your request)
// ---------------------------------------------------------------------
bool Preprocess(string &inst, double thresh) {
    N = 0;
    L = 0;
    freq.clear();
    item_dic.clear();
    item_map.clear();
    item_map_rev.clear();

    ifstream file(inst);

    if (file.good()) {
        string line;
        while (getline(file, line) && give_time(clock() - start_time) < time_limit) {
            ++N;
            vector<bool> counted(L, false);

            istringstream word(line);
            string itm;
            while (word >> itm) {
                int ditem = stoi(itm);
                if (L < abs(ditem))
                    L = abs(ditem);

                while (static_cast<int>(freq.size()) < L) {
                    freq.push_back(0);
                    counted.push_back(false);
                }

                int idx = abs(ditem) - 1;
                if (!counted[idx]) {
                    ++freq[idx];
                    counted[idx] = true;
                }
            }
        }
    } else {
        cout << "!!!!!! No such file exists: " << inst << " !!!!!!\n";
        return false;
    }

    if (thresh < 1)
        theta = static_cast<int>(ceil(thresh * N * N_mult));
    else
        theta = static_cast<int>(thresh);

    int real_L = 0;
    item_dic = vector<int>(L, -1);
    for (int i = 0; i < L; ++i) {
        if (freq[i] >= theta)
            item_dic[i] = ++real_L;
    }
    if (b_disp)
      cout << "Original number of items: " << L
         << " Reduced to: " << real_L << endl;

    L = real_L;
    N = 0; 

    return true;
}

// ---------------------------------------------------------------------
// Load_items_pre (UNCHANGED)
// ---------------------------------------------------------------------
void Load_items_pre(string &inst_name) {
    ifstream file(inst_name);

    if (file.good()) {
        string line;
        while (getline(file, line) && give_time(clock() - start_time) < time_limit) {
            istringstream word(line);
            string itm;
            vector<int> temp_vec;
            bool sgn = false;
            while (word >> itm) {
                int ditem;
                if (use_dic) {
                    auto it = item_map.find(itm);
                    if (it == item_map.end()) {
                        item_map[itm] = ++L;
                        item_map_rev[L] = itm;
                        ditem = L;
                    } else {
                        ditem = it->second;
                    }
                } else {
                    ditem = stoi(itm);
                }

                if (freq[abs(ditem) - 1] < theta) {
                    if (!sgn)
                        sgn = (ditem < 0);
                    continue;
                } else {
                    if (ditem > 0)
                        ditem = item_dic[ditem - 1];
                    else
                        ditem = -item_dic[-ditem - 1];
                }

                if (sgn) {
                    if (ditem > 0)
                        ditem = -ditem;
                    sgn = false;
                }

                temp_vec.push_back(ditem);
            }

            if (temp_vec.empty())
                continue;

            ++N;

            if (static_cast<int>(temp_vec.size()) > M)
                M = static_cast<int>(temp_vec.size());

            Build_MDD(temp_vec);
        }
    }
}

// ---------------------------------------------------------------------
// Load_items (FIXED: Moved DFS logic out of 'else')
// ---------------------------------------------------------------------
bool Load_items(string &inst_name) {
    ifstream file(inst_name);

    if (file.good()) {
        string line;
        while (getline(file, line) && give_time(clock() - start_time) < time_limit) {
            ++N;
            istringstream word(line);
            string itm;
            vector<int> temp_vec;
            while (word >> itm) {
                int ditem;
                if (use_dic) {
                    auto it = item_map.find(itm);
                    if (it == item_map.end()) {
                        item_map[itm] = ++L;
                        item_map_rev[L] = itm;
                        ditem = L;
                    } else {
                        ditem = it->second;
                    }
                    // CRASH FIX: Ensure L tracks max item even in dictionary mode
                    if (L < abs(ditem)) L = abs(ditem);
                } else {
                    ditem = stoi(itm);
                    if (L < abs(ditem)) {
                        L = abs(ditem);
                    }
                }

                // CRASH FIX: This loop used to be inside the 'else' block.
                // It must run regardless of use_dic mode, or DFS vector stays empty -> Crash.
                while (static_cast<int>(DFS.size()) < L && !just_build) {
                    DFS.reserve(L);
                    DFS.emplace_back(-((int)DFS.size()) - 1);
                }

                temp_vec.push_back(ditem);
            }

            if (static_cast<int>(temp_vec.size()) > M)
                M = static_cast<int>(temp_vec.size());

            Build_MDD(temp_vec);
        }
    } else {
        cout << "!!!!!! No such file exists: " << inst_name << " !!!!!!\n";
        return false;
    }

    return true;
}

} // namespace btminer