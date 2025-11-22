#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include "load_inst.hpp"
#include "freq_miner.hpp"
#include "utility.hpp"

namespace largepp {   
using namespace std;

/* ------------------------------------------------------------------
 * Global definitions
 * ---------------------------------------------------------------- */
unsigned int        M = 0, L = 0;
unsigned long long  N = 0, E = 0;
double              theta = 0.01;
vector<vector<int>> items;
vector<Pattern>     DFS;
vector<int>         item_dic;

/* Forward decls  */
static bool  Load_items(string& inst);
static void  Load_items_pre(string& inst);
static bool  Preprocess(string& inst, double thresh);

/* ==================================================================
 * MAIN ENTRY
 * ================================================================= */
bool Load_instance(string& items_file, double thresh)
{
    clock_t kk = clock();

    if (pre_pro) {
        if (!Preprocess(items_file, thresh)) return false;

        if (b_disp) cout << "\nPreprocess done in " << give_time(clock() - kk) << " seconds\n\n";

        DFS.clear();
        DFS.reserve(L);
        for (unsigned int i = 0; i < L; ++i)
            DFS.emplace_back(-int(i) - 1);

        kk = clock();
        Load_items_pre(items_file);
        N = items.size();
    }
    else if (!Load_items(items_file))
        return false;
    else
        theta = (thresh < 1.0) ? ceil(thresh * N) : thresh;
    
    if (b_disp) {
        cout << "\nMDD Database built in " << give_time(clock() - kk) << " seconds\n";
        cout << "Found " << N << " sequence, with max line len " << M
             << ", and " << L << " items, and " << E << " enteries\n";
    }

    return true;
}

/* ==================================================================
 * ALT ENTRY — load directly from a Python list of lists
 * ================================================================= */
void Load_py(const pybind11::object& data, double thresh)
{
    items = data.cast<vector<vector<int>>>();
    N     = items.size();

    int max_id = 0;
    M = 0;  E = 0;
    for (auto& seq : items) {
        M = max<unsigned int>(M, static_cast<unsigned int>(seq.size()));
        E += seq.size();
        for (int x : seq)
            max_id = max(max_id, abs(x));
    }
    L = static_cast<unsigned int>(max_id);
    theta = (thresh < 1.0) ? ceil(thresh * N) : thresh;

    DFS.clear();
    DFS.reserve(L);
    for (unsigned int i = 0; i < L; ++i)
        DFS.emplace_back(-int(i) - 1);
}

/* =================================================================
 * Helpers
 * ================================================================= */
static bool Preprocess(string& inst, double thresh)
{
    ifstream file(inst);
    // Use dynamic resizing instead of hardcoded 1000000 to prevent crashes on huge items
    vector<unsigned long long> freq; 
    vector<unsigned long long> counted;

    if (file.good()) {
        string line; int ditem;
        while (getline(file, line) && give_time(clock() - start_time) < time_limit) {
            ++N;
            istringstream word(line);
            string itm;
            while (word >> itm) {
                ditem = stoi(itm);
                if (L < static_cast<unsigned int>(abs(ditem)))
                    L = static_cast<unsigned int>(abs(ditem));

                // Resize logic to prevent out-of-bounds access
                if (freq.size() < L) {
                    freq.resize(L, 0);
                    counted.resize(L, 0);
                }
                
                if (counted[abs(ditem) - 1] != N) {
                    ++freq[abs(ditem) - 1];
                    counted[abs(ditem) - 1] = N;
                }
            }
        }
    } else {
        cout << "!!!!!! No such file exists: " << inst << " !!!!!!\n";
        return false;
    }

    theta = (thresh < 1.0) ? ceil(thresh * N) : thresh;

    int real_L = 0;
    item_dic.assign(L, -1);
    for (unsigned int i = 0; i < L; ++i)
        if (freq[i] >= theta) item_dic[i] = ++real_L;

    if(b_disp) 
        cout << "Original number of items: " << L << "  Reduced to: " << real_L << '\n';

    L = real_L;
    N = 0;
    return true;
}

static void Load_items_pre(string& inst)
{
    ifstream file(inst);

    if (!file.good()) return;
    string line; int size_m, ditem; bool empty_seq = false;

    while (getline(file, line) && give_time(clock() - start_time) < time_limit) {
        vector<bool> counted(L, 0);
        istringstream word(line);

        if (!empty_seq) items.emplace_back();
        string itm; size_m = 0; bool sgn = false; empty_seq = true;

        while (word >> itm) {
            ditem = stoi(itm);

            // Safety check for dictionary bounds
            if (abs(ditem)-1 >= item_dic.size() || item_dic[abs(ditem) - 1] == -1) {
                if (!sgn) sgn = ditem < 0;
                continue;
            } else {
                ditem = (ditem > 0)
                      ?  item_dic[ditem - 1]
                      : -item_dic[-ditem - 1];
            }
            empty_seq = false;

            if (sgn) { if (ditem > 0) ditem = -ditem; sgn = false; }

            items.back().push_back(ditem);

            if (!counted[abs(ditem) - 1] && !just_build) {
                DFS[abs(ditem) - 1].seq_ID.push_back(items.size() - 1);
                DFS[abs(ditem) - 1].str_pnt.push_back(items.back().size() - 1);
                ++DFS[abs(ditem) - 1].freq;
                counted[abs(ditem) - 1] = true;
            }
            ++size_m;
        }
        if (empty_seq) continue;

        ++N;  E += size_m;  M = max<unsigned int>(M, static_cast<unsigned int>(size_m));
    }
}

static bool Load_items(string& inst)
{
    ifstream file(inst);
    if (!file.good()) {
        cout << "!!!!!! No such file exists: " << inst << " !!!!!!\n";
        return false;
    }

    string line; int size_m, ditem;
    while (getline(file, line) && give_time(clock() - start_time) < time_limit) {
        ++N;
        vector<bool> counted(L, 0); // Initial size based on current L
        istringstream word(line);

        items.emplace_back();
        string itm; size_m = 0;

        while (word >> itm) {
            ditem = stoi(itm);
            if (L < static_cast<unsigned int>(abs(ditem))) {
                L = static_cast<unsigned int>(abs(ditem));
                
                // 1. Grow global DFS structure
                while (DFS.size() < L) {
                    DFS.emplace_back(-int(DFS.size()) - 1);
                }
                // 2. CRASH FIX: Grow local counted vector immediately!
                // Without this, counted[ditem-1] writes to invalid memory.
                if (counted.size() < L) {
                    counted.resize(L, 0);
                }
            }
            items.back().push_back(ditem);

            if (!counted[abs(ditem) - 1] && !just_build) {
                DFS[abs(ditem) - 1].seq_ID.push_back(items.size() - 1);
                DFS[abs(ditem) - 1].str_pnt.push_back(items.back().size() - 1);
                ++DFS[abs(ditem) - 1].freq;
                counted[abs(ditem) - 1] = true;
            }
            ++size_m;
        }
        E += size_m;
        M = max<unsigned int>(M, static_cast<unsigned int>(size_m));
    }
    return true;
}

} // namespace largepp