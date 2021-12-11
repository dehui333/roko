#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <unordered_set>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <utility>
#include <random>
#include <set>
#include <string>
#include <cmath>
#include <iostream>
#include "align_ins.h"
//#include "generate.h"
//#include "models.h"


// For reverse strand add +6
std::unordered_map<Bases, uint8_t, EnumClassHash> ENCODED_BASES = {
        {Bases::A, 0},
        {Bases::C, 1},
        {Bases::G, 2},
        {Bases::T, 3},
        {Bases::GAP, 4},
        {Bases::UNKNOWN, 5}
};


std::unique_ptr<Data> generate_features(const char* filename, const char* ref, const char* region) {
    auto bam = readBAM(filename);
	
    npy_intp dims[2];
    for (int i = 0; i < 2; i++) {
        dims[i] = dimensions[i];
    }


    std::vector<std::pair<long, long>> pos_queue;
    std::unordered_map<std::pair<long, long>, std::unordered_map<uint32_t, PosInfo>, pair_hash> align_info;
    std::unordered_map<uint32_t, std::pair<long, long>> align_bounds;
    std::unordered_map<uint32_t, bool> strand;

    auto data = std::unique_ptr<Data>(new Data());

    auto pileup_iter = bam->pileup(region);
    while (pileup_iter->has_next()) {
        auto column = pileup_iter->next();
             
        
        long rpos = column->position;
        if (rpos < pileup_iter->start()) continue;
        if (rpos >= pileup_iter->end()) break;
        std::vector<segment> ins_segments;
        //int longest_index = 0;
        bool col_has_ins = false;
        while(column->has_next()) {
            auto r = column->next();
            if (r->is_refskip()) continue;
            if (align_bounds.find(r->query_id()) == align_bounds.end()) {
                align_bounds.emplace(r->query_id(), std::make_pair(r->ref_start(), r->ref_end()));
            }
            strand.emplace(r->query_id(), !r->rev());
            std::pair<long, long> index(rpos, 0);
            if (align_info.find(index) == align_info.end()) {
                pos_queue.emplace_back(rpos, 0);
            }
            if (r->is_del()) {
                // DELETION
                align_info[index].emplace(r->query_id(), PosInfo(Bases::GAP));
            } else {
                // POSITION
                auto qbase = r->qbase(0);
                align_info[index].emplace(r->query_id(), PosInfo(qbase));
                // INSERTION
                std::string s;
                if (r->indel() > 0) {
                    s.reserve(r->indel());
                }
                for (int i = 1, n = r->indel(); i <= n; ++i) {
                  
                    qbase = r->qbase(i);
                    
                    s.push_back(base_to_char(qbase));                 
                }

                if (r->indel() > 0) {
                    col_has_ins = true;
                                 
                    ins_segments.emplace_back(s, r->indel(), r->query_id());

                    //if (ins_segments[longest_index].len < r->indel()) {

                      //  longest_index = ins_segments.size() - 1;

                    //}
                }
            }
            
            
        }
        if (col_has_ins) {
            align_ins_center_star(rpos, ins_segments, align_info, pos_queue);
            col_has_ins = false;
        }
        //BUILD FEATURE MATRIX
        while (pos_queue.size() >= dimensions[1]) {
            std::set<uint32_t> valid_aligns;
            const auto it = pos_queue.begin();

            for (auto s = 0; s < dimensions[1]; s++) {
                auto curr = it + s;
                for (auto& align : align_info[*curr]) {
                    if (align.second.base != Bases::UNKNOWN) {
                        valid_aligns.emplace(align.first);
                    }
                }
            }

            std::vector<uint32_t> valid(valid_aligns.begin(), valid_aligns.end());
            int valid_size = valid.size();

            auto X = PyArray_SimpleNew(2, dims, NPY_UINT8);
            uint8_t* value_ptr;

            // First handle assembly (REF_ROWS)
            for (auto s = 0; s < dimensions[1]; s++) {
                auto curr = it + s; uint8_t value;

                if (curr->second != 0) value = ENCODED_BASES[Bases::GAP];
                else value = ENCODED_BASES[get_base(ref[curr->first])];

                for (int r = 0; r < REF_ROWS; r++) {
                    value_ptr = (uint8_t*) PyArray_GETPTR2(X, r, s);
                    *value_ptr = value; // Forward strand - no +6
                }
            }

            for (int r = REF_ROWS; r < dimensions[0]; r++) {
                uint8_t base;
                auto random_num = rand() % valid_size;
                uint32_t query_id = valid[random_num];

                auto& fwd = strand[query_id];

                auto it = pos_queue.begin();
                for (auto s = 0; s < dimensions[1]; s++) {
                    auto curr = it + s;

                    auto pos_itr = align_info[*curr].find(query_id);
                    auto& bounds = align_bounds[query_id];
                    if (pos_itr == align_info[*curr].end()) {
                        if (curr->first < bounds.first || curr->first > bounds.second) {
                            base = ENCODED_BASES[Bases::UNKNOWN];
                        } else {
                            base = ENCODED_BASES[Bases::GAP];
                        }
                    } else {
                        base = ENCODED_BASES[pos_itr->second.base];
                    }

                    value_ptr = (uint8_t*) PyArray_GETPTR2(X, r, s);
                    *value_ptr = fwd ? base : (base + 6);
	
                }
            }
			
            data->X.push_back(X);
            data->positions.emplace_back(pos_queue.begin(), pos_queue.begin() + dimensions[1]);

            for (auto it = pos_queue.begin(), end = pos_queue.begin() + WINDOW; it != end; ++it) {
                align_info.erase(*it);
            }
            pos_queue.erase(pos_queue.begin(), pos_queue.begin() + WINDOW);
        }
    }

    return data;
}


