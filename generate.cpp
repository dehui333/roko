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


std::unordered_map<std::pair<long, long>, uint8_t, pair_hash> convert_py_labels_dict(PyObject *dict) {

    Py_ssize_t pos = 0;
    PyObject *key = NULL;
    PyObject *value = NULL;
    std::unordered_map<std::pair<long, long>, uint8_t, pair_hash> map;    

    if (! PyDict_Check(dict)) {
        PyErr_Format(PyExc_TypeError, 
                "Argument \"dict\" to %s must be dict not \"%s\"", 
                __FUNCTION__, Py_TYPE(dict)->tp_name);
        return map;	       
    }
    while (PyDict_Next(dict, &pos, &key, &value)) {
        if (! PyTuple_Check(key)) {
            PyErr_SetString(PyExc_TypeError, "A key of dict is not a tuple!");
            map.clear();
            return map;
        } 
        if (PyTuple_Size(key) != static_cast<Py_ssize_t>(2)) {
            PyErr_SetString(PyExc_TypeError, "A tuple of dict is not a pair!");
            map.clear();
            return map;
        }
        PyObject *pair_item0 = PyTuple_GetItem(key, 0);
        PyObject *pair_item1 = PyTuple_GetItem(key, 1);
        if ((!PyLong_Check(pair_item0)) || (!PyLong_Check(pair_item1))) {
            PyErr_SetString(PyExc_TypeError, "A tuple of dict does contain two longs!");
            map.clear();
            return map;
        }
        if (! PyLong_Check(value)) {
            PyErr_SetString(PyExc_TypeError, "A value of dict is not of long type!");
            map.clear();
            return map;
        }
        long pair_item0_c = PyLong_AsLong(pair_item0);
        long pair_item1_c = PyLong_AsLong(pair_item1);
        uint8_t value_c = PyLong_AsLong(value);
        if (PyErr_Occurred()) {
            map.clear();
            return map;
        }
        map.emplace(std::make_pair(pair_item0_c, pair_item1_c), value_c);

    }

    return map;
}

std::unique_ptr<Data> generate_features(const char* filename, const char* ref, const char* region,
        std::unordered_map<std::pair<long, long>, uint8_t, pair_hash>& labels, int use_labels) {
    auto bam = readBAM(filename);

    npy_intp dims[2];
    npy_intp labels_dim[1];
    labels_dim[0] = dimensions[1];
    for (int i = 0; i < 2; i++) {
        dims[i] = dimensions[i];
    }

   /* if (use_labels) {
        for (auto& item:labels) {
            auto index = item.first;
            auto label = item.second;
            std::cout << index.first << ", " << index.second << " is " << int_to_char(label) << std::endl;
        }

    }*/

    std::vector<std::pair<long, long>> pos_queue;
    std::unordered_map<std::pair<long, long>, std::unordered_map<uint32_t, PosInfo>, pair_hash> align_info;
    std::unordered_map<std::pair<long, long>, uint8_t, pair_hash> labels_info;
    std::unordered_map<uint32_t, std::pair<long, long>> align_bounds;
    std::unordered_map<uint32_t, bool> strand;

    auto data = std::unique_ptr<Data>(new Data());

    auto pileup_iter = bam->pileup(region);
    while (pileup_iter->has_next()) {
        auto column = pileup_iter->next();
        long rpos = column->position;
        if (rpos < pileup_iter->start()) continue;
        if (rpos >= pileup_iter->end()) break;
        unsigned int threshold_num = threshold_prop * column->count();
        if (threshold_num == 0) threshold_num = 1;	
        std::vector<segment> ins_segments;
        std::vector<uint32_t> no_ins_reads;
        //int longest_index = 0;
        bool col_has_enough_ins = false;
        unsigned int total_ins_len = 0;
        unsigned int max_indel = 0;
        std::string s;
        if (use_labels) {
            std::pair<long, long> index {rpos, 0};
            labels_info[index] = labels[index];
            long ins_count = 1;
            index = std::make_pair(rpos, ins_count);
            auto found = labels.find(index);
            while (found != labels.end()) {
                char c = int_to_char(labels[index]);
                s.push_back(c);
                ins_count++;
                index = std::make_pair(rpos, ins_count);
                found = labels.find(index);
            }
        }
        if (s.size() > 0) {
            ins_segments.emplace_back(s, s.size(), -1);\
        } else {
            no_ins_reads.push_back(-1);
        }
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
                if (r-> indel() > 0) {
                    std::string s;
                    s.reserve(r->indel());
                    for (int i = 1, n = r->indel(); i <= n; ++i) {
                        qbase = r->qbase(i);
                        s.push_back(base_to_char(qbase));                 
                    }
                    if (static_cast<unsigned int>(r->indel()) > max_indel) max_indel = r->indel();
                    total_ins_len += r->indel();     
                    ins_segments.emplace_back(s, r->indel(), static_cast<int>(r->query_id()));
                } else {
                    no_ins_reads.push_back(r->query_id());

                }

                /*std::string s;
                if (r->indel() > 0) {
                    s.reserve(r->indel());
                }
                for (int i = 1, n = r->indel(); i <= n; ++i) {

                    qbase = r->qbase(i);

                    s.push_back(base_to_char(qbase));                 
                }
                if (r->indel() > 0) {
                    if (static_cast<unsigned int>(r->indel()) > max_indel) max_indel = r->indel();
                    total_ins_len += r->indel();     
                    ins_segments.emplace_back(s, r->indel(), static_cast<int>(r->query_id()));

                    //if (ins_segments[longest_index].len < r->indel()) {

                    //  longest_index = ins_segments.size() - 1;

                    //}
                } else {
                    no_ins_reads.push_back(r->query_id());		
                }*/
            }


        }
        float avg_ins_len = total_ins_len/ (float) column->count();
        if (avg_ins_len   >= align_len_threshold && ins_segments.size() > 0) col_has_enough_ins = true;

        if (col_has_enough_ins) {
            align_ins_center_star(rpos, ins_segments, align_info, pos_queue, threshold_num, no_ins_reads, labels_info);
            col_has_enough_ins = false;

        } else {

            for (auto& s: ins_segments) {
                long count = 1;    
                for (auto& c: s.sequence) {		
                    std::pair<long, long> index(rpos, count);
                    align_info[index].emplace(s.index, PosInfo(char_to_base(c)));
                    if (align_info[index].size() == threshold_num) {
                        pos_queue.emplace_back(rpos, count);	
                    }
                    count++;
                }	
            }

            for (auto& s: ins_segments) {
                for (long i = s.len + 1; i <= max_indel; i++) {
                    std::pair<long, long> index(rpos, i);
                    align_info[index].emplace(s.index, PosInfo(Bases::GAP));		   

                }	

            }	    

            for (auto& id: no_ins_reads) {
                for (long i = 1; i <= max_indel; i++) {
                    std::pair<long, long> index(rpos, i);
                    align_info[index].emplace(id, PosInfo(Bases::GAP));		   

                }	
            }		    

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
            auto Y = PyArray_SimpleNew(1, labels_dim, NPY_UINT8);
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
                auto random_n = rand();
                auto random_num = random_n  % valid_size;

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

            for (auto s = 0; s < dimensions[1]; s++) {
                auto curr = it + s;
                uint8_t value = labels_info[*curr];
                value_ptr = (uint8_t*) PyArray_GETPTR1(Y, s);
                *value_ptr = value;
            }

            data->X.push_back(X);
            data->Y.push_back(Y);
            data->positions.emplace_back(pos_queue.begin(), pos_queue.begin() + dimensions[1]);

            for (auto it = pos_queue.begin(), end = pos_queue.begin() + WINDOW; it != end; ++it) {
                align_info.erase(*it);
            }
            pos_queue.erase(pos_queue.begin(), pos_queue.begin() + WINDOW);

        }
    }

    return data;
}


