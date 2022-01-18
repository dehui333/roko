#include <iostream>
#include <unordered_map>

#include "align_ins.h"
#include "edlib.h"
#include "generate.h"


Bases char_to_base(char c) {
    switch (c) {
        case 'A':
            return Bases::A;
        case 'C':
            return Bases::C;
        case 'G':
            return Bases::G;
        case 'T':
            return Bases::T;
        case '*':
            return Bases::GAP;
        case 'N':
            std::cout << "Unknown base!" << std::endl;
            return Bases::UNKNOWN;
        default:
            std::cout << "Non N unknown!" << std::endl;
            return Bases::UNKNOWN;
    }
}

char base_to_char(Bases b) {
    switch (b) {
        case Bases::A:
            return 'A';
        case Bases::C:
            return 'C';
        case Bases::G:
            return 'G';
        case Bases::T:
            return 'T';
        case Bases::GAP:
            std::cout << "Gap!" << std::endl;
            return '*';
        default:
            std::cout << "Unknown base!" << std::endl;
            return 'N';
    }
}

char int_to_char(uint8_t i) {
    switch (i) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        case 4:
            return '*';
        default:
            return 'N';
    }
}

Bases int_to_base(uint8_t i) {

    switch (i) {
        case 0:
            return Bases::A;
        case 1:
            return Bases::C;
        case 2:
            return Bases::G;
        case 3:
            return Bases::T;
        case 4:
            return Bases::GAP;
        default:
            return Bases::UNKNOWN;
    }
}

uint8_t char_to_int(char c) {
    switch (c) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        case '*':
            return 4;
        default:
            return 5;
    }
}


void align_center_star(long base_index, std::vector<segment>& segments, int star_index, 
        std::unordered_map<std::pair<long, long>, std::unordered_map<uint32_t, PosInfo>, pair_hash>& align_info, 
        std::vector<std::pair<long, long>>& pos_queue, unsigned int& threshold_num, std::vector<uint32_t>& no_ins_reads, 
        std::unordered_map<std::pair<long, long>, uint8_t, pair_hash>& labels_info) {

    std::vector<uint32_t> seq_indices(segments.size());
    segment star = segments[star_index];
    std::unordered_map<uint32_t, PosInfo> star_positions[star.len]; //stores bases aligned to original positions on the star
    std::vector<std::unordered_map<uint32_t, PosInfo>> ins_positions[star.len+1]; //stores bases aligned to gaps inserted into the original seq of star
    uint8_t star_positions_labels[star.len];
    for (int i = 0; i < star.len; i ++) {
        star_positions_labels[i] = 4;
    }
    std::vector<uint8_t> ins_positions_labels[star.len+1];
    int total_ins_pos = 0;
    for (auto s: segments) {
        //std::cout << s.index << " " << s.sequence << std::endl;
        if (s.index != -1) seq_indices.push_back(s.index);
        if (s.index != star.index) {
            EdlibAlignResult result = edlibAlign(s.sequence.c_str(), s.len, star.sequence.c_str(),
                    star.len, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
            int ref_pos = -1; // pointing to before next to read ref base
            int query_pos = -1; // pointing to before next to read query base 
            unsigned int ins_index = 0; // index of next insertion, 0-based
            char char_at_pos;
            Bases base_at_pos;
            for (int i = 0; i < result.alignmentLength; i++) {
                switch (result.alignment[i]) {
                    case 0: // match
                        ins_index = 0;	      
                        char_at_pos = s.sequence[++query_pos];
                        base_at_pos = char_to_base(char_at_pos);
                        ref_pos++;
                        if (s.index == -1) {
                            star_positions_labels[ref_pos] = char_to_int(char_at_pos);
                        } else {
                            star_positions[ref_pos].emplace(s.index, PosInfo(base_at_pos));
                        }
                        break;
                    case 1: // insertion

                        char_at_pos = s.sequence[++query_pos];
                        base_at_pos = char_to_base(char_at_pos);
                        if (ins_positions[ref_pos+1].size() < ins_index + 1) { // if not enough maps to record bases in that position
                            ins_positions[ref_pos+1].push_back(std::unordered_map<uint32_t, PosInfo>{});
                            ins_positions_labels[ref_pos+1].push_back(4);
                            total_ins_pos++;
                        }
                        if (s.index == -1) {
                            ins_positions_labels[ref_pos+1][ins_index] = char_to_int(char_at_pos);
                        } else {
                            ins_positions[ref_pos+1][ins_index].emplace(s.index, PosInfo(base_at_pos));
                        }
                        ins_index++;
                        break;
                    case 2: // deletion

                        ins_index = 0;
                        ref_pos++;
                        break;
                    case 3: // mismatch

                        ins_index = 0; 
                        char_at_pos = s.sequence[++query_pos];
                        base_at_pos = char_to_base(char_at_pos);
                        ref_pos++;
                        if (s.index == -1) {
                            star_positions_labels[ref_pos] = char_to_int(char_at_pos);
                        } else {
                            star_positions[ref_pos].emplace(s.index, PosInfo(base_at_pos));
                        }
                        break;
                    default:
                        std::cout << "Uknown alignment result!\n";


                }
            }

            edlibFreeAlignResult(result);

        } else {
            // record bases on the star    
            for (int i = 0; i < s.len; i++) {
                const char char_at_pos = s.sequence[i];
                Bases base_at_pos = char_to_base(char_at_pos);               
                if (s.index == -1) {
                    star_positions_labels[i] = char_to_int(char_at_pos);
                } else {
                    star_positions[i].emplace(s.index, PosInfo(base_at_pos));
                }
            }

        }


    }

    long count = 1;
    for (unsigned int i = 0; i < ins_positions[0].size(); i++) {
        auto& map = ins_positions[0][i];
        auto index = std::pair<long, long>(base_index, count);
        if (map.size() >= threshold_num) {
            pos_queue.emplace_back(base_index, count);
            count++;
        }

        for (auto& id: seq_indices) {
            if (map.find(id) == map.end()) {
                map.emplace(id, PosInfo(Bases::GAP));	
            }	  
        }
        for (auto& id: no_ins_reads) {
            
            map.emplace(id, PosInfo(Bases::GAP));	  
            
        }

        align_info[index] = map;
        labels_info[index] = ins_positions_labels[0][i];
        
    }
    for (int i = 0; i < star.len; i++) {
        auto index = std::pair<long, long>(base_index, count);

        if (star_positions[i].size() >= threshold_num) {
            pos_queue.emplace_back(base_index, count);
            count++;
        }

        for (auto& id: seq_indices) {
            if (star_positions[i].find(id) == star_positions[i].end()) {
                star_positions[i].emplace(id, PosInfo(Bases::GAP));	
            }	  
        }
        for (auto& id: no_ins_reads) {
            star_positions[i].emplace(id, PosInfo(Bases::GAP));	  
        }

        align_info[index] = star_positions[i];
        labels_info[index] = star_positions_labels[i]; 

        for (unsigned int j = 0; j < ins_positions[i+1].size(); j++) {
            auto& map = ins_positions[i+1][j];
            auto index = std::pair<long, long>(base_index, count);

            if (map.size() >= threshold_num) {
                pos_queue.emplace_back(base_index, count);
                count++;
            }



            for (auto& id: seq_indices) {
                if (map.find(id) == map.end()) {
                    map.emplace(id, PosInfo(Bases::GAP));	
                }	  
            }
            for (auto& id: no_ins_reads) {
                map.emplace(id, PosInfo(Bases::GAP));	  
            }

            align_info[index] = map;
            labels_info[index] = ins_positions_labels[i+1][j];
        }
    }


    /*
       std::cout << "pos -1" << std::endl;
       for (int i = 0; i < ins_positions[0].size(); i++) {
       std::cout << "ins " << i << ": ";
       for (auto p: ins_positions[0][i]) {
       std::cout << p.first << ":" << base_to_char(p.second.base) << " ";
       }
       std::cout << std::endl;

       }
       std::cout << std::endl;
       for (int i = 0; i < star.len; i++) {
       std::cout << "pos " << i << std::endl;
       for (auto p: star_positions[i]) {
       std::cout << p.first << ":" << base_to_char(p.second.base) << " ";              
       }
       std::cout << std::endl;

       for (int j = 0; j < ins_positions[i+1].size(); j++) {
       std::cout << "ins " << j << ": ";
       for (auto p: ins_positions[i+1][j]) {
       std::cout << p.first << ":" << base_to_char(p.second.base) << " ";
       }
       std::cout << std::endl;

       }
       std::cout << std::endl;
       } */

}

int find_center(std::vector<segment>& segments) {
    int dists[segments.size()]{0};
    for (unsigned int i = 0; i < segments.size(); i++) {
        for (unsigned int j = i + 1; j < segments.size(); j++) {

            EdlibAlignResult result = edlibAlign(segments[i].sequence.c_str(), segments[i].len, segments[j].sequence.c_str(),
                    segments[j].len, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
            dists[i] += result.editDistance;
            dists[j] += result.editDistance;
            edlibFreeAlignResult(result);
        }
    }
    int best_pos_index = 0;
    for (unsigned int i = 0; i < segments.size(); i++) {
        if (dists[i] < dists[best_pos_index]) {
            best_pos_index = i;
        }
    }
    return best_pos_index;    

}

void align_ins_longest_star(long base_index, std::vector<segment>& ins_segments, int longest_index,
        std::unordered_map<std::pair<long, long>, std::unordered_map<uint32_t, PosInfo>, pair_hash>& align_info, 
        std::vector<std::pair<long, long>>& pos_queue, unsigned int& threshold_num, std::vector<uint32_t>& no_ins_reads) {
    std::unordered_map<std::pair<long, long>, uint8_t, pair_hash> dummy;
    align_center_star(base_index, ins_segments, longest_index, align_info, pos_queue, threshold_num, no_ins_reads, dummy);

}
void align_ins_center_star(long base_index, std::vector<segment>& ins_segments,
        std::unordered_map<std::pair<long, long>, std::unordered_map<uint32_t, PosInfo>, pair_hash>& align_info, 
        std::vector<std::pair<long, long>>& pos_queue, unsigned int& threshold_num, std::vector<uint32_t>& no_ins_reads,
        std::unordered_map<std::pair<long, long>, uint8_t, pair_hash>& labels_info) {
    int center_index = find_center(ins_segments);
    align_center_star(base_index, ins_segments, center_index, align_info, pos_queue, threshold_num, no_ins_reads, labels_info);

}

