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
        default:
            std::cout << "Unknown base type!" << std::endl;
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
        default:
            std::cout << "Unknown base type!" << std::endl;
            return '*';
    }
}

void align_center_star(long base_index, std::vector<segment>& segments, int star_index, 
        std::unordered_map<std::pair<long, long>, std::unordered_map<uint32_t, PosInfo>, pair_hash>& align_info, 
        std::vector<std::pair<long, long>>& pos_queue) {

    
    segment star = segments[star_index];
    std::unordered_map<uint32_t, PosInfo> star_positions[star.len];
       // Insertions at position x means they follows position x of star. 
    // Insertions at position -1 means they are left of position 0 of star. 
    std::vector<std::unordered_map<uint32_t, PosInfo>> ins_positions[star.len+1]; 
    int total_ins_pos = 0;
    for (auto s: segments) {
        //std::cout << s.index << " " << s.sequence << std::endl;
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
                  case 0:
	              ins_index = 0;	      
                      char_at_pos = s.sequence[++query_pos];
                      base_at_pos = char_to_base(char_at_pos);
                      star_positions[++ref_pos].emplace(s.index, PosInfo(base_at_pos));                   
                      break;
                  case 1: 
                      char_at_pos = s.sequence[++query_pos];
                      base_at_pos = char_to_base(char_at_pos);
                      if (ins_positions[ref_pos+1].size() < ins_index + 1) {
			      
                          ins_positions[ref_pos+1].push_back(std::unordered_map<uint32_t, PosInfo>{});
                          total_ins_pos++;
                      }

                      ins_positions[ref_pos+1][ins_index++].emplace(s.index, PosInfo(base_at_pos));
                      break;
                  case 2: 
		      ins_index = 0;
                      ref_pos++;
                      break;
                  case 3:
		      ins_index = 0; 
                      char_at_pos = s.sequence[++query_pos];
                      base_at_pos = char_to_base(char_at_pos);
                      star_positions[++ref_pos].emplace(s.index, PosInfo(base_at_pos));
                      break;
                  default:
                    std::cout << "Uknown alignment result!\n";
                  
                  
              }
            }
            
            edlibFreeAlignResult(result);
                
      } else {
	      
          for (int i = 0; i < s.len; i++) {
              const char char_at_pos = s.sequence[i];
              Bases base_at_pos = char_to_base(char_at_pos);
              star_positions[i].emplace(s.index, PosInfo(base_at_pos));
          }
          
      }
      
      
  }
  
  long count = 0;
  for (auto map: ins_positions[0]) {
      auto index = std::pair<long, long>(base_index, ++count);
      pos_queue.emplace_back(base_index, count);
      align_info[index] = map;
  }
  for (int i = 0; i < star.len; i++) {
      auto index = std::pair<long, long>(base_index, ++count);
      pos_queue.emplace_back(base_index, count);
      align_info[index] = star_positions[i];
      for (auto map:ins_positions[i+1]) {
          auto index = std::pair<long, long>(base_index, ++count);
          pos_queue.emplace_back(base_index, count);
          align_info[index] = map;
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
        std::vector<std::pair<long, long>>& pos_queue) {
    align_center_star(base_index, ins_segments, longest_index, align_info, pos_queue);
    
}
void align_ins_center_star(long base_index, std::vector<segment>& ins_segments,
        std::unordered_map<std::pair<long, long>, std::unordered_map<uint32_t, PosInfo>, pair_hash>& align_info, 
        std::vector<std::pair<long, long>>& pos_queue) {

    int center_index = find_center(ins_segments);
    align_center_star(base_index, ins_segments, center_index, align_info, pos_queue);
    
}

