#include<vector>
#include<string>

#include "models.h"
#include "generate.h"


struct segment {
    std::string sequence;
    int len;
    uint32_t index;
    segment(std::string seq, int l, uint32_t id) : sequence(seq), len(l), index(id) {};
    
};

Bases char_to_base(char c);
char base_to_char(Bases b);
void align_center_star(long base_index, std::vector<segment>& segments, int star_index,
    std::unordered_map<std::pair<long, long>, std::unordered_map<uint32_t, PosInfo>, pair_hash>& align_info, 
    std::vector<std::pair<long, long>>& pos_queue);
void align_ins_longest_star(long base_index, std::vector<segment>& ins_segments, int longest_index, 
    std::unordered_map<std::pair<long, long>, std::unordered_map<uint32_t, PosInfo>, pair_hash>& align_info, 
    std::vector<std::pair<long, long>>& pos_queue);
void align_ins_center_star(long base_index, std::vector<segment>& ins_segments,
    std::unordered_map<std::pair<long, long>, std::unordered_map<uint32_t, PosInfo>, pair_hash>& align_info, 
    std::vector<std::pair<long, long>>& pos_queue);
int find_center(std::vector<segment>& segments);
