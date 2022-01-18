#include<vector>
#include<string>

#include "models.h"
#include "generate.h"


struct segment {
    std::string sequence;
    int len;
    int index;
    segment(std::string seq, int l, uint32_t id) : sequence(seq), len(l), index(id) {};
    
};

Bases char_to_base(char c);

char base_to_char(Bases b);

char int_to_char(uint8_t i);

Bases int_to_base(uint8_t i);

uint8_t char_to_int(char c);

void align_center_star(long base_index, std::vector<segment>& segments, int star_index,
    std::unordered_map<std::pair<long, long>, std::unordered_map<uint32_t, PosInfo>, pair_hash>& align_info, 
    std::vector<std::pair<long, long>>& pos_queue, unsigned int& threshold_num, std::vector<uint32_t>& no_ins_reads,
    std::unordered_map<std::pair<long, long>, uint8_t, pair_hash>& labels_info);

void align_ins_longest_star(long base_index, std::vector<segment>& ins_segments, int longest_index, 
    std::unordered_map<std::pair<long, long>, std::unordered_map<uint32_t, PosInfo>, pair_hash>& align_info, 
    std::vector<std::pair<long, long>>& pos_queue, unsigned int& threshold_num, std::vector<uint32_t>& no_ins_reads);

void align_ins_center_star(long base_index, std::vector<segment>& ins_segments,
    std::unordered_map<std::pair<long, long>, std::unordered_map<uint32_t, PosInfo>, pair_hash>& align_info, 
    std::vector<std::pair<long, long>>& pos_queue, unsigned int& threshold_num, std::vector<uint32_t>& no_ins_reads, 
    std::unordered_map<std::pair<long, long>, uint8_t, pair_hash>& labels_info);

int find_center(std::vector<segment>& segments);
