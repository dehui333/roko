#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <type_traits>
#include <functional>
#include <algorithm>
#include <iostream>
#include <utility>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <set>

#include "roko/generate.hpp"

namespace roko {

namespace detail {

template <class T>
auto constexpr to_underlying(T t) noexcept {
  return static_cast<std::underlying_type_t<T>>(t);
}

} // namespace detail
/*
Inputs:
    filename: path to bam file mapping reads to ref(draft sequence)
    ref: draft sequence
    region: A string demarcating a region on ref. Format: NAME:START-END or NAME:START. The indices are 1-based and both ends inclusive.
            if END is not included, goes till end of ref.
Output:
    struct Data. Contains feature matrices and corresponding position indices indicating what positions the columns are associated with.

*/
Data generate_features(const char* filename, const char* ref,
                       const char* region) {
  auto bam = read_bam(filename); // Abstraction over the bam file 

  npy_intp dims[2]; // Required (int) type to create the python matrix object
  for (int i = 0; i < 2; i++) {
    dims[i] = dimensions[i];
  }

  std::vector<std::pair<std::int64_t, std::int64_t>> pos_queue; // Position indices on the reference

  std::unordered_map<std::pair<std::int64_t, std::int64_t>, 
                     std::unordered_map<std::uint32_t, BaseType>, pair_hash>
      align_info; // Map each position onto another map, which maps query/read id to bases(which are aligned to the position on the ref)

  std::unordered_map<uint32_t, std::pair<std::int64_t, std::int64_t>>
      align_bounds; // Map each query id to its first and last aligned position on the ref
  std::unordered_map<uint32_t, std::uint8_t> strand; // Map each query id to a bool indicating reverse strand or not

  auto data = Data();

  auto pileup_iter = bam->pileup(region); // iterator over positions on ref
  
  // Iterate through the positions, when associating each position with the bases aligned to it and add to queue to process.
  // When there is enough positions in the queue, build a feature matrix for one window of positions.
  while (pileup_iter->has_next()) {
    auto column = pileup_iter->next(); // A position object. 

    
    const std::int64_t rpos = column->position; // position index on the reference 
    
    //I think the checks are not required? Already included in next()?..
    if (rpos < pileup_iter->begin())
      // Not in the specified region yet
      continue;
    if (rpos >= pileup_iter->end())
      // Region has ended  
      break;


    // Associate position with the aligned bases of the reads at the position
    while (column->has_next()) { // Iterating the reads aligned to the region
      auto r = column->next(); // An alignment object

      if (r->is_refskip()) // For this read, no position on the read is aligned to the position on the ref
        continue;

      if (align_bounds.find(r->query_id()) == align_bounds.end()) {
        align_bounds.emplace(r->query_id(),
                             std::make_pair(r->ref_start(), r->ref_end()));
      } // Record the 1st and last aligned positions on the ref for the read, if not done yet.
      strand.emplace(r->query_id(), !r->rev()); // May be one read can have both directions being aligned?

      std::pair<std::int64_t, std::int64_t> index(rpos, 0); // The first element is the ref position of the Position object
      if (align_info.find(index) == align_info.end()) {
        pos_queue.emplace_back(rpos, 0);
      } // If this position on the ref has not been processed, add to queue? I think it makes sense to move this out of this loop?
        // It seems to be related to position on the ref, not the reads, which this loop is about.
      
      if (r->is_del()) {
        // DELETION
        align_info[index].emplace(r->query_id(), BaseType::GAP); // This position has a read indicating GAP  
      } else {
        // POSITION
        auto qbase = r->qbase(0);
        align_info[index].emplace(r->query_id(), qbase); 

        // INSERTION
        for (int i = 1, n = std::min(r->indel(), MAX_INS); i <= n; ++i) {// Take the number of insertions, capped at MAX_INS 
          index = std::pair<std::int64_t, std::int64_t>(rpos, i);        // add to allowance for insertions 

          if (align_info.find(index) == align_info.end()) {
            pos_queue.emplace_back(rpos, i); // Once one read has x indels, allocate allowance for x insertions, if <= MAX_INS
          }

          qbase = r->qbase(i);
          align_info[index].emplace(r->query_id(), qbase); // A read indicates there's an insertion at this position
        }
      }
    }

    // BUILD FEATURE MATRIX
    while (pos_queue.size() >= dimensions[1]) { // While there is enough positions for one matrix
      std::vector<std::uint32_t> valid_aligns; // contains read/query id
      const auto it = pos_queue.begin();

      for (auto s = 0; s < dimensions[1]; s++) {
        auto curr = it + s;
        for (auto& align : align_info[*curr]) { // For all the reads that is aligned with some base(s) in the region with the reference 
          if (align.second != BaseType::UNKNOWN) { // If has at least one not UNKNOWN
            valid_aligns.emplace_back(align.first); // Consider the read valid
          }
        }
      }
      
      // Removes duplicate valid read ids:
      // In containers that support equivalent keys,         
      // elements with equivalent keys are adjacent to each other in the iteration order of the container.
      valid_aligns.erase(std::unique(valid_aligns.begin(), valid_aligns.end()),
                         valid_aligns.end());
      valid_aligns.shrink_to_fit();

      auto X = PyArray_SimpleNew(2, dims, NPY_UINT8); //pointer to uninitialized array of 2 dimensions, the sides given by dims. Type UINT8
      uint8_t* value_ptr;

      // First handle assembly (REF_ROWS)
      // For one window worth of positions in the queue
      for (auto s = 0; s < dimensions[1]; s++) { // be fancy, decltype?
        auto curr = it + s;
        uint8_t value; // Base value on in the position for the ref

        if (curr->second != 0)
          value = detail::to_underlying(BaseType::GAP); // Positions after (n, 0) = gap
        else
          value = detail::to_underlying(get_base(ref[curr->first])); // (n, 0) is the nth base on the reference.

        for (int r = 0; r < REF_ROWS; r++) { // REF_ROWS = 0, this loop not used?
          value_ptr = (uint8_t*)PyArray_GETPTR2(X, r, s);
          *value_ptr = value; // Forward strand - no +6
        }
      }

      for (int r = REF_ROWS; r < dimensions[0]; r++) {
        // Sample reads that has at least one non unknown overlap in the window
        std::uint8_t base;
        auto const random_num = rand() % valid_aligns.size();
        std::uint32_t query_id = valid_aligns[random_num];

        auto& fwd = strand[query_id];

        auto it = pos_queue.begin();
        for (auto s = 0; s < dimensions[1]; s++) {
          auto curr = it + s;

          auto pos_itr = align_info[*curr].find(query_id); //find - Is the sampled read aligned at this position?
          auto& bounds = align_bounds[query_id];
          if (pos_itr == align_info[*curr].end()) { // If the sampled read is not aligned at this position
            if (curr->first < bounds.first || curr->first > bounds.second) {  // Aligned in the region but this position is not covered
              base = detail::to_underlying(BaseType::UNKNOWN); 
            } else { // Read has bases aligned to somewhere to the left and right of the current position, but not the current.
              base = detail::to_underlying(BaseType::GAP);
            }
          } else {
            // If aligned, get the base on the read at the aligned base
            base = detail::to_underlying(pos_itr->second);
          }

          value_ptr = (uint8_t*)PyArray_GETPTR2(X, r, s);
          *value_ptr = fwd ? base : (base + 6); // Insert value into matrix
        }
      }

      data.X.push_back(X); // Matrix for multiple windows
      data.positions.emplace_back(pos_queue.begin(),
                                  pos_queue.begin() + dimensions[1]); // positions the current matrix's columns correspond to 

      for (auto it = pos_queue.begin(), end = pos_queue.begin() + WINDOW; // Shift the sliding window covering positions on the ref 
           it != end; ++it) {                                             // may overlap with the previous
        align_info.erase(*it);
      }
      pos_queue.erase(pos_queue.begin(), pos_queue.begin() + WINDOW);
    }
  }

  return data;
}

} // namespace roko
