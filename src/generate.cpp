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

Data generate_features(const char* filename, const char* ref,
                       const char* region) {
  auto bam = read_bam(filename);

  npy_intp dims[2];
  for (int i = 0; i < 2; i++) {
    dims[i] = dimensions[i];
  }

  std::vector<std::pair<std::int64_t, std::int64_t>> pos_queue;

  // positions? alignment for a position?
  std::unordered_map<std::pair<std::int64_t, std::int64_t>, // ?!?!?
                     std::unordered_map<std::uint32_t, BaseType>, pair_hash>
      align_info;

  std::unordered_map<uint32_t, std::pair<std::int64_t, std::int64_t>>
      align_bounds;
  std::unordered_map<uint32_t, std::uint8_t> strand;

  auto data = Data();

  auto pileup_iter = bam->pileup(region);
  while (pileup_iter->has_next()) {
    auto column = pileup_iter->next();

    
    const std::int64_t rpos = column->position;
    if (rpos < pileup_iter->begin())
      continue;
    if (rpos >= pileup_iter->end())
      break;

    while (column->has_next()) {
      auto r = column->next();

      if (r->is_refskip())
        continue;

      if (align_bounds.find(r->query_id()) == align_bounds.end()) {
        align_bounds.emplace(r->query_id(),
                             std::make_pair(r->ref_start(), r->ref_end()));
      }
      strand.emplace(r->query_id(), !r->rev());

      std::pair<std::int64_t, std::int64_t> index(rpos, 0);
      if (align_info.find(index) == align_info.end()) {
        pos_queue.emplace_back(rpos, 0);
      }

      if (r->is_del()) {
        // DELETION
        align_info[index].emplace(r->query_id(), BaseType::GAP);
      } else {
        // POSITION
        auto qbase = r->qbase(0);
        align_info[index].emplace(r->query_id(), qbase);

        // INSERTION
        for (int i = 1, n = std::min(r->indel(), MAX_INS); i <= n; ++i) {
          index = std::pair<std::int64_t, std::int64_t>(rpos, i);

          if (align_info.find(index) == align_info.end()) {
            pos_queue.emplace_back(rpos, i);
          }

          qbase = r->qbase(i);
          align_info[index].emplace(r->query_id(), qbase);
        }
      }
    }

    // BUILD FEATURE MATRIX
    while (pos_queue.size() >= dimensions[1]) {
      std::vector<std::uint32_t> valid_aligns;
      const auto it = pos_queue.begin();

      for (auto s = 0; s < dimensions[1]; s++) {
        auto curr = it + s;
        for (auto& align : align_info[*curr]) {
          if (align.second != BaseType::UNKNOWN) {
            valid_aligns.emplace_back(align.first);
          }
        }
      }

      valid_aligns.erase(std::unique(valid_aligns.begin(), valid_aligns.end()),
                         valid_aligns.end());
      valid_aligns.shrink_to_fit();

      auto X = PyArray_SimpleNew(2, dims, NPY_UINT8);
      uint8_t* value_ptr;

      // First handle assembly (REF_ROWS)
      for (auto s = 0; s < dimensions[1]; s++) { // be fancy, decltype?
        auto curr = it + s;
        uint8_t value;

        if (curr->second != 0)
          value = detail::to_underlying(BaseType::GAP);
        else
          value = detail::to_underlying(get_base(ref[curr->first]));

        for (int r = 0; r < REF_ROWS; r++) {
          value_ptr = (uint8_t*)PyArray_GETPTR2(X, r, s);
          *value_ptr = value; // Forward strand - no +6
        }
      }

      for (int r = REF_ROWS; r < dimensions[0]; r++) {
        std::uint8_t base;
        auto const random_num = rand() % valid_aligns.size();
        std::uint32_t query_id = valid_aligns[random_num];

        auto& fwd = strand[query_id];

        auto it = pos_queue.begin();
        for (auto s = 0; s < dimensions[1]; s++) {
          auto curr = it + s;

          auto pos_itr = align_info[*curr].find(query_id);
          auto& bounds = align_bounds[query_id];
          if (pos_itr == align_info[*curr].end()) {
            if (curr->first < bounds.first || curr->first > bounds.second) {
              base = detail::to_underlying(BaseType::UNKNOWN);
            } else {
              base = detail::to_underlying(BaseType::GAP);
            }
          } else {
            base = detail::to_underlying(pos_itr->second);
          }

          value_ptr = (uint8_t*)PyArray_GETPTR2(X, r, s);
          *value_ptr = fwd ? base : (base + 6);
        }
      }

      data.X.push_back(X);
      data.positions.emplace_back(pos_queue.begin(),
                                  pos_queue.begin() + dimensions[1]);

      for (auto it = pos_queue.begin(), end = pos_queue.begin() + WINDOW;
           it != end; ++it) {
        align_info.erase(*it);
      }
      pos_queue.erase(pos_queue.begin(), pos_queue.begin() + WINDOW);
    }
  }

  return data;
}

} // namespace roko
