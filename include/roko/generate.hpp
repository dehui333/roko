#ifndef ROKO_GENERATE_HPP_
#define ROKO_GENERATE_HPP_

#include <cstdint>
#include <unordered_map>
#include <memory>
#include <string>
#include <vector>

#include <Python.h>

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL gen_ARRAY_API
#include "numpy/arrayobject.h"

#include "roko/models.hpp"

namespace roko {

constexpr int dimensions[] = {200, 90};
constexpr int CENTER = dimensions[1] / 2;
constexpr int WINDOW = dimensions[1] / 3;
constexpr int MAX_INS = 3;
constexpr int REF_ROWS = 0;

struct Data {
  std::vector<std::vector<std::pair<std::int64_t, std::int64_t>>> positions; // these are the positions
                                                                             // the matrix columns correspond to
  std::vector<PyObject*> X; //Each is a pointer to an array of 2 dimensions(matrix), the sides given by dimensions. Type UINT8
};

struct pair_hash {
  template <class T1, class T2>
  std::size_t operator()(const std::pair<T1, T2>& p) const {
    const auto h1 = std::hash<T1>{}(p.first);
    const auto h2 = std::hash<T2>{}(p.second);

    return h1 ^ h2;
  }
};

// decouple from python binding
Data generate_features(const char*, const char*, const char*);

} // namespace roko

#endif // ROKO_GENERATE_HPP_
