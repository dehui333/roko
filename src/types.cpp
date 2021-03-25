#include <stdexcept>

#include "roko/types.hpp"

namespace roko {

unique_htsFile create_unique_htsFile(const char* htsfile_name) {
  const auto file_handle = ::hts_open(htsfile_name, "rb");
  if (file_handle == nullptr) {
    throw std::runtime_error("[roko::create_unique_htsFile] failed to open: " +
                             std::string(htsfile_name));
  }

  return unique_htsFile(file_handle, &::hts_close);
}

unique_hts_idx_t create_unique_hts_idx_t(::htsFile* htsfile_handle,
                                         const char* htsfile_path) {
  const auto index_handle = ::sam_index_load(htsfile_handle, htsfile_path);
  if (index_handle == nullptr) {
    throw std::runtime_error(
        "[roko::create_unique_hts_idx_t] failed to index: " +
        std::string(htsfile_path));
  }

  return unique_hts_idx_t(index_handle, &::hts_idx_destroy);
}

unique_bam_hdr_t create_unique_bam_hdr_t(::htsFile* htsfile_handle, const char* htsfile_path) {
  const auto header_handle = ::sam_hdr_read(htsfile_handle);
  if (header_handle == nullptr) {
    throw std::runtime_error(
        "[roko::create_unique_bam_hdr_t] failed to read bam header from:"
        + std::string(htsfile_path));
  }

  return unique_bam_hdr_t(header_handle, &bam_hdr_destroy);
}

} // namespace roko