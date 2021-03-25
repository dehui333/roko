// author: tbrekalo

#ifndef ROKO_TYPES_HPP_
#define ROKO_TYPES_HPP_

/**
 * @brief support for htslib types
 */

#include <memory>

#include "htslib/sam.h"

namespace roko {

/* clang-format off */
using unique_htsFile = std::unique_ptr<::htsFile, decltype(&::hts_close)>;
using unique_hts_idx_t = std::unique_ptr<::hts_idx_t, decltype(&::hts_idx_destroy)>;
using unique_bam_hdr_t = std::unique_ptr<::bam_hdr_t, decltype(&::bam_hdr_destroy)>;
/* clang-format on */

/**
 * @brief Create a unique htsFile object
 * 
 * @param htsfile_path path to targeted file
 *
 * @throws std::runtime_error on a failed read
 *
 * @return unique_htsFile wrapped htsFile* on a successful read
 */
unique_htsFile create_unique_htsFile(const char* htsfile_path);

/**
 * @brief Create a unique hts idx t object
 * 
 * @param htsfile_handle htsFile* for a targeted file
 * @param htsfile_path  path to a tergeted file
 *
 * @throws std::runtime_error on a failed read
 *
 * @return unique_hts_idx_t wrapped hts_idx_t* on a successful read
 */
unique_hts_idx_t create_unique_hts_idx_t(::htsFile* htsfile_handle,
                                         const char* htsfile_path);

/**
 * @brief Create a unique bam hdr t object
 * 
 * @param htsfile_handle htsFile* for a targeted file
 * @param htsfile_path path to a targeted file
 *
 * @throws std::runtime_error on a failed read
 *
 * @return unique_bam_hdr_t wrapped bam_hdr_t on a successful read
 */
unique_bam_hdr_t create_unique_bam_hdr_t(::htsFile* htsfile_handle,
                                         const char* htsfile_path);

} // namespace roko

#endif /* ROKO_TYPES_HPP_ */