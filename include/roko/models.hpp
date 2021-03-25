//
// Created by dominik on 10. 10. 2019..
//

#ifndef ROKO_MODELS_HPP_
#define ROKO_MODELS_HPP_

#include <cstdint>
#include <memory>
#include <string>

#include "htslib/sam.h"
#include "types.hpp"

namespace roko {

struct PileupData {
  htsFile* file;
  bam_hdr_t* header;
  hts_itr_t* iter;
};

std::int32_t iter_bam(void* data, bam1_t* b);

// what does min mapping quality mean?
constexpr std::uint8_t MIN_MAPPING_QUALITY = 10;
constexpr std::uint16_t FILTER_FLAG =
    BAM_FUNMAP | BAM_FDUP | BAM_FQCFAIL | BAM_FSUPPLEMENTARY | BAM_FSECONDARY;

enum class BaseType : std::uint8_t { A, C, G, T, GAP, UNKNOWN };

BaseType get_base(char b) noexcept;

// TODO: meaing of the members?
struct RegionInfo {
  RegionInfo(std::string, std::int32_t, std::int32_t);

  const std::string name;
  const std::int32_t begin;
  const std::int32_t end;
};

class Position;
class Alignment;
class PositionIterator;
std::unique_ptr<RegionInfo> parse_region(const std::string&);

/**
 * @brief hts file wrapper with index information and bam header
 */
class BAMFile {
public:
  BAMFile(unique_htsFile, unique_hts_idx_t, unique_bam_hdr_t);
  ~BAMFile() = default;

  BAMFile(BAMFile const&) = delete;
  BAMFile& operator=(BAMFile const&) = delete;

  BAMFile(BAMFile&&) = default;
  BAMFile& operator=(BAMFile&&) = default;

  std::unique_ptr<PositionIterator> pileup(const std::string&);

private:
  unique_htsFile bam_;
  unique_hts_idx_t bam_idx_;
  unique_bam_hdr_t header_;
};

std::unique_ptr<BAMFile> read_bam(const char*);

class PositionIterator {
public:
  friend std::unique_ptr<PositionIterator> BAMFile::pileup(const std::string&);
  std::unique_ptr<Position> next();
  bool has_next();
  // TODO: refactor into iterators
  int begin() { return region_->begin; };
  int end() { return region_->end; };

protected:
  std::unique_ptr<PileupData> pileup_data_;
  std::unique_ptr<bam_mplp_t, decltype(&bam_mplp_destroy)> mplp_iter_;
  std::shared_ptr<const bam_pileup1_t*> pileup_; // sniff... sniff...
  std::unique_ptr<RegionInfo> region_;

  int pos_, tid_, count_, current_next_;
  bool processed_ = true;

  PositionIterator(std::unique_ptr<PileupData>,
                   std::unique_ptr<bam_mplp_t, decltype(&bam_mplp_destroy)>,
                   std::shared_ptr<const bam_pileup1_t*>,
                   std::unique_ptr<RegionInfo>);
};

class Position {
public:
  int position;

  friend std::unique_ptr<Position> PositionIterator::next();
  std::unique_ptr<Alignment> next();
  bool has_next();
  int count() { return count_; };

protected:
  std::string contig_;
  int count_;
  int current_ = 0;
  std::shared_ptr<const bam_pileup1_t*> data_;

  Position(std::string, int, int, std::shared_ptr<const bam_pileup1_t*>);
};

class Alignment {
public:
  friend std::unique_ptr<Alignment> Position::next();
  int is_refskip() { return read_->is_refskip; };
  int is_del() { return read_->is_del; };
  uint32_t query_id() { return read_->b->id; };
  int indel() { return read_->indel; };
  BaseType qbase(const std::int32_t) const noexcept;
  uint8_t qqual(int) const noexcept;
  long ref_start() { return read_->b->core.pos; };
  long ref_end() { return bam_endpos(read_->b); };
  bool rev() { return bam_is_rev(read_->b); };

protected:
  const bam_pileup1_t* read_;

  explicit Alignment(const bam_pileup1_t* read);
};

} // namespace roko

#endif /* ROKO_MODELS_HPP_ */
