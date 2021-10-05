//
// Created by dominik on 10. 10. 2019..
//
#include <memory>
#include <stdexcept>
#include <iostream>
#include <cstring>
#include <cassert>
#include <utility>

#include "htslib/sam.h"
#include "roko/models.hpp"
#include "roko/types.hpp"

namespace roko {

// what's the purpose of this?!?!
// Takes in a PileupData object and the information for an alignment and returns a status number.
// Used in iterating over positions?
int iter_bam(void* data, bam1_t* b) {
  int status;
  PileupData* plp_data = (PileupData*)data;

  // TODO: idiomatic c++ iterators
  while (1) {
    if (plp_data->iter) {
      status = sam_itr_next(plp_data->file, plp_data->iter, b);
    } else {
      status = sam_read1(plp_data->file, plp_data->header, b);
    }
    if (status < 0)
      break; // failed with an error

    if (b->core.flag & FILTER_FLAG)
      continue;
    if (b->core.flag & BAM_FPAIRED && ((b->core.flag & BAM_FPROPER_PAIR) == 0))
      continue;
    if (b->core.qual < MIN_MAPPING_QUALITY)
      continue;
    break;
  }

  return status;
}

RegionInfo::RegionInfo(std::string name, std::int32_t begin, std::int32_t end)
    : name(std::move(name)), begin(begin), end(end) {}

std::unique_ptr<BAMFile> read_bam(const char* filename) {
  auto bam = create_unique_htsFile(filename);
  auto idx = create_unique_hts_idx_t(bam.get(), filename);
  auto header = create_unique_bam_hdr_t(bam.get(), filename);

  return std::make_unique<BAMFile>(std::move(bam), std::move(idx),
                                   std::move(header));
}

BAMFile::BAMFile(unique_htsFile bam, unique_hts_idx_t idx,
                 unique_bam_hdr_t header)
    : bam_(std::move(bam)), bam_idx_(std::move(idx)),
      header_(std::move(header)) {}

std::unique_ptr<RegionInfo> parse_region(const std::string& region) {
  std::int32_t start, end;

  // 13 145
  const auto end_name = hts_parse_reg(region.c_str(), &start, &end);
  const auto len = end_name - region.c_str(); // 132

  // TODO: sniff sniff
  // 13 146
  std::string contig(region, 0, len);

  return std::unique_ptr<RegionInfo>(
      new RegionInfo(std::move(contig), start, end));
}

std::unique_ptr<PositionIterator> BAMFile::pileup(const std::string& region) {
  std::unique_ptr<PileupData> data(new PileupData);
   
  //Filling up the PileupData struct with content 
  data->file = this->bam_.get();
  data->header = this->header_.get();
  data->iter =
      bam_itr_querys(this->bam_idx_.get(), this->header_.get(), region.c_str());

  // Creating multi-iterator
  auto data_raw = data.get(); // The wrapped pointer to the PileupData struct
  bam_mplp_t mplp = bam_mplp_init(1, iter_bam, (void**)&data_raw);
  std::unique_ptr<bam_mplp_s, decltype(&bam_mplp_destroy)> mplp_iter(
      mplp, bam_mplp_destroy);

  // Pointer to data array for one position
  std::shared_ptr<const bam_pileup1_t*> pileup(
      const_cast<const bam_pileup1_t**>(new bam_pileup1_t*));

  // Region info
  auto region_info = parse_region(region);

  return std::unique_ptr<PositionIterator>(
      new PositionIterator(std::move(data), std::move(mplp_iter),
                           std::move(pileup), std::move(region_info)));
}

PositionIterator::PositionIterator(
    std::unique_ptr<PileupData> pileup_data,
    std::unique_ptr<bam_mplp_t, decltype(&bam_mplp_destroy)> mplp_iter,
    std::shared_ptr<const bam_pileup1_t*> pileup,
    std::unique_ptr<RegionInfo> region)
    : pileup_data_(std::move(pileup_data)), mplp_iter_(std::move(mplp_iter)),
      pileup_(std::move(pileup)), region_(std::move(region)) {}

bool PositionIterator::has_next() {
  if (processed_) {
    current_next_ =
        bam_mplp_auto(mplp_iter_.get(), &tid_, &pos_, &count_, pileup_.get());
    processed_ = false;
  }

  return current_next_ > 0;
}

std::unique_ptr<Position> PositionIterator::next() {
  if (!has_next()) {
    throw std::runtime_error("No more positions to iterate.");
  }

  const char* contig_name = pileup_data_->header->target_name[tid_];
  assert(region_->name == contig_name);
  assert(pos_ >= region_->begin);
  assert(pos_ < region_->end);

  processed_ = true;
  return std::unique_ptr<Position>(
      new Position(region_->name, pos_, count_, pileup_));
}

Position::Position(std::string contig, int pos, int count,
                   std::shared_ptr<const bam_pileup1_t*> data)
    : position(pos), contig_(std::move(contig)), count_(count),
      data_(std::move(data)) {}

bool Position::has_next() { return current_ < count_; }

std::unique_ptr<Alignment> Position::next() {
  if (current_ >= count_) {
    throw std::runtime_error("No more reads in position to iterate.");
  }

  const bam_pileup1_t* read = *data_ + current_;
  current_++;

  return std::unique_ptr<Alignment>(new Alignment(read));
}

Alignment::Alignment(const bam_pileup1_t* read) : read_(read) {}

BaseType Alignment::qbase(const std::int32_t offset) const noexcept {
  const auto seq = bam_get_seq(read_->b);
  const auto base = bam_seqi(seq, read_->qpos + offset);

  switch (base) {
  case 1:
    return BaseType::A;
    break;
  case 2:
    return BaseType::C;
    break;
  case 4:
    return BaseType::G;
    break;
  case 8:
    return BaseType::T;
    break;
  }

  return BaseType::UNKNOWN;
}

std::uint8_t Alignment::qqual(int offset) const noexcept {
  const auto qual_string = bam_get_qual(read_->b);
  if (qual_string) {
    return qual_string[read_->qpos + offset];
  }

  return 10;
}
BaseType get_base(const char base_raw) noexcept {
  const char base_low = std::tolower(base_raw);
  // TODO: ditch switch case
  switch (base_raw) {
  case 'a':
    return BaseType::A;
  case 'c':
    return BaseType::C;
  case 'g':
    return BaseType::G;
  case 't':
    return BaseType::T;
  case 'n':
  case '-':
    return BaseType::UNKNOWN;
  case '*':
    return BaseType::GAP;
  }
}

} // namespace roko