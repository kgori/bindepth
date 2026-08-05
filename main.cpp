#include "CLI11.hpp"
#include "htslib/sam.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

struct Config {
  std::string bam;
  int window_size = 1000;
  int min_mapq = 20;
  int n_threads = 1;
};

Config parse_args(int argc, char **argv) {
  Config config;
  CLI::App app{"My CLI App"};

  app.add_option("-b,--bam", config.bam, "Path to BAM file")->required();
  app.add_option("-w,--window", config.window_size,
                 "Window size for depth calculation (default: 1000)");
  app.add_option("-q,--min-mapq", config.min_mapq,
                 "Minimum mapping quality for reads (default: 20)");
  app.add_option("-t,--threads", config.n_threads,
                 "Number of threads to use (default: 1)");
  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    std::exit(app.exit(e));
  }

  return config;
}

constexpr int OP_TYPE_CONSUMES_REF = 2;

int main(int argc, char **argv) {
  // Decouple C++ streams from C streams for performance
  std::ios_base::sync_with_stdio(false);

  Config config = parse_args(argc, argv);

  samFile *fp = sam_open(config.bam.c_str(), "r");

  if (fp == nullptr) {
    std::cerr << "Error opening BAM file: " << config.bam << std::endl;
    return 1;
  }

  sam_hdr_t *hdr = sam_hdr_read(fp);
  if (hdr == nullptr) {
    std::cerr << "Error reading BAM header: " << config.bam << std::endl;
    sam_close(fp);
    return 1;
  }

  int window_size = config.window_size;
  int min_mapq = config.min_mapq;
  unsigned int threads =
      std::min(std::thread::hardware_concurrency(),
               static_cast<unsigned int>(std::max(1, config.n_threads)));
  hts_set_threads(fp, threads);
  int nref = sam_hdr_nref(hdr);

  int cur_tid = -1;
  std::vector<uint32_t> bases;
  bam1_t *b = bam_init1();

  auto flush_contig = [&](int tid) {
    if (tid < 0 || tid >= nref) {
      return;
    }
    int len = sam_hdr_tid2len(hdr, tid);
    int n_windows = (len + window_size - 1) / window_size;
    const char *chrom = sam_hdr_tid2name(hdr, tid);

    bool has_data = (tid == cur_tid) && !bases.empty();

    for (int w = 0; w < n_windows; ++w) {
      int start = w * window_size;
      int end = std::min((w + 1) * window_size, len);
      int width = end - start;

      uint32_t covered = has_data ? bases[w] : 0;
      double depth = static_cast<double>(covered) / width;
      std::cout << chrom << "\t" << start << "\t" << end << "\t" << depth
                << '\n';
    }
  };

  long long reads_seen = 0;
  const long long PROGRESS_INTERVAL = 1000000;

  std::cout << "# Alignment file: " << config.bam << '\n';
  std::cout << "# Window size: " << window_size
            << " (final bin in contigs may be truncated)\n";
  std::cout << "# Minimum mapping quality: " << min_mapq << '\n';
  std::cout << "# DEPTH=Number of mapped bases within a window, divided by "
               "window width (adjustment made for narrower bins at ends of "
               "contigs).\n";
  std::cout << "CHROM\tSTART\tEND\tDEPTH\n";
  while (sam_read1(fp, hdr, b) >= 0) {
    ++reads_seen;
    if (reads_seen % PROGRESS_INTERVAL == 0) {
      std::cerr << "processed " << reads_seen << " reads. Current contig = "
                << sam_hdr_tid2name(hdr, b->core.tid) << '\r' << std::flush;
    }
    if (b->core.flag &
        (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FDUP)) {
      continue;
    }
    if (b->core.qual < min_mapq) {
      continue;
    }
    if (b->core.tid < 0 || b->core.tid >= nref) {
      continue;
    }

    if (b->core.tid != cur_tid) {
      for (int cur = cur_tid; cur < b->core.tid; ++cur) {
        flush_contig(cur);
      }
      cur_tid = b->core.tid;
      int len = sam_hdr_tid2len(hdr, cur_tid);
      int n_windows = (len + window_size - 1) / window_size;
      bases.assign(n_windows, 0);
    }

    int ref = b->core.pos;
    uint32_t *cigar = bam_get_cigar(b);
    for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
      int op = bam_cigar_op(cigar[i]);
      int len = bam_cigar_oplen(cigar[i]);

      if (bam_cigar_type(op) &
          OP_TYPE_CONSUMES_REF) { // 2 means consumes reference
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
          int block_start = ref;
          int block_end = ref + len;

          int first = std::max(0, block_start / window_size);
          int last =
              std::min((int)bases.size() - 1, (block_end - 1) / window_size);

          for (int w = first; w <= last; ++w) {
            int win_start = w * window_size;
            int win_end = win_start + window_size;
            int overlap =
                std::min(block_end, win_end) - std::max(block_start, win_start);

            if (overlap > 0)
              bases[w] += overlap;
          }
        }
        ref += len;
      }
    }
  }
  std::cerr << '\n';
  for (int cur = cur_tid; cur < nref; ++cur) {
    flush_contig(cur);
  }

  bam_destroy1(b);
  sam_hdr_destroy(hdr);
  sam_close(fp);
  return 0;
}
