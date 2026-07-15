#pragma once
#include <string>
#include <vector>
#include <cstdint>
#include <ostream>

struct PreparedMotifs {
    std::vector<std::string>          motifs;
    std::vector<std::vector<uint8_t>> masks;
};

struct ScanResult {
    std::vector<int> counts;
};

PreparedMotifs prepare_motifs(const std::vector<std::string>& raw_motifs, int motif_len);

ScanResult scan_motifs(
    const std::vector<std::string>& sequences,
    const PreparedMotifs& motifs,
    int threads,
    int mpi_rank,
    int mpi_size
);

void print_results_tsv(
    const PreparedMotifs& motifs,
    const ScanResult& result,
    int total_sequences,
    std::ostream& out
);

void print_usage(std::ostream& out);
