#include "scan.h"
#include "motiv.h"
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace std;

PreparedMotifs prepare_motifs(const vector<string>& raw_motifs, int motif_len) {
    PreparedMotifs out;
    for (const auto& raw : raw_motifs) {
        if (motif_len <= 0) continue;
        if (static_cast<int>(raw.size()) < motif_len) continue;
        string m = raw.substr(0, motif_len);
        out.motifs.push_back(m);
        out.masks.push_back(compile_motif_masks(m));
    }
    return out;
}

static int count_matches_serial(
    const vector<string>& sequences,
    const vector<uint8_t>& masks
) {
    int count = 0;
    for (const auto& s : sequences)
        if (contains_motif_compiled(s, masks)) ++count;
    return count;
}

static int count_matches_omp(
    const vector<string>& sequences,
    const vector<uint8_t>& masks,
    int threads
) {
#ifndef _OPENMP
    (void)threads;
    return count_matches_serial(sequences, masks);
#else
    int count = 0;
    const int N = static_cast<int>(sequences.size());
    omp_set_num_threads(threads);
#pragma omp parallel for reduction(+:count) schedule(static)
    for (int si = 0; si < N; ++si)
        if (contains_motif_compiled(sequences[si], masks)) ++count;
    return count;
#endif
}

static int count_for_motif(
    const vector<string>& sequences,
    const vector<uint8_t>& masks,
    int threads
) {
    return threads > 1
        ? count_matches_omp(sequences, masks, threads)
        : count_matches_serial(sequences, masks);
}

ScanResult scan_motifs(
    const vector<string>& sequences,
    const PreparedMotifs& motifs,
    int threads,
    int mpi_rank,
    int mpi_size
) {
    ScanResult result;
    const int M = static_cast<int>(motifs.motifs.size());
    if (M == 0) return result;

#ifdef USE_MPI
    if (mpi_size > 1) {
        const int base = M / mpi_size;
        const int rem  = M % mpi_size;
        const int local_n     = base + ((mpi_rank < rem) ? 1 : 0);
        const int local_start = mpi_rank * base + ((mpi_rank < rem) ? mpi_rank : rem);

        vector<int> local_counts(local_n, 0);
        for (int k = 0; k < local_n; ++k) {
            const int mi = local_start + k;
            local_counts[k] = count_for_motif(sequences, motifs.masks[mi], threads);
        }

        vector<int> recv_counts;
        vector<int> recvcounts, displs;
        if (mpi_rank == 0) {
            recvcounts.resize(mpi_size);
            displs.resize(mpi_size);
            int disp = 0;
            for (int r = 0; r < mpi_size; ++r) {
                const int rn = base + ((r < rem) ? 1 : 0);
                recvcounts[r] = rn;
                displs[r] = disp;
                disp += rn;
            }
            recv_counts.resize(M);
        }

        MPI_Gatherv(
            local_counts.data(), local_n, MPI_INT,
            mpi_rank == 0 ? recv_counts.data() : nullptr,
            mpi_rank == 0 ? recvcounts.data() : nullptr,
            mpi_rank == 0 ? displs.data() : nullptr,
            MPI_INT, 0, MPI_COMM_WORLD
        );

        if (mpi_rank == 0)
            result.counts = std::move(recv_counts);
        return result;
    }
#else
    (void)mpi_rank;
    (void)mpi_size;
#endif

    result.counts.resize(M, 0);
    for (int mi = 0; mi < M; ++mi)
        result.counts[mi] = count_for_motif(sequences, motifs.masks[mi], threads);
    return result;
}

void print_results_tsv(
    const PreparedMotifs& motifs,
    const ScanResult& result,
    int total_sequences,
    ostream& out
) {
    out << "motif\tmatches\tsequences\tfrequency\n";
    const int M = static_cast<int>(motifs.motifs.size());
    for (int i = 0; i < M; ++i) {
        const int count = (i < static_cast<int>(result.counts.size())) ? result.counts[i] : 0;
        const double freq = total_sequences > 0
            ? static_cast<double>(count) / total_sequences
            : 0.0;
        out << motifs.motifs[i] << '\t' << count << '\t'
            << total_sequences << '\t' << freq << '\n';
    }
}

void print_usage(ostream& out) {
    out <<
        "Motif scanner\n\n"
        "Usage:\n"
        "  ./main --fasta <file> --motifs <file> --threads <N> [options]\n"
        "  mpirun -np <N> ./main --fasta <file> --motifs <file> --threads <N>\n\n"
        "Options:\n"
        "  --fasta <file>       FASTA input (required)\n"
        "  --motifs <file>      motif list input (required)\n"
        "  --threads <N>        OpenMP threads per process (required)\n"
        "  --motif-len <N>      motif length, default 8\n"
        "  -h, --help           show this help\n\n"
        "Output (TSV to stdout, rank 0 only):\n"
        "  motif  matches  sequences  frequency\n\n"
        "Examples:\n"
        "  ./main --fasta data.fst --motifs motifs.mot --threads 4\n"
        "  mpirun -np 4 ./main --fasta data.fst --motifs motifs.mot --threads 2\n";
}
