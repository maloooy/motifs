#include <iostream>
#include <string>
#include <cstdlib>
#ifdef USE_MPI
#include <mpi.h>
#endif
#include "io.h"
#include "scan.h"

using namespace std;

static int fail(const string& message, int code = 2) {
    cerr << message << '\n';
    return code;
}

int main(int argc, char** argv) {
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
    int mpi_rank = 0, mpi_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#else
    const int mpi_rank = 0, mpi_size = 1;
#endif

    string fasta_file;
    string motif_file;
    int motif_len = 8;
    int threads = 0;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        auto need_value = [&](const char* name) -> string {
            if (i + 1 >= argc) throw runtime_error(string("missing value for ") + name);
            return argv[++i];
        };

        try {
            if (arg == "--fasta") fasta_file = need_value("--fasta");
            else if (arg == "--motifs") motif_file = need_value("--motifs");
            else if (arg == "--motif-len") motif_len = stoi(need_value("--motif-len"));
            else if (arg == "--threads") threads = stoi(need_value("--threads"));
            else if (arg == "--help" || arg == "-h") {
                if (mpi_rank == 0) print_usage(cout);
#ifdef USE_MPI
                MPI_Finalize();
#endif
                return 0;
            } else {
                if (mpi_rank == 0) return fail("unknown argument: " + arg);
#ifdef USE_MPI
                MPI_Abort(MPI_COMM_WORLD, 2);
#endif
                return 2;
            }
        } catch (const exception& e) {
            if (mpi_rank == 0) return fail(string("invalid value for ") + arg + ": " + e.what());
#ifdef USE_MPI
            MPI_Abort(MPI_COMM_WORLD, 2);
#endif
            return 2;
        }
    }

    auto validate = [&](const string& message, int code = 2) -> int {
        if (mpi_rank == 0) return fail(message, code);
#ifdef USE_MPI
        MPI_Abort(MPI_COMM_WORLD, code);
#endif
        return code;
    };

    if (fasta_file.empty() || motif_file.empty())
        return validate("both --fasta and --motifs are required\n\nTry: ./main --help");
    if (threads <= 0)
        return validate("--threads is required and must be positive");
    if (motif_len <= 0)
        return validate("--motif-len must be positive");

    auto sequences = read_fasta(fasta_file);
    auto motifs_raw = read_motifs(motif_file);
    if (sequences.empty())
        return validate("no sequences found in " + fasta_file, 1);
    if (motifs_raw.empty())
        return validate("no motifs found in " + motif_file, 1);

    PreparedMotifs motifs;
    try {
        motifs = prepare_motifs(motifs_raw, motif_len);
    } catch (const exception& e) {
        return validate(e.what());
    }
    if (motifs.motifs.empty())
        return validate("no motifs are at least --motif-len characters long", 1);

    const int total = static_cast<int>(sequences.size());
    ScanResult result = scan_motifs(sequences, motifs, threads, mpi_rank, mpi_size);

    if (mpi_rank == 0)
        print_results_tsv(motifs, result, total, cout);

#ifdef USE_MPI
    MPI_Finalize();
#endif
    return 0;
}
