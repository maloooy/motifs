#include "io.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cctype>

using namespace std;

static inline void rtrim_inplace(string& s) {
    while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '\n' || s.back() == '\r'))
        s.pop_back();
}

static inline bool is_all_whitespace(const string& s) {
    for (unsigned char c : s) {
        if (!std::isspace(c)) return false;
    }
    return true;
}

vector<string> read_fasta(const string& filename) {
    ifstream file(filename);
    vector<string> sequences;

    if (!file) {
        cerr << "Error opening FASTA file: " << filename << "\n";
        return sequences;
    }

    string line, current;

    while (getline(file, line)) {
        if (line.empty() || is_all_whitespace(line)) continue;

        if (line[0] == '>') {
            if (!current.empty()) {
                sequences.push_back(current);
                current.clear();
            }
        } else {
            rtrim_inplace(line);
            if (!line.empty()) {
                for (unsigned char ch : line) {
                    current.push_back(static_cast<char>(std::toupper(ch)));
                }
            }
        }
    }

    if (!current.empty())
        sequences.push_back(current);

    return sequences;
}


vector<string> read_motifs(const string& filename) {
    ifstream file(filename);
    vector<string> motifs;

    if (!file) {
        cerr << "Error opening MOT file: " << filename << "\n";
        return motifs;
    }

    string line;
    while (getline(file, line)) {
        if (line.empty() || is_all_whitespace(line)) continue;
        std::istringstream iss(line);
        std::string motif;
        if (iss >> motif) {
            for (char& ch : motif) {
                ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
            }
            motifs.push_back(motif);
        }
    }

    return motifs;
}
