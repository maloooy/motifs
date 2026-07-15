#include "motiv.h"
#include <cctype>

static inline uint8_t dna_mask(char dna) {
    switch (dna) {
        case 'A': return 0x1;
        case 'C': return 0x2;
        case 'G': return 0x4;
        case 'T': return 0x8;
        default: return 0x0;
    }
}

static inline uint8_t iupac_mask(char motif) {
    switch (motif) {
        case 'A': return 0x1;
        case 'C': return 0x2;
        case 'G': return 0x4;
        case 'T': return 0x8;
        case 'R': return 0x1 | 0x4;
        case 'Y': return 0x2 | 0x8;
        case 'M': return 0x1 | 0x2;
        case 'K': return 0x4 | 0x8;
        case 'W': return 0x1 | 0x8;
        case 'S': return 0x2 | 0x4;
        case 'B': return 0x2 | 0x4 | 0x8;
        case 'V': return 0x1 | 0x2 | 0x4;
        case 'H': return 0x1 | 0x2 | 0x8;
        case 'D': return 0x1 | 0x4 | 0x8;
        case 'N': return 0xF;
        default: return 0x0;
    }
}

std::vector<uint8_t> compile_motif_masks(const std::string& motif) {
    std::vector<uint8_t> masks;
    masks.reserve(motif.size());
    for (unsigned char ch : motif) {
        masks.push_back(iupac_mask(static_cast<char>(std::toupper(ch))));
    }
    return masks;
}

bool contains_motif_compiled(const std::string& s, const std::vector<uint8_t>& motif_masks) {
    const int n = static_cast<int>(s.size());
    const int m = static_cast<int>(motif_masks.size());
    if (m == 0 || n < m) return false;

    for (int i = 0; i <= n - m; ++i) {
        bool ok = true;
        for (int j = 0; j < m; ++j) {
            const uint8_t dm = dna_mask(s[i + j]);
            if ((dm & motif_masks[j]) == 0) {
                ok = false;
                break;
            }
        }
        if (ok) return true;
    }
    return false;
}
