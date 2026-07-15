#pragma once
#include <string>
#include <vector>
#include <cstdint>

std::vector<uint8_t> compile_motif_masks(const std::string& motif);
bool contains_motif_compiled(const std::string& s, const std::vector<uint8_t>& motif_masks);
