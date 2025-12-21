#ifndef HEX_MATH_CORE_H
#define HEX_MATH_CORE_H

#include <string>

#include "hex_math/hex_constants.hpp"
#include "hex_math/hex_types.hpp"
#include "hex_math/hex_conversions.hpp"
#include "hex_math/hex_directions.hpp"
#include "hex_math/hex_relations.hpp"
#include "hex_math/hex_dda.hpp"
#include "hex_math/hex_neighbors.hpp"
#include "hex_math/hex_ranges.hpp"
#include "hex_math/hex_transformations.hpp"

namespace hexmath {

    // Helper : strings out whatever I want to see in the editor
    inline std::string helper() {
        std::string res = "GDE Hex Math Core Module\n";
        return res;

    }
} // namespace hexmath

#endif // HEX_MATH_CORE_H