// hex_directions.hpp
// This file defines the six primary directions and corner offsets for hexagonal grids.
// Includes:
// - `DIRECTIONS` array for axial directions
// - `CORNERS` array for corner offsets
#ifndef HEX_DIRECTIONS_HPP
#define HEX_DIRECTIONS_HPP

#include "hex_types.hpp"

namespace hexmath {

    inline const Vec2i DIRECTIONS[6] = {Vec2i(1,0), Vec2i(1,-1), Vec2i(0,-1), Vec2i(-1,0), Vec2i(-1,1),Vec2i(0,1)};

    inline const Vec2i CORNER_DIRECTIONS[6] = {
        Vec2i(DIRECTIONS[0] + DIRECTIONS[1]),
        Vec2i(DIRECTIONS[1] + DIRECTIONS[2]),
        Vec2i(DIRECTIONS[2] + DIRECTIONS[3]),
        Vec2i(DIRECTIONS[3] + DIRECTIONS[4]),
        Vec2i(DIRECTIONS[4] + DIRECTIONS[5]),
        Vec2i(DIRECTIONS[5] + DIRECTIONS[0])
    };

    inline const Vec2f CORNERS[6] = {
        Vec2f(CORNER_DIRECTIONS[0]) * (1.0f/3),
        Vec2f(CORNER_DIRECTIONS[1]) * (1.0f/3),
        Vec2f(CORNER_DIRECTIONS[2]) * (1.0f/3),
        Vec2f(CORNER_DIRECTIONS[3]) * (1.0f/3),
        Vec2f(CORNER_DIRECTIONS[4]) * (1.0f/3),
        Vec2f(CORNER_DIRECTIONS[5]) * (1.0f/3)
    };

} // namespace hexmath

#endif // HEX_DIRECTIONS_HPP