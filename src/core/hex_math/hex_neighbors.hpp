// hex_neighbors.hpp
// This file provides functions for finding neighboring hexes in a grid.
#ifndef HEX_NEIGHBORS_HPP
#define HEX_NEIGHBORS_HPP

#include <array>

#include "hex_types.hpp"
#include "hex_directions.hpp"

namespace hexmath {

    /* ###############################################
        Neighbor methods
    ############################################### */
    
    inline Vec2i get_neighbor(Vec2i coord, int direction) { 
        return coord + DIRECTIONS[direction % 6];
    }

    inline std::array<Vec2i, 6> get_neighbors(Vec2i coord) {
        std::array<Vec2i, 6> neighbors;
        for (int i = 0 ; i < 6 ; i++) {
            neighbors[i] = coord + DIRECTIONS[i];
        }
        return neighbors;
    }

} // namespace hexmath

#endif // HEX_NEIGHBORS_HPP