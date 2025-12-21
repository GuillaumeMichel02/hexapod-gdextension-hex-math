// hex_transformations.hpp
// This file provides functions for transforming hex coordinates.
#ifndef HEX_TRANSFORMATIONS_HPP
#define HEX_TRANSFORMATIONS_HPP

#include <cstdlib>
#include <cmath>
#include <array>

#include "hex_types.hpp"

namespace hexmath {

    /* ###############################################
        Transform methods
    ############################################### */
    inline Vec2i rotate_once(const Vec2i& coord, int direction = 1) {
        if ( direction > 0 ) {
            return Vec2i( - coord.r, coord.q + coord.r );
        }
        else {
            return Vec2i( coord.q + coord.r, - coord.q );
        }
    }

    inline Vec2i rotate(Vec2i coord, int steps, Vec2i center = Vec2i()) {
        const int mod_steps = steps % 6;
        const int abs_steps = std::abs(mod_steps);
        Vec2i rotated_coords = coord - center;
        for (int i = 0; i < abs_steps; i++) {
            rotated_coords = rotate_once(rotated_coords, mod_steps);
        }
        return rotated_coords + center;
    }

    // Symmetric reflection across an axis. 1 = q, 2 = r, 3 = s.
    // Opposite for symmetric axis: -1 = -q, -2 = -r, -3 = -s.
    inline Vec2i reflect(Vec2i coord, int axis, Vec2i center = Vec2i()) {
        axis = axis % 4; // Normalize axis to -3 to 3

        if (axis == 0) {
            return coord;
        }

        int abs_axis = std::abs(axis);
        Vec2i rc = coord - center;
        std::array<int,3> cube = {rc.q, rc.r, -rc.q - rc.r};
        // Handle permutation depending on axis.
        // Goes in this order : q = (q, s, -q - r), r = (-q - r, r, q).
        rc = Vec2i(cube[(4 - abs_axis) % 3], cube[3 - (abs_axis)]);

        // Apply sign
        if (axis < 0) {
            rc = rc * -1;
        }

        return rc + center;
    }

} // namespace hexmath

#endif // HEX_TRANSFORMATIONS_HPP