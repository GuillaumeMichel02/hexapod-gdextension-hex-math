// hex_conversions.hpp
// This file provides functions for converting between different coordinate systems.

#ifndef HEX_CONVERSIONS_HPP
#define HEX_CONVERSIONS_HPP

#include <cstdlib>
#include <cmath>

#include "hex_constants.hpp"
#include "hex_types.hpp"

namespace hexmath {

    /* ###############################################
        Coordinate conversion methods
    ############################################### */

    inline Vec2f hex_to_pixel(Vec2i coord, float size) {
        float x = size * SQRT_3 * (coord.q + coord.r / 2.0f);
        float y = size * 1.5f * coord.r;
        return Vec2f(x, y);
    }

    inline Vec2i round_hex(float q, float r) {
        float s = -q - r;

        int round_q = std::round(q);
        int round_r = std::round(r);
        int round_s = std::round(s);

        float q_diff = std::fabs(q - round_q);
        float r_diff = std::fabs(r - round_r);
        float s_diff = std::fabs(s - round_s);

        if (q_diff > r_diff && q_diff > s_diff) {
            round_q = -round_r - round_s;
        }
        else if (r_diff > s_diff) {
            round_r = -round_q - round_s;
        }

        return Vec2i(round_q, round_r);
    }

    inline Vec2i pixel_to_hex(Vec2f pixel, float size) {
        float q = (SQRT_3/3 * pixel.q - 1.0f/3 * pixel.r) / size;
        float r = (2.0f/3 * pixel.q) / size;
        return round_hex(q, r);
    }

} // namespace hexmath

#endif // HEX_CONVERSIONS_HPP