// hex_ranges.hpp
// This file contains functions for range-based queries on hexagonal grids.
#ifndef HEX_RANGES_HPP
#define HEX_RANGES_HPP

#include <vector>
#include <cmath>

#include "hex_constants.hpp"
#include "hex_types.hpp"
#include "hex_relations.hpp"
#include "hex_directions.hpp"
#include "hex_dda.hpp"

namespace hexmath {

    /* ###############################################
        Range methods
    ############################################### */

    inline bool is_in_range(Vec2i coord, int radius, Vec2i center = Vec2i()) {
        Vec2i rel_coord = coord - center;
        return rel_coord.length() <= radius;
    }

    inline bool is_in_cone(Vec2i coord, int direction, int angle_width, Vec2i center = Vec2i()) {
        Vec2i rel_coord = coord - center;
        float angle_to_coord = get_angle(DIRECTIONS[direction % 6], rel_coord);
        return std::fabs(angle_to_coord) <= (angle_width / 2.0f) * (PI / 180.0f);
    }

    inline bool is_in_ring(Vec2i coord, int radius, Vec2i center = Vec2i()) {
        Vec2i rel_coord = coord - center;
        return rel_coord.length() == radius;
    }

    inline bool is_in_supercover(Vec2i coord, Vec2i from, Vec2i to = Vec2i()) {
        return distance_to_line_border(coord, from, to) < EPSILON;
    }

    // inline bool is_crossed_by_line(Vec2i coord, Vec2i from, Vec2i to = Vec2i()) {
    //     Vec2i direction = to - from;
    //     Vec2i rel_coord = coord - from;

    // }

    inline std::vector<Vec2i> get_all_in_range(int radius, Vec2i center = Vec2i()) {
        std::vector<Vec2i> in_radius;
        for (int q = -radius; q <= radius; q++) {
            for (int r = -radius; r <= radius; r++) {
                if (std::abs(q + r) > radius) continue;
                in_radius.push_back(Vec2i(q, r) + center);
            }
        }
        return in_radius;
    }

    // width : Size of the ring segment in number of directions (6 = full ring)
    inline std::vector<Vec2i> get_all_in_ring(int radius, int width = 6, int init_direction = 0, Vec2i center = Vec2i()) {
        std::vector<Vec2i> in_ring;

        if (radius == 0) {
            in_ring.push_back(center);
            return in_ring;
        }

        for (int i = 0; i < width; i++) {
                int dir_index = (i + init_direction) % 6;
                Vec2i corner = center + DIRECTIONS[dir_index] * radius;
                Vec2i direction = DIRECTIONS[(dir_index+2)%6];
            for (int j = 0; j < radius; j++) {
                in_ring.push_back(corner + direction * j);
            }
        }
        return in_ring;
    }

    inline std::vector<Vec2i> get_all_in_spiral(int radius, int width = 6, int init_direction = 0, Vec2i center = Vec2i()) {
        std::vector<Vec2i> in_spiral;
        for (int i = 0; i <= radius; i++) {
            for (const auto& hex : get_all_in_ring(i, width, init_direction, center)) {
                in_spiral.push_back(hex);
            }
        }
        return in_spiral;
    }

    inline std::vector<Vec2i> get_all_in_line(Vec2i from, Vec2i to = Vec2i()) {
        std::vector<Vec2i> in_line;

        int N = distance(from, to);
        Vec2f current_interpolation;
        for (int i = 0; i < N; i++) {
            // Add epsilon to avoid edge cases where the lerp is exactly on a hex boundary.
            current_interpolation = from.lerp(to, i * (1.0f / N) + EPSILON);
            in_line.push_back(current_interpolation.round());
        }

        // Pushing back "to" out of the loop fixes the edge case where N = 0.
        in_line.push_back(to);

        return in_line;
    }

    inline std::vector<Vec2i> get_all_in_supercover(Vec2i from, Vec2i to = Vec2i()) {
        return hex_dda(from, to);
    }

} // namespace hexmath

#endif // HEX_RANGES_HPP