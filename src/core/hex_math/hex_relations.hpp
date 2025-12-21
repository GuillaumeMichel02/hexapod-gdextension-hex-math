// hex_relations.hpp
// This file provides functions for calculating relationships between hexes.
#ifndef HEX_RELATIONS_HPP
#define HEX_RELATIONS_HPP

#include "hex_types.hpp"
#include "hex_directions.hpp"

namespace hexmath {

    /* ###############################################
        Relation methods
    ############################################### */

    inline int distance(Vec2i a, Vec2i b) {
        return (a - b).length();
    }

    inline float distance_to_line(Vec2i coord, Vec2i from, Vec2i to) {
        Vec2f projection = (coord - from).project_onto(to - from);
        Vec2f perp = Vec2f(coord.q - from.q, coord.r - from.r) - projection;
        return perp.length();
    }

    inline float distance_to_line_border(Vec2i coord, Vec2i from, Vec2i to) {
        Vec2f border;
        Vec2f projection;
        Vec2f perp;
        float min_distance = std::numeric_limits<float>::max();
        for (int i = 0; i < 6; i++) {
            border = Vec2f(coord) + CORNERS[i];
            projection = (border - Vec2f(from)).project_onto(Vec2f(to) - Vec2f(from));
            perp = Vec2f(border.q - from.q, border.r - from.r) - projection;
            float dist = perp.length();
            if (dist < min_distance) {
                min_distance = dist;
            }
        }
        return min_distance;
    }

    inline Vec2f project(Vec2i coord, Vec2i line) {
        return coord.project_onto(line);
    }

    inline int cross_product(const Vec2i& a, const Vec2i& b) {
        return 3* (a.q * b.r - a.r * b.q);
    }

    inline float dot_product(const Vec2i& a, const Vec2i& b) {
        // In axial coordinates with 120° angle between axes:
        return a.q * b.q + a.r * b.r + 0.5f * (a.q * b.r + a.r * b.q);
    }

    inline int cross_product(const Vec2f& a, const Vec2f& b) {
        return 3* (a.q * b.r - a.r * b.q);
    }

    inline float dot_product(const Vec2f& a, const Vec2f& b) {
        // In axial coordinates with 120° angle between axes:
        return a.q * b.q + a.r * b.r + 0.5f * (a.q * b.r + a.r * b.q);
    }

    inline float get_angle(const Vec2i& a, const Vec2i& b) {
        float dot = dot_product(a, b);
        float len_a = a.length();
        float len_b = b.length();
        if (len_a == 0 || len_b == 0) return 0.0f;
        float cos_theta = dot / (len_a * len_b);
        // Clamp to avoid numerical issues
        if (cos_theta > 1.0f) cos_theta = 1.0f;
        if (cos_theta < -1.0f) cos_theta = -1.0f;
        return std::acos(cos_theta); // Return angle in degrees
    }

} // namespace hexmath

#endif // HEX_MATH_HPP