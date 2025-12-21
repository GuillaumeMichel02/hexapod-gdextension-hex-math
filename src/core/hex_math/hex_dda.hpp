// hex_dda.hpp
// This file implements the Digital Differential Analyzer (DDA) algorithm for hexagonal grids.
// Includes functions for:
// - Tracing a path between two hexes
// - Finding boundary crossings
// - Generating intermediate hexes along a line

#ifndef HEX_DDA_HPP
#define HEX_DDA_HPP

#include <cstdlib>
#include <cmath>
#include <array>
#include <vector>
#include <set>

#include "hex_constants.hpp"
#include "hex_types.hpp"
#include "hex_relations.hpp"
#include "hex_conversions.hpp"

namespace hexmath {

    inline bool is_straight_line(const Vec3f& direction) {
        // One component is 0
        int zero_count = 0;
        if (std::fabs(direction.q) < EPSILON) zero_count++;
        if (std::fabs(direction.r) < EPSILON) zero_count++;
        if (std::fabs(direction.s) < EPSILON) zero_count++;
        return zero_count == 1;
    }

    inline bool is_corner_aligned(const Vec3f& direction) {
        // Two components have equal absolute value

        return (std::fabs(direction.q - direction.r) < EPSILON) ||
            (std::fabs(direction.q - direction.s) < EPSILON) ||
            (std::fabs(direction.r - direction.s) < EPSILON);
    }

    inline std::set<float> get_crossing_intervals(float start, float delta) {
        std::set<float> crossing_t_values;

        if (delta == 0) {
            return crossing_t_values; // No crossings if no movement
        }
        
        float min_val = std::min(start, start + delta);
        float max_val = std::max(start, start + delta);

        for (int n = std::ceil(min_val - 0.5f); n <= std::floor(max_val - 0.5f); n++) {
            float t = (n + 0.5f - start) / delta;
            if (t <= 0.0f || t > 1.0f) continue;
            crossing_t_values.insert(t);
        }

        return crossing_t_values;

    }

    inline std::set<float> solve_abs_eq(float start_a, float delta_a, float start_b, float delta_b, float start_c, float delta_c) {
        // Solve |start_a + delta_a * t| = |start_b + delta_b * t|
        // This gives two linear equations to solve:
        // 1. start_a + delta_a * t = start_b + delta_b * t
        // 2. start_a + delta_a * t = - (start_b + delta_b * t)

        std::set<float> solutions;

        if (std::fabs(delta_a - delta_b) > EPSILON) {
            float t1 = (start_b - start_a) / (delta_a - delta_b);
            if (t1 >= 0.0f && t1 <= 1.0f) {
            solutions.insert(t1);
            }
        }

        if (std::fabs(delta_a + delta_b) > EPSILON) {
            float t2 = (-start_b - start_a) / (delta_a + delta_b);
            if (t2 >= 0.0f && t2 <= 1.0f) {
            solutions.insert(t2);
            }
        }

        // Check if solutions are valid. A solution is valid if the third equation result is less than the other two.
        std::set<float> valid_solutions;
        for (float t : solutions) {
            float solution_value = std::fabs(start_a + delta_a * t);
            float third_value = std::fabs(start_c + delta_c * t);
            if (third_value <= solution_value + EPSILON) {
                valid_solutions.insert(t);
            }
        }
        
        return valid_solutions;
    }

    inline float find_first_boundary_crossing(const Vec3f& start, const Vec3f& end) {
        Vec3f delta = end - start;

        // Step 1: Find all t values where a coordinate crosses n + 0.5
        std::set<float> crossing_t_values;
        crossing_t_values.insert(0.0f); // Start at t=0
        for (float t : get_crossing_intervals(start.q, delta.q)) {
            crossing_t_values.insert(t);
        }
        for (float t : get_crossing_intervals(start.r, delta.r)) {
            crossing_t_values.insert(t);
        }
        for (float t : get_crossing_intervals(start.s, delta.s)) {
            crossing_t_values.insert(t);
        }
        crossing_t_values.insert(1.0f); // End at t=1
    
        // Step 2: For each interval, find boundary crossings
        // The first crossing will be the first hex boundary crossed
        auto it = crossing_t_values.begin();
        while (it != crossing_t_values.end()) {
            float t_start = *it;
            ++it;
            if (it == crossing_t_values.end()) break;
            float t_end = *it;

            float t_middle = (t_start + t_end) / 2.0f;

            Vec3f middle(
                start.q + delta.q * t_middle,
                start.r + delta.r * t_middle,
                start.s + delta.s * t_middle
            );

            // Compute rounded values, constant within this interval
            Vec3f rounded(
                std::round(middle.q),
                std::round(middle.r),
                std::round(middle.s)
            );

            // Compute round differences
            Vec3f round_diff(
                start.q - rounded.q,
                start.r - rounded.r,
                start.s - rounded.s
            );

            if (round_diff == Vec3f(0,0,0)) {
                continue; // It's the middle of a hex, no crossing here
            }

            // Solve for boundary crossings in this interval
            std::set<float> t_candidates;
            std::set<float> t_qr = solve_abs_eq(
                round_diff.q, delta.q,
                round_diff.r, delta.r,
                round_diff.s, delta.s
            );
            std::set<float> t_qs = solve_abs_eq(
                round_diff.q, delta.q,
                round_diff.s, delta.s,
                round_diff.r, delta.r
            );
            std::set<float> t_rs = solve_abs_eq(
                round_diff.r, delta.r,
                round_diff.s, delta.s,
                round_diff.q, delta.q
            );
            t_candidates.insert(t_qr.begin(), t_qr.end());
            t_candidates.insert(t_qs.begin(), t_qs.end());
            t_candidates.insert(t_rs.begin(), t_rs.end());
            
            Vec2i next_cell;
            for (float t : t_candidates) {
                if (t < t_start || t >= t_end) continue;
                next_cell = round_hex(
                    start.q + delta.q * (t + EPSILON),
                    start.r + delta.r * (t + EPSILON)
                    );
                if (next_cell != round_hex(start.q, start.r)) {
                    // Valid crossing
                    return t;
                    }
                }
        }
        return 1.0; // No crossing found, return destination.
    }

    inline std::vector<Vec2i> get_next_hexes(const Vec3f& start, const Vec3f& end, const Vec2i current, float t) {
        std::vector<Vec2i> next_hexes;
        next_hexes.reserve(2);
        
        Vec3f direction = end - start;

        Vec3f pos = start + direction * t;

        Vec3f round_pos(
            std::round(pos.q),
            std::round(pos.r),
            std::round(pos.s)
        );

        float round_diff_q = pos.q - round_pos.q;
        float round_diff_r = pos.r - round_pos.r;
        float round_diff_s = pos.s - round_pos.s;
        // If the three are not equal, return the next hex.
        if (std::abs(std::abs(round_diff_q) - std::abs(round_diff_r)) > EPSILON || 
            std::abs(std::abs(round_diff_q) - std::abs(round_diff_s)) > EPSILON || 
            std::abs(std::abs(round_diff_r) - std::abs(round_diff_s)) > EPSILON) {
            
            Vec2i hex = round_hex(
                pos.q + direction.q * EPSILON,
                pos.r + direction.r * EPSILON
            );
            //if (hex == current) return next_hexes;
            next_hexes.push_back(hex);
            return next_hexes;
        }

        // Determine the three hexes that shares the corner
        std::array<Vec2i, 3> corner_hexes;
        corner_hexes[0] = round_hex(
            pos.q - round_diff_q,
            pos.r - round_diff_r
        );
        corner_hexes[1] = round_hex(
            pos.q + round_diff_s + round_diff_r,
            pos.r - round_diff_r
        );
        corner_hexes[2] = round_hex(
            pos.q - round_diff_q,
            pos.r + round_diff_q + round_diff_s
        );
        for (const auto& hex : corner_hexes) {
            if (hex != current) {
                next_hexes.push_back(hex);
            }
        }

        return next_hexes;


    }

    inline std::vector<Vec2i> get_straight_line(Vec2i from, Vec2i to) {
        std::vector<Vec2i> result;
        result.reserve(distance(from, to));
        int N = distance(from, to);
        for (int i = 0; i <= N; i++) {
            Vec2i hex = from.lerp(to, i * (1.0f / N)).round();
            result.push_back(hex);
        }
        return result;
    }

    inline std::vector<Vec2i> get_corner_aligned_line(Vec2i from, Vec2i to) {
        std::vector<Vec2i> result;
        result.reserve(distance(from, to) * 2);
        result.push_back(from);
        Vec3f direction(to - from);

        // Bias along the perpendicular direction
        int N = distance(from, to);
        Vec2f middle;
        Vec2f bias;
        if (std::fabs(direction.q - direction.r) < EPSILON) {
            bias = Vec2f(EPSILON, -EPSILON);
        } 
        else if (std::fabs(direction.q - direction.s) < EPSILON) {
            bias = Vec2f(EPSILON, 0);
        } 
        else {
            bias = Vec2f(0, -EPSILON);
        }
        for (int i = 1; i <= N; i++) {
            middle = from.lerp(to, i * (1.0f / N));
            if (i % 2 == 0) {
                result.push_back(middle.round());
            }
            else {
                result.push_back((middle+bias).round());
                result.push_back((middle-bias).round());
            }
        }
        return result;
    }

    inline std::vector<Vec2i> hex_dda(Vec2i from, Vec2i to) {
        std::vector<Vec2i> result;
        result.push_back(from);
        
        if (from == to) return result;
        
        // Convert to cube floats for precise math
        Vec3f start(from);
        Vec3f end(to);
        Vec3f direction(to - from);

        if (direction == Vec3f(0,0,0)) {
            return result; // No movement
        }

        if (is_straight_line(direction)) {
            // Straight line optimization
            return get_straight_line(from, to);
        }

        if (is_corner_aligned(direction)) {
            // Corner-aligned optimization
            return get_corner_aligned_line(from, to);
        }
        
        Vec2i current = from;
        Vec3f pos = start;
        // Step through until we reach 'to'
        while (current != to) {
            float t = find_first_boundary_crossing(pos, end);
            
            if (t == 1.0f) {
                result.push_back(to);
                break;
            }

            std::vector<Vec2i> next_hexes = get_next_hexes(pos, end, current, t);
            for (const auto& hex : next_hexes) {
                if (hex != current)
                    result.push_back(hex);
            }

            pos = Vec3f(
                pos.q + (end.q - pos.q) * (t + EPSILON),
                pos.r + (end.r - pos.r) * (t + EPSILON)
            );

            current = round_hex(pos.q, pos.r);
        }
        
        return result;
    }

} // namespace hexmath

#endif // HEX_DDA_HPP