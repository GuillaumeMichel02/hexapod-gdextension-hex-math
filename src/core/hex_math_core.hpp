#ifndef HEX_MATH_CORE_H
#define HEX_MATH_CORE_H

#include <cstdlib>
#include <cmath>
#include <array>
#include <vector>
#include <limits>
#include <set>
#include <string>

namespace hexmath {

    constexpr float EPSILON = 0.000001f;
    constexpr float SQRT_3 = 1.73205080757f;

    /* ###############################################
        Forward Declarations
    ############################################### */
    //struct Vec3i; // Forward declaration
    struct Vec2f;
    struct Vec2i;
    inline Vec2i round_hex(float q, float r);

    /* ###############################################
        Type Definitions
    ############################################### */


    struct Vec2i { 
        int q, r; 
        Vec2i() : q(0), r(0) {}
        Vec2i(int q, int r) : q(q), r(r) {}
        //Vec2i(Vec3i a) : q(a.q), r(a.r) {}

        Vec2i operator+(const Vec2i& other) const {
            return Vec2i(q + other.q, r + other.r);
        }

        Vec2i operator-(const Vec2i& other) const {
            return Vec2i(q - other.q, r - other.r);
        }

        Vec2i operator*(int scalar) const {
            return Vec2i(q * scalar, r * scalar);
        }

        bool operator==(const Vec2i& other) const {
            return q == other.q && r == other.r;
        }

        bool operator!=(const Vec2i& other) const {
            return !(*this == other);
        }

        int length() const {
            return (std::abs(q) + std::abs(r) + std::abs(-q - r)) / 2;
        }

        Vec2f lerp(const Vec2i other, float t) const;

    };

    struct Vec2f { 
        float q, r; 
        Vec2f() : q(0), r(0) {}
        Vec2f( float q, float r ) : q(q), r(r) {}
        Vec2f( Vec2i a) : q(a.q), r(a.r) {}

        Vec2f operator+(const Vec2f& other) const {
            return Vec2f(q + other.q, r + other.r);
        }

        Vec2f operator-(const Vec2f& other) const {
            return Vec2f(q - other.q, r - other.r);
        }

        Vec2i round() const {
            return round_hex(q, r);
        }
    };

    struct Vec3f {
        float q, r, s;
        Vec3f() : q(0), r(0), s(0) {}
        Vec3f( float q, float r ) : q(q), r(r), s(-q - r) {}
        Vec3f( Vec2f a ) : q(a.q), r(a.r), s(-a.q - a.r) {}
        Vec3f( Vec2i a ) : q(a.q), r(a.r), s(-a.q - a.r) {}
        Vec3f(float q, float r, float s) : q(q), r(r), s(s) {}

        Vec3f operator+(const Vec3f& other) const {
            return Vec3f(q + other.q, r + other.r, s + other.s);
        }

        Vec3f operator-(const Vec3f& other) const {
            return Vec3f(q - other.q, r - other.r, s - other.s);
        }

        Vec3f operator*(float scalar) const {
            return Vec3f(q * scalar, r * scalar, s * scalar);
        }

        bool operator==(const Vec3f& other) const {
            return q == other.q && r == other.r && s == other.s;
        }

        bool operator!=(const Vec3f& other) const {
            return !(*this == other);
        }
    };

    /* ###############################################
        Vec2i Methods
    ############################################### */

    inline Vec2f Vec2i::lerp(const Vec2i other, float t) const {
        return Vec2f(q + (other.q - q) * t, r + (other.r - r) * t);
    }

    // struct Vec3i {
    //     int q, r, s;
    //     Vec3i() : q(0), r(0), s(0) {}
    //     Vec3i(int q, int r) : q(q), r(r), s(-q - r) {}
    //     Vec3i(int q, int r, int s) : q(q), r(r), s(s) {}
    //     Vec3i(Vec2i a) : q(a.q), r(a.r), s(-a.q - a.r) {}

    // };

    /* ###############################################
        Relation methods
    ############################################### */

    inline int distance(Vec2i a, Vec2i b) {
        return (a - b).length();
    }

    /* ###############################################
        Hex algorithm utility methods
    ############################################### */

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
        // Step through until we reach 'to'
        while (current != to) {
            float t = find_first_boundary_crossing(pos, end);
            
            if (t == 1.0f) {
                result.push_back(to);
                break;
            }

            Vec2i hex = round_hex(
                pos.q + (end.q - pos.q) * (t + EPSILON),
                pos.r + (end.r - pos.r) * (t + EPSILON)
            );

            pos = Vec3f(
                pos.q + (end.q - pos.q) * (t + EPSILON),
                pos.r + (end.r - pos.r) * (t + EPSILON)
            );

            if (hex != current) {
                result.push_back(hex);
                current = hex;
            }
        }
        
        return result;
    }

    /* ###############################################
        Neighbor methods
    ############################################### */

    inline const Vec2i DIRECTIONS[6] = {Vec2i(1,0), Vec2i(1,-1), Vec2i(0,-1), Vec2i(-1,0), Vec2i(-1,1),Vec2i(0,1)};

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
        return round_hex(q, r); // Placeholder implementation
    }

    /* ###############################################
        Range methods
    ############################################### */

    inline bool is_in_range(Vec2i coord, int radius) {
        return coord.length() <= radius;
    }

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

    inline std::vector<Vec2i> get_all_in_ring(int radius, Vec2i center = Vec2i()) {
        std::vector<Vec2i> in_ring;

        if (radius == 0) {
            in_ring.push_back(center);
            return in_ring;
        }

        for (int i = 0; i < 6; i++) {
                Vec2i corner = center + DIRECTIONS[i] * radius;
                Vec2i direction = DIRECTIONS[(i+2)%6];
            for (int j = 0; j < radius; j++) {
                in_ring.push_back(corner + direction * j);
            }
        }
        return in_ring;
    }

    inline std::vector<Vec2i> get_all_in_spiral(int radius, Vec2i center = Vec2i()) {
        std::vector<Vec2i> in_spiral;
        for (int i = 0; i <= radius; i++) {
            for (const auto& hex : get_all_in_ring(i, center)) {
                in_spiral.push_back(hex);
            }
        }
        return in_spiral;
    }

    inline std::vector<Vec2i> get_line(Vec2i from, Vec2i to = Vec2i()) {
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

    // Helper : strings out whatever I want to see in the editor
    inline std::string helper() {
        Vec2i from(0,0);
        Vec2i to(2,6);

        std::vector<Vec2i> result;
        result.push_back(from);
        
        // Convert to cube floats for precise math
        Vec3f start(from);
        Vec3f end(to);
        Vec3f direction(to - from);
        
        // Current position (starts at from)
        Vec3f pos = start;
        Vec2i current = from;
        std::string res;
        // Step through until we reach 'to'
        while (current != to) {
            float t = find_first_boundary_crossing(pos, end);
            
            if (t == 1.0f) {
                result.push_back(to);
                break;
            }

            Vec2i hex = round_hex(
                pos.q + (end.q - pos.q) * (t + EPSILON),
                pos.r + (end.r - pos.r) * (t + EPSILON)
            );

            pos = Vec3f(
                pos.q + (end.q - pos.q) * (t + EPSILON),
                pos.r + (end.r - pos.r) * (t + EPSILON)
            );



            if (hex != current) {
                result.push_back(hex);
                current = hex;
            }
            res += "pos: (" + std::to_string(pos.q) + ", " + std::to_string(pos.r) + ", " + std::to_string(pos.s) + ")\n";
            res += "end: (" +  std::to_string(end.q) + ", " + std::to_string(end.r) + ", " + std::to_string(end.s) + ")\n";
            res += "t: " + std::to_string(t) + "\n";
            res += "hex: (" + std::to_string(hex.q) + ", " + std::to_string(hex.r) + ")\n";
            res+= "----\n";

        }


        return res;
    }

}

#endif // HEX_MATH_CORE_H