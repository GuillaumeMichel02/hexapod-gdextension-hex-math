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
    constexpr float PI = 3.14159265359f;

    /* ###############################################
        Forward Declarations
    ############################################### */
    //struct Vec3i; // Forward declaration
    struct Vec2f;
    struct Vec2i;
    inline Vec2i round_hex(float q, float r);
    inline float dot_product(const Vec2i& a, const Vec2i& b);
    inline float dot_product(const Vec2f& a, const Vec2f& b);
    /* ###############################################
        Type Definitions
    ############################################### */


    struct Vec2i { 
        int q, r; 
        Vec2i() : q(0), r(0) {}
        Vec2i(int q, int r) : q(q), r(r) {}

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
        Vec2f normalized() const;
        Vec2f project_onto(const Vec2i& other) const;

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

        Vec2f operator*(float scalar) const {
            return Vec2f(q * scalar, r * scalar);
        }

        Vec2i round() const {
            return round_hex(q, r);
        }

        float length() const {
            return (std::abs(q) + std::abs(r) + std::abs(-q - r)) / 2;
        }

        Vec2f project_onto(const Vec2f& other) const {
            float dot = dot_product(*this, other);
            float other_len_sq = other.length() * other.length();
            if (other_len_sq == 0) return Vec2f(0,0);
            float scalar = dot / other_len_sq;
            return Vec2f(other.q * scalar, other.r * scalar);
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

    inline Vec2f Vec2i::normalized() const {
        float len = this->length();
        if (len == 0) return Vec2f(0,0);
        return Vec2f(q * (1.0f / len), r * (1.0f / len));
    }

    inline Vec2f Vec2i::project_onto(const Vec2i& other) const {
        return Vec2f(Vec2f(*this).project_onto(Vec2f(other)));
    }

    inline const Vec2i DIRECTIONS[6] = {Vec2i(1,0), Vec2i(1,-1), Vec2i(0,-1), Vec2i(-1,0), Vec2i(-1,1),Vec2i(0,1)};
    inline const Vec2f CORNERS[6] = {
        Vec2f(DIRECTIONS[0] + DIRECTIONS[1]) * (1.0f/3),
        Vec2f(DIRECTIONS[1] + DIRECTIONS[2]) * (1.0f/3),
        Vec2f(DIRECTIONS[2] + DIRECTIONS[3]) * (1.0f/3),
        Vec2f(DIRECTIONS[3] + DIRECTIONS[4]) * (1.0f/3),
        Vec2f(DIRECTIONS[4] + DIRECTIONS[5]) * (1.0f/3),
        Vec2f(DIRECTIONS[5] + DIRECTIONS[0]) * (1.0f/3)
    };
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


    /* ###############################################
        Helper methods
    ############################################### */

    // Helper : strings out whatever I want to see in the editor
    inline std::string helper() {
        Vec2i from(0,0);
        Vec2i to(2,5);
        std::string res = "Hex DDA from (0,0) to (2,6):\n";
        std::vector<Vec2i> result;
        result.push_back(from);
        
        // Convert to cube floats for precise math
        Vec3f start(from);
        Vec3f end(to);
        Vec3f direction(to - from);
        
        Vec2i current = from;
        Vec3f pos = start;
        // Step through until we reach 'to'
        while (current != to) {
            float t = find_first_boundary_crossing(pos, end);
            //get string of current, pos_ante, end, t
            res += "Current Hex: (" + std::to_string(current.q) + "," + std::to_string(current.r) + ")\n";
            res += "Position: (" + std::to_string(pos.q) + "," + std::to_string(pos.r) + ")\n";
            res += "End: (" + std::to_string(end.q) + "," + std::to_string(end.r) + ")\n";
            res += "t: " + std::to_string(t) + "\n";
            if (t == 1.0f) {
                result.push_back(to);
                break;
            }

            std::vector<Vec2i> next_hexes = get_next_hexes(pos, end, current, t);
            for (const auto& hex : next_hexes) {
                // get string of hex
                res += "Next Hex: (" + std::to_string(hex.q) + "," + std::to_string(hex.r) + ")\n";
                if (hex != current)
                    result.push_back(hex);
            }

            pos = Vec3f(
                pos.q + (end.q - pos.q) * (t + EPSILON),
                pos.r + (end.r - pos.r) * (t + EPSILON)
            );
            res += "New Position: (" + std::to_string(pos.q) + "," + std::to_string(pos.r) + ")\n";

            current = round_hex(pos.q, pos.r);
            res += "New Current Hex: (" + std::to_string(current.q) + "," + std::to_string(current.r) + ")\n\n";
        }
        
        return res;

    }
} // namespace hexmath

#endif // HEX_MATH_CORE_H