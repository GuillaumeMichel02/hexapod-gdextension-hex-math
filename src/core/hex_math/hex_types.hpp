// hex_types.hpp
// This file defines the primary data structures for hexagonal grids.
#ifndef HEX_TYPES_HPP
#define HEX_TYPES_HPP

#include <cstdlib>
#include <cmath>

namespace hexmath {

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

} // namespace hexmath

#endif // HEX_TYPES_HPP