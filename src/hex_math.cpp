#include "hex_math.h"
#include "core/hex_math_core.hpp"
#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/vector2i.hpp>
#include <godot_cpp/variant/typed_array.hpp>

using namespace godot;

namespace {
    inline Vector2i gd(const hexmath::Vec2i& v) {
        return Vector2i(v.q, v.r);
    }

    inline Vector2 gd(const hexmath::Vec2f& v) {
        return Vector2(v.q, v.r);
    }

    inline hexmath::Vec2i hx(const Vector2i& v) {
        return hexmath::Vec2i(v.x, v.y);
    }

    inline hexmath::Vec2f hx(const Vector2& v) {
        return hexmath::Vec2f(v.x, v.y);
    }

    inline TypedArray<Vector2i> to_typed_array(const std::vector<hexmath::Vec2i>& hexes) {
        TypedArray<Vector2i> result;
        result.resize(hexes.size());
        int i = 0;
        for (const auto& hex : hexes) {
            result[i++] = gd(hex);
        }
        return result;
    }
} // namespace

HexMath::HexMath() {
    // Constructor implementation (if needed)
}

HexMath::~HexMath() {
    // Destructor implementation (if needed)
}

int HexMath::distance(Vector2i a, Vector2i b) {
    return hexmath::distance({a.x, a.y}, {b.x, b.y});
}

Vector2i HexMath::get_neighbor(Vector2i coord, int direction) {
    return gd(hexmath::get_neighbor(hx(coord), direction)); // Convert from and to Vector2i
}

TypedArray<Vector2i> HexMath::get_neighbors(Vector2i coord) {
    TypedArray<Vector2i> result;
    result.resize(6);
    for (int i = 0; i < 6; i++) {
        result[i] = gd(hexmath::get_neighbor(hx(coord), i)); // Convert from and to Vector2i
    }
    return result;
}

Vector2 HexMath::hex_to_pixel(Vector2i hex, float size) {
    return gd(hexmath::hex_to_pixel(hx(hex), size)); // Convert from and to Vector2
}

Vector2i HexMath::pixel_to_hex(Vector2 pixel, float size) {
    return gd(hexmath::pixel_to_hex(hx(pixel), size)); // Convert from and to Vector2i
}

bool HexMath::is_in_range(Vector2i coord, int radius) {
    return hexmath::is_in_range(hx(coord), radius);
}

TypedArray<Vector2i> HexMath::get_all_in_range(int radius, Vector2i center = {0,0}) {
    return to_typed_array(hexmath::get_all_in_range(radius, hx(center)));
}

TypedArray<Vector2i> HexMath::get_all_in_ring(int radius, Vector2i center = {0,0}) {
    return to_typed_array(hexmath::get_all_in_ring(radius, hx(center)));
}

TypedArray<Vector2i> HexMath::get_all_in_spiral(int radius, Vector2i center = {0,0}) {
    return to_typed_array(hexmath::get_all_in_spiral(radius, hx(center)));
}

TypedArray<Vector2i> HexMath::get_line(Vector2i from, Vector2i to = {0,0}) {
    return to_typed_array(hexmath::get_line(hx(from), hx(to)));
}

TypedArray<Vector2i> HexMath::hex_dda(Vector2i from, Vector2i to) {
    return to_typed_array(hexmath::hex_dda(hx(from), hx(to)));
}

String HexMath::helper_method() {
    return String::utf8(hexmath::helper().c_str());
}

void HexMath::_bind_methods() {
    ClassDB::bind_static_method("HexMath", D_METHOD("distance", "a", "b"), &HexMath::distance);
    ClassDB::bind_static_method("HexMath", D_METHOD("get_neighbor", "coord", "direction"), &HexMath::get_neighbor);
    ClassDB::bind_static_method("HexMath", D_METHOD("get_neighbors", "coord"), &HexMath::get_neighbors);
    ClassDB::bind_static_method("HexMath", D_METHOD("hex_to_pixel", "hex", "size"), &HexMath::hex_to_pixel);
    ClassDB::bind_static_method("HexMath", D_METHOD("pixel_to_hex", "pixel", "size"), &HexMath::pixel_to_hex);
    ClassDB::bind_static_method("HexMath", D_METHOD("is_in_range", "coord", "radius"), &HexMath::is_in_range);
    ClassDB::bind_static_method("HexMath", D_METHOD("get_all_in_range", "radius", "center"), &HexMath::get_all_in_range);
    ClassDB::bind_static_method("HexMath", D_METHOD("get_all_in_ring", "radius", "center"), &HexMath::get_all_in_ring);
    ClassDB::bind_static_method("HexMath", D_METHOD("get_all_in_spiral", "radius", "center"), &HexMath::get_all_in_spiral);
    ClassDB::bind_static_method("HexMath", D_METHOD("get_line", "from", "to"), &HexMath::get_line);
    ClassDB::bind_static_method("HexMath", D_METHOD("hex_dda", "from", "to"), &HexMath::hex_dda);
    ClassDB::bind_static_method("HexMath", D_METHOD("helper_method"), &HexMath::helper_method);
}

