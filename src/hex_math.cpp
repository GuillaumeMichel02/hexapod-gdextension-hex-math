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
    return hexmath::distance(hx(a), hx(b));
}

float HexMath::distance_to_line(Vector2i coord, Vector2i from, Vector2i to = {0,0}) {
    return hexmath::distance_to_line(hx(coord), hx(from), hx(to));
}

float HexMath::distance_to_line_border(Vector2i coord, Vector2i from, Vector2i to = {0,0}) {
    return hexmath::distance_to_line_border(hx(coord), hx(from), hx(to));
}

Vector2 HexMath::project(Vector2i coord, Vector2i line) {
    return gd(hexmath::project(hx(coord), hx(line)));
}

int HexMath::cross_product(Vector2i a, Vector2i b) {
    return hexmath::cross_product(hx(a), hx(b));
}

float HexMath::dot_product(Vector2i a, Vector2i b) {
    return hexmath::dot_product(hx(a), hx(b));
}

float HexMath::get_angle(Vector2i a, Vector2i b) {
    return hexmath::get_angle(hx(a), hx(b));
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

std::array<Vector2, 6> HexMath::get_hex_corners(Vector2i center, float size) {
    auto corners = hexmath::get_hex_corners(hx(center), size);
    std::array<Vector2, 6> gd_corners;
    for (size_t i = 0; i < corners.size(); i++) {
        gd_corners[i] = gd(corners[i]);
    }
    return gd_corners;
}

bool HexMath::is_in_range(Vector2i coord, int radius) {
    return hexmath::is_in_range(hx(coord), radius);
}

bool HexMath::is_in_cone(Vector2i coord, int direction, int angle_width, Vector2i center) {
    return hexmath::is_in_cone(hx(coord), direction, angle_width, hx(center));
}

bool HexMath::is_in_ring(Vector2i coord, int radius, Vector2i center) {
    return hexmath::is_in_ring(hx(coord), radius, hx(center));
}

bool HexMath::is_in_supercover(Vector2i coord, Vector2i from, Vector2i to) {
    return hexmath::is_in_supercover(hx(coord), hx(from), hx(to));
}

TypedArray<Vector2i> HexMath::get_all_in_range(int radius, Vector2i center = {0,0}) {
    return to_typed_array(hexmath::get_all_in_range(radius, hx(center)));
}

TypedArray<Vector2i> HexMath::get_all_in_ring(int radius, int width = 6, int init_direction = 0, Vector2i center = {0,0}) {
    return to_typed_array(hexmath::get_all_in_ring(radius, width, init_direction, hx(center)));
}

TypedArray<Vector2i> HexMath::get_all_in_spiral(int radius, int width = 6, int init_direction = 0, Vector2i center = {0,0}) {
    return to_typed_array(hexmath::get_all_in_spiral(radius, width, init_direction, hx(center)));
}

TypedArray<Vector2i> HexMath::get_all_in_line(Vector2i from, Vector2i to = {0,0}) {
    return to_typed_array(hexmath::get_all_in_line(hx(from), hx(to)));
}

TypedArray<Vector2i> HexMath::get_all_in_supercover(Vector2i from, Vector2i to) {
    return to_typed_array(hexmath::get_all_in_supercover(hx(from), hx(to)));
}

Vector2i HexMath::rotate(Vector2i coord, int steps, Vector2i center) {
    return gd(hexmath::rotate(hx(coord), steps, hx(center)));
}

Vector2i HexMath::reflect(Vector2i coord, int axis, Vector2i center) {
    return gd(hexmath::reflect(hx(coord), axis, hx(center)));
}

String HexMath::helper_method() {
    return String::utf8(hexmath::helper().c_str());
}

void HexMath::_bind_methods() {
    ClassDB::bind_static_method("HexMath", D_METHOD("distance", "a", "b"), &HexMath::distance);
    ClassDB::bind_static_method("HexMath", D_METHOD("distance_to_line", "coord", "from", "to"), &HexMath::distance_to_line);
    ClassDB::bind_static_method("HexMath", D_METHOD("distance_to_line_border", "coord", "from", "to"), &HexMath::distance_to_line_border); 
    ClassDB::bind_static_method("HexMath", D_METHOD("project", "coord", "line"), &HexMath::project);
    ClassDB::bind_static_method("HexMath", D_METHOD("cross_product", "a", "b"), &HexMath::cross_product);
    ClassDB::bind_static_method("HexMath", D_METHOD("dot_product", "a", "b"), &HexMath::dot_product);
    ClassDB::bind_static_method("HexMath", D_METHOD("get_angle", "a", "b"), &HexMath::get_angle);
    ClassDB::bind_static_method("HexMath", D_METHOD("get_neighbor", "coord", "direction"), &HexMath::get_neighbor);
    ClassDB::bind_static_method("HexMath", D_METHOD("get_neighbors", "coord"), &HexMath::get_neighbors);
    ClassDB::bind_static_method("HexMath", D_METHOD("hex_to_pixel", "hex", "size"), &HexMath::hex_to_pixel);
    ClassDB::bind_static_method("HexMath", D_METHOD("pixel_to_hex", "pixel", "size"), &HexMath::pixel_to_hex);
    ClassDB::bind_static_method("HexMath", D_METHOD("get_hex_corners", "center", "size"), &HexMath::get_hex_corners);
    ClassDB::bind_static_method("HexMath", D_METHOD("is_in_range", "coord", "radius"), &HexMath::is_in_range);
    ClassDB::bind_static_method("HexMath", D_METHOD("is_in_cone", "coord", "direction", "angle_width", "center"), &HexMath::is_in_cone);
    ClassDB::bind_static_method("HexMath", D_METHOD("is_in_ring", "coord", "radius", "center"), &HexMath::is_in_ring);
    ClassDB::bind_static_method("HexMath", D_METHOD("is_in_supercover", "coord", "from", "to"), &HexMath::is_in_supercover);
    ClassDB::bind_static_method("HexMath", D_METHOD("get_all_in_range", "radius", "center"), &HexMath::get_all_in_range);
    ClassDB::bind_static_method("HexMath", D_METHOD("get_all_in_ring", "radius", "width", "init_direction", "center"), &HexMath::get_all_in_ring);
    ClassDB::bind_static_method("HexMath", D_METHOD("get_all_in_spiral", "radius", "width", "init_direction", "center"), &HexMath::get_all_in_spiral);
    ClassDB::bind_static_method("HexMath", D_METHOD("get_all_in_line", "from", "to"), &HexMath::get_all_in_line);
    ClassDB::bind_static_method("HexMath", D_METHOD("get_all_in_supercover", "from", "to"), &HexMath::get_all_in_supercover);
    ClassDB::bind_static_method("HexMath", D_METHOD("rotate", "coord", "steps", "center"), &HexMath::rotate);
    ClassDB::bind_static_method("HexMath", D_METHOD("reflect", "coord", "axis", "center"), &HexMath::reflect);
    ClassDB::bind_static_method("HexMath", D_METHOD("helper_method"), &HexMath::helper_method);
}

