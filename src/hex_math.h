#ifndef HEX_MATH_H
#define HEX_MATH_H

#include <godot_cpp/classes/object.hpp>
#include <godot_cpp/variant/vector2i.hpp>
#include <godot_cpp/variant/vector2.hpp>
#include <godot_cpp/variant/typed_array.hpp>
#include <godot_cpp/variant/string.hpp>
#include "core/hex_math_core.hpp"

namespace godot {

class HexMath : public Object {
    GDCLASS(HexMath, Object);

protected:
    static void _bind_methods();

public:
    HexMath();
    ~HexMath();
    static int distance(Vector2i a, Vector2i b);
    static float distance_to_line(Vector2i coord, Vector2i from, Vector2i to);
    static float distance_to_line_border(Vector2i coord, Vector2i from, Vector2i to);
    static Vector2 project(Vector2i coord, Vector2i line);
    static int cross_product(Vector2i a, Vector2i b);
    static float dot_product(Vector2i a, Vector2i b);
    static float get_angle(Vector2i a, Vector2i b);
    static Vector2i get_neighbor(Vector2i coord, int direction);
    static TypedArray<Vector2i> get_neighbors(Vector2i coord);
    static Vector2 hex_to_pixel(Vector2i coord, float size);
    static Vector2i pixel_to_hex(Vector2 pixel, float size);
    static TypedArray<Vector2> get_hex_corners(Vector2i center, float size);
    static bool is_in_range(Vector2i coord, int radius);
    static bool is_in_cone(Vector2i coord, int direction, int angle_width, Vector2i center);
    static bool is_in_ring(Vector2i coord, int radius, Vector2i center);
    static bool is_in_supercover(Vector2i coord, Vector2i from, Vector2i to);
    static TypedArray<Vector2i> get_all_in_range(int radius, Vector2i center);
    static TypedArray<Vector2i> get_all_in_ring(int radius, int width, int init_direction, Vector2i center);
    static TypedArray<Vector2i> get_all_in_spiral(int radius, int width, int init_direction, Vector2i center);
    static TypedArray<Vector2i> get_all_in_line(Vector2i from, Vector2i to);
    static TypedArray<Vector2i> get_all_in_supercover(Vector2i from, Vector2i to);
    static Vector2i rotate(Vector2i coord, int steps, Vector2i center);
    static Vector2i reflect(Vector2i coord, int axis, Vector2i center);

    static String helper_method();
};

} // namespace godot

#endif // HEX_MATH_H