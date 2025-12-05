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
    static Vector2i get_neighbor(Vector2i coord, int direction);
    static TypedArray<Vector2i> get_neighbors(Vector2i coord);
    static Vector2 hex_to_pixel(Vector2i coord, float size);
    static Vector2i pixel_to_hex(Vector2 pixel, float size);
    static bool is_in_range(Vector2i coord, int radius);
    static TypedArray<Vector2i> get_all_in_range(int radius, Vector2i center);
    static TypedArray<Vector2i> get_all_in_ring(int radius, Vector2i center);
    static TypedArray<Vector2i> get_all_in_spiral(int radius, Vector2i center);
    static TypedArray<Vector2i> get_line(Vector2i from, Vector2i to);
    static TypedArray<Vector2i> hex_dda(Vector2i from, Vector2i to);


    static String helper_method();
};

} // namespace godot

#endif // HEX_MATH_H