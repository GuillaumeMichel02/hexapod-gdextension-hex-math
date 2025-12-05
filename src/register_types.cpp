#include "register_types.h"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/defs.hpp>
#include <godot_cpp/godot.hpp>
#include "hex_math.h"


using namespace godot;

void initialize_hex_math_module(ModuleInitializationLevel p_level) {
    ClassDB::register_class<HexMath>();
}

void uninitialize_hex_math_module(ModuleInitializationLevel p_level) {
    // Cleanup if needed (usually empty)
}

extern "C" {
    // This is the entry point Godot looks for
    // The name MUST match entry_symbol in the .gdextension file
    GDExtensionBool GDE_EXPORT hex_math_library_init(
        GDExtensionInterfaceGetProcAddress p_get_proc_address,
        const GDExtensionClassLibraryPtr p_library,
        GDExtensionInitialization *r_initialization
    ) {
        godot::GDExtensionBinding::InitObject init_obj(p_get_proc_address, p_library, r_initialization);

        init_obj.register_initializer(initialize_hex_math_module);
        init_obj.register_terminator(uninitialize_hex_math_module);
        init_obj.set_minimum_library_initialization_level(MODULE_INITIALIZATION_LEVEL_SCENE);

        return init_obj.init();
    }
}