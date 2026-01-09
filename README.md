# GDE Hex Math

A GDExtension for Godot 4.5 providing comprehensive hexagonal grid mathematical operations. It is designed to handle hex grid calculations in C++ for efficiency and reusability, based on the axial coordinates system.

It is currently only intended for pointy-top hexagons, and specifically tailored for the Hexapod project.

The `main` branch is the host of the most up-to-date stable version of the extension.

## Project architecture

- `demo` : Small Godot project intended to test the extension. The built library must be included manually.
- `godot-cpp` : GDExtension submodule, linked to the official godot extension repository.
- `scripts` : Utilitary scripts to manually build the library during development and help remove the macos quarantine.
- `src` :
  - `register_types.cpp register_types.h` : Mandatory files to define the extension.
  - `core` : Defines all the mathematical methods for the extension, with no link to Godot. The `core` may be used as-is in other C++ libraries that need hexaxonal grid mathematical operations.
  - `hex_math.cpp hex_math.h` : Defines the HexMath class, that creates a bridge between the core methods to the extension methods.
  - `register_types.cpp register_types.h` : Register the HexMath class to the GDExtension.

## Todos

- **Continuous development** **:** New functions or optimization of old ones might occur as the project for which this extension was made is still in development.
- **CI/CD** **:** Releases in the `main` branch should initiate a GitHub Action to compile the extension into a ready-to-use hexmath addon folder.
- **Settings :** The extension will get parameters to handle flat-top orientation, and inverted axis.
- **Documentation :** The methods will get a Godot compatible documentation to help use the extension.
- **Unit testing :** For future releases, tests the methods to ensure that nothing breaks before pushing.
