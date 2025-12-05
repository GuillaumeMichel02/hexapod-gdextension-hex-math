#!/usr/bin/env python
import os
import sys

################################
#### Check dependencies
################################

if not (os.path.isdir("godot-cpp") and os.listdir("godot-cpp")):
    print("ERROR: godot-cpp submodule not initialized.")
    print("Run: git submodule update --init --recursive")
    sys.exit(1)

################################
#### Env configuration
################################

# Load godot-cpp's build environment
env = SConscript("godot-cpp/SConstruct", variant_dir="godot-cpp-build", duplicate=0)

env.Append(CPPPATH=["src/"])

# Gather source files
sources = Glob("src/*.cpp")

################################
#### Build the shared library
################################

# Define our library name
libname = "hex_math"

# Get the suffix from godot-cpp (includes platform, target, arch)
suffix = env["suffix"]

# Example: libhex_math.macos.template_debug.arm64.dylib
lib_filename = "{}{}{}{}".format(
    env.subst("$SHLIBPREFIX"),  # "lib" on Unix, "" on Windows
    libname,
    suffix,                     # ".macos.template_debug.arm64"
    env.subst("$SHLIBSUFFIX")   # ".dylib", ".dll", ".so", etc.
)

# Create platform-specific output directory
bin_dir = "bin/{}".format(env["platform"])

# Build the shared library
library = env.SharedLibrary(
    "{}/{}".format(bin_dir, lib_filename),
    source=sources,
)

################################
#### macOS code signing
################################

################################

# Tell SCons this is the default target to build
Default(library)

# Print build info (helpful for debugging)
print("Building: {} for {}, target: {}".format(
    lib_filename, env["platform"], env["target"]
))