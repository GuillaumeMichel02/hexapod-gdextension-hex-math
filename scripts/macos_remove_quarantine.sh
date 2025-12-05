#!/bin/bash
# Run this if macOS blocks the library with "unidentified developer" warning
xattr -dr com.apple.quarantine bin/
xattr -dr com.apple.quarantine demo/addons/hex_math/bin/
echo "Quarantine attributes removed"