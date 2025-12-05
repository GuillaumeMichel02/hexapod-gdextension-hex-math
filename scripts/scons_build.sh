scons platform=macos target=template_debug
# Replace the bin folder content in the demo
rm -rf demo/addons/hex_math/bin/*
cp -r bin/macos/* demo/addons/hex_math/bin/