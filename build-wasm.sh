#!/bin/bash
set -e

echo "Building brainwires-compute-engine for WebAssembly..."

# Check if wasm-pack is installed
if ! command -v wasm-pack &> /dev/null; then
    echo "wasm-pack not found. Installing..."
    cargo install wasm-pack
fi

# Build for different targets
echo ""
echo "Building for bundler (webpack, rollup, vite - for client-side)..."
wasm-pack build --target bundler --features wasm

echo ""
echo "Building for Node.js (for server-side)..."
# Build to temp directory then move to pkg/nodejs/
wasm-pack build --target nodejs --features wasm

# Move nodejs build to subdirectory
echo "Organizing Node.js output..."
mkdir -p pkg/nodejs
cp pkg/computational_engine.js pkg/nodejs/
cp pkg/computational_engine.d.ts pkg/nodejs/
cp pkg/computational_engine_bg.wasm pkg/nodejs/
cp pkg/computational_engine_bg.wasm.d.ts pkg/nodejs/
cp pkg/package.json pkg/nodejs/
cp pkg/README.md pkg/nodejs/ 2>/dev/null || true

# Restore bundler target (it was overwritten)
echo "Restoring bundler target..."
wasm-pack build --target bundler --features wasm

echo ""
echo "✓ WASM build complete!"
echo ""
echo "Output directories:"
echo "  - pkg/          (bundler - for client-side webpack/rollup/vite)"
echo "  - pkg/nodejs/   (nodejs - for server-side API routes)"
echo ""
echo "Usage in Next.js:"
echo "  Client: import from '../rust/brainwires-compute-engine/pkg/computational_engine.js'"
echo "  Server: import from '../rust/brainwires-compute-engine/pkg/nodejs/computational_engine.js'"
