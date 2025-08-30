#!/bin/bash

# CircleSeeker2 Local Test Script
# Tests the conda package locally before uploading

set -euo pipefail

RECIPE_DIR="$(dirname "$0")"

echo "🧪 Testing CircleSeeker conda package locally..."

# Prefer boa/mambabuild if available for speed
if command -v conda-mambabuild &>/dev/null; then
  BUILD_CMD="conda mambabuild"
else
  BUILD_CMD="conda build"
fi

CHANNELS=(-c conda-forge -c bioconda)

# Build and test
${BUILD_CMD} "$RECIPE_DIR" "${CHANNELS[@]}" --python=3.9 --no-anaconda-upload

# Get the built package path
BUILT_PACKAGE=$(${BUILD_CMD/ mambabuild/ build} "$RECIPE_DIR" --output --python=3.9 "${CHANNELS[@]}")

echo "📦 Built package: $BUILT_PACKAGE"

# Create test environment and install using conda run (no shell activation needed)
echo "🔬 Creating test environment..."
conda create -n circleseeker2_test python=3.9 -y "${CHANNELS[@]}"

echo "📥 Installing package in test environment..."
conda run -n circleseeker2_test conda install -y "$BUILT_PACKAGE" "${CHANNELS[@]}"

# Test installation
echo "🔍 Testing installation..."
conda run -n circleseeker2_test circleseeker --help
conda run -n circleseeker2_test python -c "import circleseeker; print('✅ Import successful, version:', circleseeker.__version__)"

# Cleanup
conda remove -n circleseeker2_test --all -y

echo "✅ Local test completed successfully!"
