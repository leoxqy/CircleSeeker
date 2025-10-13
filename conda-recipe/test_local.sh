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
echo "🔨 Building package..."
${BUILD_CMD} "$RECIPE_DIR" "${CHANNELS[@]}" --python=3.9 --no-anaconda-upload

# Get the built package path
BUILT_PACKAGE=$(${BUILD_CMD/ mambabuild/ build} "$RECIPE_DIR" --output --python=3.9 "${CHANNELS[@]}")

echo "📦 Built package: $BUILT_PACKAGE"

# Index the local build directory to ensure dependency resolution
echo "📇 Indexing local build directory..."
conda index "${CONDA_PREFIX}/conda-bld"

# Clean up old test environment if it exists
echo "🧹 Cleaning up old test environment..."
conda env remove -n circleseeker2_test -y 2>/dev/null || true

# Create test environment with the package and all dependencies
echo "🔬 Creating test environment and installing CircleSeeker..."
# Method 1: Direct creation with package (recommended)
conda create -n circleseeker2_test circleseeker=2.1.1 python=3.9 -y \
    -c local "${CHANNELS[@]}"

# Alternative Method 2 (if above fails): Install from file with channels
# conda create -n circleseeker2_test python=3.9 -y "${CHANNELS[@]}"
# conda install -n circleseeker2_test "$BUILT_PACKAGE" -y \
#     -c local "${CHANNELS[@]}"

# Test installation
echo "🔍 Testing installation..."
echo "✅ Testing command line tools:"
conda run -n circleseeker2_test circleseeker --help
conda run -n circleseeker2_test CircleSeeker --help

echo "✅ Testing Python import:"
conda run -n circleseeker2_test python -c "import circleseeker; print('Version:', circleseeker.__version__)"
conda run -n circleseeker2_test python -c "from circleseeker.modules import ecc_summary; print('Module test passed')"

echo "✅ Checking pip dependencies:"
conda run -n circleseeker2_test python -m pip check

# Optional: List installed packages for verification
echo ""
echo "📋 Key packages installed:"
conda run -n circleseeker2_test conda list | grep -E "(circleseeker|click|pandas|numpy|biopython|pysam)" || true

# Cleanup
echo ""
echo "🧹 Cleaning up test environment..."
conda remove -n circleseeker2_test --all -y

echo ""
echo "✅ Local test completed successfully!"
echo ""
echo "Package is ready for distribution: $BUILT_PACKAGE"
