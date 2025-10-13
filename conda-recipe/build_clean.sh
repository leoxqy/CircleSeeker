#!/bin/bash
# Clean build script for CircleSeeker - suppresses numpy warning for noarch packages

set -euo pipefail

RECIPE_DIR="$(dirname "$0")"

echo "🧪 Building CircleSeeker conda package..."

# Clean previous builds
conda build purge

# Build with explicit numpy version to suppress warning
# Note: For noarch packages, numpy version doesn't matter
CONDA_NPY=126 conda build "$RECIPE_DIR" \
    -c conda-forge \
    -c bioconda \
    --python=3.9 \
    --no-anaconda-upload \
    2>&1 | grep -v "WARNING: No numpy version specified"

# Get the built package path
BUILT_PACKAGE=$(conda build "$RECIPE_DIR" --output --python=3.9 -c conda-forge -c bioconda)

echo "📦 Built package: $BUILT_PACKAGE"
echo "✅ Build completed successfully!"
echo ""
echo "To test the package:"
echo "  conda install $BUILT_PACKAGE -c conda-forge -c bioconda"