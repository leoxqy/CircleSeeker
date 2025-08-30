#!/bin/bash

# CircleSeeker Conda Build and Upload Script
# Usage: ./build_and_upload.sh [conda-channel]

set -euo pipefail

PACKAGE_NAME="circleseeker"
VERSION="2.0.1"
RECIPE_DIR="$(dirname "$0")"
CONDA_CHANNEL="${1:-local}"  # Bioconda uses PRs; default to local/no upload

echo "🔨 Building CircleSeeker conda package..."
echo "Package: $PACKAGE_NAME"
echo "Version: $VERSION"
echo "Recipe: $RECIPE_DIR"
echo "Channel: $CONDA_CHANNEL"
echo

# Check dependencies
echo "🔍 Checking conda-build availability..."
if ! command -v conda-build &> /dev/null; then
    echo "❌ conda-build not found. Installing..."
    conda install conda-build -c conda-forge -y
fi

if ! command -v anaconda &> /dev/null; then
    echo "❌ anaconda-client not found. Installing..."
    conda install anaconda-client -c conda-forge -y
fi

echo

# Clean previous builds
echo "🧹 Cleaning previous builds..."
conda build purge
echo

CHANNELS=(-c conda-forge -c bioconda)

# Build the package (pick a single Python for local sanity check)
echo "🔨 Building conda package (py 3.9)..."
conda build "$RECIPE_DIR" --no-test --python=3.9 "${CHANNELS[@]}"

# Find the built package
BUILT_PACKAGE=$(conda build "$RECIPE_DIR" --output)
echo "📦 Built package: $BUILT_PACKAGE"

# Test the package
echo "🧪 Testing conda package..."
conda build "$RECIPE_DIR" --test "${CHANNELS[@]}"

echo
echo "📤 Upload step: skipped. Bioconda packages are published via PRs to bioconda-recipes, not manual uploads."

echo
echo "🎉 CircleSeeker conda build process completed!"
echo
echo "To install the package locally:"
echo "  conda install $BUILT_PACKAGE"
echo
echo "To install from a channel (after Bioconda merge):"
echo "  conda install -c conda-forge -c bioconda $PACKAGE_NAME"
