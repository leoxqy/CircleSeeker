# Contributing to CircleSeeker

Thank you for your interest in contributing to CircleSeeker! This document provides guidelines and instructions for contributing to the project.

## Code of Conduct

By participating in this project, you agree to abide by our Code of Conduct:
- Be respectful and inclusive
- Welcome newcomers and help them get started
- Focus on constructive criticism
- Respect differing viewpoints and experiences

## How to Contribute

### Reporting Issues

1. **Check existing issues** first to avoid duplicates
2. **Use the issue template** when available
3. **Provide detailed information**:
   - CircleSeeker version (`circleseeker --version`)
   - Python version (`python --version`)
   - Operating system and version
   - Complete error messages and stack traces
   - Minimal reproducible example if possible

### Suggesting Features

1. Open an issue with the "enhancement" label
2. Clearly describe the feature and its benefits
3. Provide use cases and examples
4. Be open to discussion and feedback

### Contributing Code

#### Setting Up Development Environment

```bash
# Fork and clone the repository
git clone https://github.com/YOUR_USERNAME/circleseeker2.git
cd circleseeker2

# Create a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode with all dependencies
pip install -e ".[dev,test,docs]"

# Install pre-commit hooks
pre-commit install
```

#### Development Workflow

1. **Create a branch** for your feature or fix:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes** following our coding standards

3. **Write tests** for new functionality:
   ```bash
   # Run tests
   pytest tests/
   
   # Run with coverage
   pytest --cov=circleseeker tests/
   ```

4. **Format and lint** your code:
   ```bash
   # Format with black
   black src/circleseeker tests/
   
   # Sort imports
   isort src/circleseeker tests/
   
   # Check with flake8
   flake8 src/circleseeker tests/
   
   # Type checking
   mypy src/circleseeker
   ```

5. **Update documentation** if needed

6. **Commit your changes**:
   ```bash
   git add .
   git commit -m "feat: add new feature X"
   ```

7. **Push and create a pull request**:
   ```bash
   git push origin feature/your-feature-name
   ```

### Commit Message Guidelines

We follow [Conventional Commits](https://www.conventionalcommits.org/):

- `feat:` New features
- `fix:` Bug fixes
- `docs:` Documentation changes
- `test:` Test additions or modifications
- `refactor:` Code refactoring
- `perf:` Performance improvements
- `style:` Code style changes (formatting, etc.)
- `chore:` Maintenance tasks

Examples:
```
feat: add support for paired-end reads
fix: correct coverage calculation in Contortionist module
docs: update installation instructions for macOS
test: add unit tests for Carousel module
```

### Code Style Guidelines

1. **Follow PEP 8** with these modifications:
   - Line length: 100 characters
   - Use double quotes for strings

2. **Use type hints** for all functions:
   ```python
   def process_reads(input_file: Path, min_length: int = 100) -> List[Dict[str, Any]]:
       """Process reads from input file."""
       pass
   ```

3. **Write comprehensive docstrings**:
   ```python
   def calculate_coverage(
       bam_file: Path,
       bed_file: Path,
       window_size: int = 50
   ) -> pd.DataFrame:
       """Calculate coverage statistics for regions.
       
       Args:
           bam_file: Path to sorted BAM file
           bed_file: Path to BED file with regions
           window_size: Window size for coverage calculation
           
       Returns:
           DataFrame with coverage statistics per region
           
       Raises:
           FileNotFoundError: If input files don't exist
           ValueError: If window_size is not positive
       """
   ```

4. **Use meaningful variable names** and avoid abbreviations

5. **Add logging** for important operations:
   ```python
   logger.info(f"Processing {len(reads)} reads")
   logger.debug(f"Read {read_id} passed quality filter")
   ```

### Testing Guidelines

1. **Write tests for all new features**
2. **Aim for >80% code coverage**
3. **Use pytest fixtures** for test data
4. **Test edge cases** and error conditions
5. **Mock external tool calls** when appropriate

Example test:
```python
def test_carousel_process_simple_reads(tmp_path):
    """Test Carousel processing of simple tandem repeats."""
    input_file = tmp_path / "input.txt"
    input_file.write_text(">read1\nACGT\n2\t4\t10\t0.95\n")
    
    carousel = Carousel()
    result = carousel.process(input_file)
    
    assert len(result) == 1
    assert result[0]["read_id"] == "read1"
```

### Documentation

1. **Update relevant documentation** when changing functionality
2. **Add docstrings** to all public functions and classes
3. **Include examples** in docstrings when helpful
4. **Update README.md** for user-facing changes

## Pull Request Process

1. **Ensure all tests pass** locally
2. **Update documentation** as needed
3. **Add entry to CHANGELOG.md** for significant changes
4. **Request review** from maintainers
5. **Address review feedback** promptly
6. **Squash commits** if requested before merging

## Questions?

If you have questions about contributing:
1. Check our [documentation](https://circleseeker2.readthedocs.io)
2. Open a [discussion](https://github.com/circleseeker/circleseeker2/discussions)
3. Contact the maintainers

Thank you for contributing to CircleSeeker!