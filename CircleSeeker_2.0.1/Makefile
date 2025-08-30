# CircleSeeker convenience Makefile

.PHONY: help install dev test unit integration coverage lint format clean distclean

help:
	@echo "Common targets:"
	@echo "  make install      # pip install -e ."
	@echo "  make dev          # install dev extras (tests, lint)"
	@echo "  make test         # run all tests"
	@echo "  make unit         # run unit tests"
	@echo "  make integration  # run integration tests"
	@echo "  make coverage     # run tests with coverage"
	@echo "  make lint         # run flake8 + mypy (if available)"
	@echo "  make format       # run black + isort (if available)"
	@echo "  make clean        # remove caches and reports"

install:
	pip install -e .

dev:
	pip install -e .[dev]

test:
	./scripts/run_tests.sh

unit:
	./scripts/run_tests.sh unit

integration:
	./scripts/run_tests.sh integration

coverage:
	./scripts/run_tests.sh coverage

lint:
	@flake8 || true
	@mypy src || true

format:
	@black src tests || true
	@isort src tests || true

clean:
	rm -rf .pytest_cache __pycache__ */__pycache__ htmlcov .coverage coverage.xml
	find . -name "*.pyc" -delete
	find . -name ".DS_Store" -delete

