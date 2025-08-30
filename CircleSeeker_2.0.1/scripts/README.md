# Scripts

Developer and testing utilities live here to keep the project root clean.

- `run_tests.sh`: Convenience runner around pytest (unit, integration, coverage, parallel).
- `migrate_existing_tests.sh`: Migrate legacy `test_modules` data into the pytest layout.

Usage examples:

- Run all tests: `./scripts/run_tests.sh`
- Coverage report: `./scripts/run_tests.sh coverage`
- Migrate legacy tests: `./scripts/migrate_existing_tests.sh -s /path/to/test_modules`

