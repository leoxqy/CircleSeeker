#!/bin/bash
# CircleSeeker 测试运行脚本

echo "=== CircleSeeker 测试运行器 ==="

# Ensure src/ is importable without editable install
export PYTHONPATH="$(pwd)/src:${PYTHONPATH:-}"

case "$1" in
    "demo")
        echo "运行演示测试..."
        pytest -c scripts/pytest.quick.ini tests/test_demo.py -v
        ;;
    "unit")
        echo "运行单元测试..."
        pytest -c scripts/pytest.quick.ini tests/test_*.py -v -m "not integration"
        ;;
    "integration") 
        echo "运行集成测试..."
        pytest -c scripts/pytest.quick.ini tests/test_*.py -v -m integration
        ;;
    "coverage")
        echo "运行测试并生成覆盖率报告..."
        pytest -c scripts/pytest.quick.ini tests/test_demo.py --cov=circleseeker --cov-report=html --cov-report=term-missing
        echo "覆盖率报告生成在: htmlcov/index.html"
        ;;
    "quick")
        echo "快速运行（跳过慢速测试）..."
        pytest -m "not slow" -v
        ;;
    "parallel")
        echo "并行运行测试..."
        pytest -c scripts/pytest.quick.ini tests/test_demo.py -n auto -v
        ;;
    "check")
        echo "检查测试发现..."
        pytest -c scripts/pytest.quick.ini --collect-only
        ;;
    *)
        echo "运行演示测试..."
        echo "用法: $0 [demo|unit|integration|coverage|parallel|quick|check]"
        echo ""
        pytest -c scripts/pytest.quick.ini tests/test_demo.py -v
        ;;
esac
