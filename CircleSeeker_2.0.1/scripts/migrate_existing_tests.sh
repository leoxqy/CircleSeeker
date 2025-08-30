#!/bin/bash
# CircleSeeker 测试迁移脚本
# 将现有的 test_modules 迁移到 pytest 格式

set -e  # 出错时退出

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 项目路径（自动定位为本脚本所在目录的上级）
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
TARGET_DIR="$PROJECT_ROOT/tests"

# 源数据目录（通过参数 -s 指定），例如：scripts/migrate_existing_tests.sh -s /path/to/test_modules
SOURCE_DIR=""
while getopts ":s:" opt; do
  case $opt in
    s) SOURCE_DIR="$OPTARG" ;;
    *) ;;
  esac
done

echo -e "${BLUE}=== CircleSeeker 测试迁移工具 ===${NC}"
echo ""

if [ -z "$SOURCE_DIR" ]; then
    echo -e "${RED}错误: 未指定源目录。用法:${NC} $0 -s /path/to/test_modules"
    exit 1
fi
if [ ! -d "$SOURCE_DIR" ]; then
    echo -e "${RED}错误: 源目录不存在: $SOURCE_DIR${NC}"
    exit 1
fi

cd "$PROJECT_ROOT"

echo -e "${YELLOW}步骤 1: 检查项目结构...${NC}"
if [ ! -f "pyproject.toml" ]; then
    echo -e "${RED}错误: 不在 CircleSeeker 项目根目录${NC}"
    exit 1
fi
echo -e "${GREEN}✓ 项目目录确认${NC}"

echo -e "${YELLOW}步骤 2: 检查测试框架...${NC}"
if [ ! -f "pytest.ini" ]; then
    echo -e "${RED}错误: pytest.ini 不存在，请先运行基础框架设置${NC}"
    exit 1
fi
echo -e "${GREEN}✓ pytest 框架已安装${NC}"

echo -e "${YELLOW}步骤 3: 创建测试数据目录结构...${NC}"
mkdir -p tests/test_data
mkdir -p tests/unit
mkdir -p tests/integration

# 获取所有步骤模块
MODULES=$(find "$SOURCE_DIR" -name "step*" -type d | sort)

echo -e "${YELLOW}步骤 4: 迁移各模块测试数据...${NC}"
for module_path in $MODULES; do
    module_name=$(basename "$module_path")
    clean_name=${module_name#step??_}  # 移除 stepNN_ 前缀
    
    echo -e "  处理模块: ${BLUE}$clean_name${NC}"
    
    # 创建模块测试数据目录
    mkdir -p "tests/test_data/$clean_name"
    
    # 复制输入数据
    if [ -d "$module_path/input" ]; then
        cp -r "$module_path/input" "tests/test_data/$clean_name/"
        echo -e "    ${GREEN}✓ 复制输入数据${NC}"
    fi
    
    # 复制期望输出
    if [ -d "$module_path/expected" ]; then
        cp -r "$module_path/expected" "tests/test_data/$clean_name/"
        echo -e "    ${GREEN}✓ 复制期望输出${NC}"
    fi
    
    # 创建基础测试文件（如果不存在）
    test_file="tests/unit/test_$clean_name.py"
    if [ ! -f "$test_file" ]; then
        cat > "$test_file" << EOF
"""Unit tests for $clean_name module."""

import pytest
from pathlib import Path
import pandas as pd
import subprocess

# TODO: Import actual module when ready
# from circleseeker.modules.$clean_name import ${clean_name^}


class Test${clean_name^}:
    """Test cases for $clean_name module."""
    
    @pytest.fixture
    def test_data_dir(self):
        """Get test data directory."""
        return Path(__file__).parent.parent / "test_data" / "$clean_name"
    
    @pytest.mark.skip(reason="Module import not ready")
    def test_basic_functionality(self, test_data_dir, temp_dir):
        """Test basic $clean_name functionality."""
        # Load test input
        input_dir = test_data_dir / "input"
        expected_dir = test_data_dir / "expected"
        
        if not input_dir.exists():
            pytest.skip(f"No test data found in {input_dir}")
        
        # TODO: Implement actual test when module is ready
        assert True, "Placeholder test"
    
    def test_test_data_exists(self, test_data_dir):
        """Verify test data is available."""
        assert test_data_dir.exists(), f"Test data directory missing: {test_data_dir}"
        
        input_dir = test_data_dir / "input"
        expected_dir = test_data_dir / "expected"
        
        if input_dir.exists():
            input_files = list(input_dir.glob("*"))
            print(f"Input files found: {len(input_files)}")
            
        if expected_dir.exists():
            expected_files = list(expected_dir.glob("*"))
            print(f"Expected files found: {len(expected_files)}")
EOF
        echo -e "    ${GREEN}✓ 创建测试文件: $test_file${NC}"
    fi
done

echo -e "${YELLOW}步骤 5: 复制通用测试工具...${NC}"
if [ -f "$SOURCE_DIR/_common/compare_utils.py" ]; then
    # 内容已经集成到 tests/utils.py 中
    echo -e "${GREEN}✓ 测试工具已集成${NC}"
fi

echo -e "${YELLOW}步骤 6: 创建测试运行脚本...${NC}"
mkdir -p "$PROJECT_ROOT/scripts"
cat > "$PROJECT_ROOT/scripts/run_tests.sh" << 'EOF'
#!/bin/bash
# CircleSeeker 测试运行脚本

echo "=== CircleSeeker 测试运行器 ==="

# 检查是否安装了测试依赖
if ! python -c "import pytest" 2>/dev/null; then
    echo "安装测试依赖..."
    pip install pytest pytest-cov pytest-xdist hypothesis
fi

case "$1" in
    "unit")
        echo "运行单元测试..."
        pytest tests/unit/ -v
        ;;
    "integration") 
        echo "运行集成测试..."
        pytest tests/integration/ -v
        ;;
    "coverage")
        echo "运行测试并生成覆盖率报告..."
        pytest --cov=circleseeker --cov-report=html --cov-report=term-missing
        echo "覆盖率报告生成在: htmlcov/index.html"
        ;;
    "parallel")
        echo "并行运行测试..."
        pytest -n auto -v
        ;;
    "quick")
        echo "快速测试（跳过慢速测试）..."
        pytest -m "not slow" -v
        ;;
    *)
        echo "运行所有测试..."
        echo "用法: $0 [unit|integration|coverage|parallel|quick]"
        pytest -v
        ;;
esac
EOF

chmod +x "$PROJECT_ROOT/scripts/run_tests.sh"
echo -e "${GREEN}✓ 创建测试运行脚本: scripts/run_tests.sh${NC}"

echo ""
echo -e "${GREEN}=== 迁移完成! ===${NC}"
echo ""
echo -e "${BLUE}下一步操作:${NC}"
echo "1. 安装测试依赖:"
echo -e "   ${YELLOW}pip install pytest pytest-cov pytest-xdist hypothesis${NC}"
echo ""
echo "2. 运行测试:"
echo -e "   ${YELLOW}./scripts/run_tests.sh${NC}                    # 运行所有测试"
echo -e "   ${YELLOW}./scripts/run_tests.sh unit${NC}               # 只运行单元测试"
echo -e "   ${YELLOW}./scripts/run_tests.sh coverage${NC}           # 生成覆盖率报告"
echo ""
echo "3. 查看生成的文件:"
echo -e "   ${YELLOW}ls tests/unit/${NC}                    # 查看单元测试文件"
echo -e "   ${YELLOW}ls tests/test_data/${NC}               # 查看测试数据"
echo ""
echo -e "${BLUE}注意: 测试文件中的模块导入需要根据实际代码结构调整${NC}"
