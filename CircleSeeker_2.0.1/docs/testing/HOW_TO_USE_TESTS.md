# 🧪 CircleSeeker 测试系统使用指南

## 📋 概述

你的CircleSeeker2.0.1项目现在拥有了一个完整的现代化pytest测试框架！所有文件都已准备就绪。

## 🎯 立即开始（5分钟快速上手）

### 第一步：安装依赖（在项目根目录执行）

```bash
cd <你的项目根目录>
# 安装测试依赖
pip install pytest pytest-cov pytest-xdist hypothesis
```

### 第二步：运行示例测试

```bash
# 运行所有可用测试
pytest -v

# 应该看到类似输出：
# ======= test session starts =======
# tests/test_carousel.py::TestCarousel::test_carousel_initialization PASSED
# tests/test_carousel.py::TestCarousel::test_process_simple_reads PASSED
# ...
```

### 第三步：迁移现有测试数据（推荐，在项目根目录执行）

```bash
# 运行迁移脚本
./scripts/migrate_existing_tests.sh

# 这会自动：
# - 复制 test_modules 中的所有测试数据
# - 为每个模块创建基础测试文件
# - 设置完整的测试结构
```

## 📁 文件结构说明

迁移完成后，你会得到以下结构：

```
CircleSeeker2.0.1/
├── tests/                          # 测试目录 ✅
│   ├── conftest.py                 # pytest 配置和公共fixtures
│   ├── utils.py                    # 测试工具函数
│   ├── test_data/                  # 测试数据 📁
│   │   ├── carousel/              # Carousel模块测试数据
│   │   │   ├── input/             # 输入文件
│   │   │   └── expected/          # 期望输出
│   │   ├── gatekeeper/            # Gatekeeper模块测试数据
│   │   └── ...                    # 其他模块
│   ├── unit/                      # 单元测试 📁
│   │   ├── test_carousel.py       # Carousel模块测试
│   │   ├── test_gatekeeper.py     # Gatekeeper模块测试
│   │   └── ...                    # 其他模块测试
│   └── integration/               # 集成测试
│       └── test_pipeline_integration.py
├── pytest.ini                     # pytest 配置文件 ✅
├── scripts/                        # 脚本目录 ✅
│   ├── migrate_existing_tests.sh   # 迁移脚本
│   └── run_tests.sh                # 测试运行器
└── TESTING.md                      # 详细测试文档 ✅
```

## 🚀 常用操作

### 运行测试的不同方式

```bash
# 基本用法
./scripts/run_tests.sh                      # 运行所有测试
./scripts/run_tests.sh unit                 # 只运行单元测试
./scripts/run_tests.sh coverage             # 生成覆盖率报告
./scripts/run_tests.sh parallel             # 并行运行测试
./scripts/run_tests.sh quick                # 跳过慢速测试

# 直接使用 pytest
pytest tests/unit/test_carousel.py  # 运行特定模块
pytest -k "carousel"                # 运行包含"carousel"的测试
pytest -m "not slow"                # 跳过标记为slow的测试
```

### 查看测试覆盖率

```bash
./scripts/run_tests.sh coverage

# 然后打开报告
open htmlcov/index.html  # macOS
```

### 调试测试问题

```bash
# 显示详细错误信息
pytest --tb=long -v

# 显示print输出
pytest -s

# 进入调试器
pytest --pdb
```

## 🔧 自定义和扩展

### 为新模块添加测试

1. **创建测试文件**：
```bash
# 例如为 Harmonizer 模块添加测试
cp tests/unit/test_carousel.py tests/unit/test_harmonizer.py
```

2. **修改测试内容**：
```python
# 编辑 test_harmonizer.py
from circleseeker.modules.harmonizer import Harmonizer

class TestHarmonizer:
    # 你的测试代码
```

3. **添加测试数据**：
```bash
mkdir -p tests/test_data/harmonizer/{input,expected}
# 添加测试文件到这些目录
```

### 添加集成测试

```python
# tests/integration/test_my_pipeline.py
@pytest.mark.integration
def test_complete_workflow(test_data_dir, temp_dir):
    """测试完整工作流程"""
    # 你的集成测试代码
```

## 💡 最佳实践

### 1. 测试驱动开发
```bash
# 先写测试
vim tests/unit/test_new_feature.py

# 然后实现功能
vim src/circleseeker/modules/new_feature.py

# 运行测试验证
pytest tests/unit/test_new_feature.py
```

### 2. 使用fixtures
```python
@pytest.fixture
def sample_data():
    return {"key": "value"}

def test_with_fixture(sample_data):
    assert sample_data["key"] == "value"
```

### 3. 参数化测试
```python
@pytest.mark.parametrize("input,expected", [
    ("ACGT", 4),
    ("ATCG", 4),
    ("", 0)
])
def test_sequence_length(input, expected):
    assert len(input) == expected
```

## 📊 监控测试质量

### 检查测试覆盖率
```bash
# 生成覆盖率报告
pytest --cov=circleseeker --cov-report=term-missing

# 目标：>80% 覆盖率
```

### 测试性能
```bash
# 查看最慢的10个测试
pytest --durations=10

# 标记慢速测试
@pytest.mark.slow
def test_large_dataset():
    # 耗时测试
```

## 🔄 CI/CD 集成

可以添加GitHub Actions：

```yaml
# .github/workflows/tests.yml
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v3
    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: '3.9'
    - name: Install dependencies
      run: pip install -e ".[test]"
    - name: Run tests
      run: pytest --cov=circleseeker
```

## ❓ 常见问题解答

### Q: 测试找不到模块？
```bash
# 确保安装了开发版本
pip install -e .
```

### Q: 外部工具测试失败？
```bash
# 跳过需要外部工具的测试
pytest -m "not external"
```

### Q: 测试太慢？
```bash
# 并行运行
pytest -n auto

# 跳过慢速测试
pytest -m "not slow"
```

### Q: 想要重新迁移？
```bash
# 删除现有测试数据
rm -rf tests/test_data tests/unit

# 重新运行迁移
./scripts/migrate_existing_tests.sh
```

## 📞 获取帮助

- 详细文档：`TESTING.md`
- 快速指南：`PYTEST_QUICKSTART.md`
- 问题反馈：在项目中创建issue

---

🎉 现在你就拥有了一个现代化、可扩展的测试框架！开始编写高质量的测试吧！
