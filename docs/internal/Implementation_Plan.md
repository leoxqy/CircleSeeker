# CircleSeeker 项目改进计划（更可落地版本）

> **当前版本**: 0.10.0  
> **目标版本**: 0.10.0（本轮里程碑）  
> **创建日期**: 2026-01-13  
> **最后更新**: 2026-01-13

## 概述

本计划聚焦 4 个方向：类型检查、测试、版本一致性、文档站点。相比“一次到位”的严格策略，本版本更偏向**渐进落地**：先让 CI 可稳定通过，再逐步提高质量门槛。

**执行原则**：
1. 先稳定可用（CI/测试/构建），再提高门槛（strict typing / 大迁移）。
2. mypy 先覆盖“核心路径”，全仓清零作为里程碑而非前置门槛。
3. 集成测试尽量不依赖外部二进制；需要外部工具的 e2e 以 `slow/optional` 形式单独跑。
4. 文档先“能构建”，再做目录重构与迁移。

**建议执行顺序**：
1. 集成测试补齐（快速回归基础）
2. mypy 渐进式（核心模块严格）
3. Sphinx 文档站点（先可构建）
4. 版本发布（根据变更决定 0.9.x 或 0.10.0）

**每步最小化验证命令**（建议每完成一个小任务就跑一次）：

```bash
# 1) 回归测试（优先跑已有 integration，快且能防回归）
pytest tests/integration/ -v -m integration

# 2) 类型门禁（只检查核心范围）
mypy src/circleseeker/cli src/circleseeker/config.py src/circleseeker/core --show-error-codes

# 3) 文档构建（需要已安装 docs 依赖）
pip install -e ".[docs]"
cd docs && make html
```

---

## 任务 1：渐进式引入 mypy（先核心、后全仓）

### 1.0 小任务拆分（建议逐步完成）

- [x] 1A. 更新 `pyproject.toml` 的 mypy 基础规则（阶段 A）
- [x] 1B. 明确“门禁范围”与执行命令（先 `cli/config/core`）
- [x] 1C. 修复核心范围 mypy → 0 errors（只要求核心门禁通过）
- [x] 1D. （可选）用 `[[tool.mypy.overrides]]` 对核心模块加严
- [x] 1E. 分批把 `utils/`、`modules/` 纳入门禁范围

### 1.1 当前状态

- **配置位置**: `pyproject.toml` 的 `[tool.mypy]`
- **现状**: 配置较宽松；全仓类型错误尚未清零（以 `mypy` 输出为准）
- **目标**: 不阻塞日常开发的前提下，把核心链路（CLI/config/core）先做到高置信度

### 1.2 分阶段配置建议

**阶段 A：全仓启用“低风险高收益”规则（不强制全量注解）**

```toml
[tool.mypy]
python_version = "3.9"
ignore_missing_imports = true
warn_return_any = true
warn_unused_configs = true

# 渐进式开启的规则（通常改动小但收益大）
check_untyped_defs = true
no_implicit_optional = true
warn_unreachable = true
warn_unused_ignores = true
warn_redundant_casts = true

# 暂不强制所有函数都写类型
disallow_untyped_defs = false
disallow_incomplete_defs = false
disallow_untyped_calls = false
disallow_any_generics = false
```

**阶段 B：对核心模块执行“严格门槛”（推荐用 CI 命令控制范围）**

- 先在 CI 中只对以下路径运行 mypy（可视作“类型门禁”）：
  - `src/circleseeker/cli/`
  - `src/circleseeker/config.py`
  - `src/circleseeker/core/`
- 等核心稳定后，再逐步把 `src/circleseeker/utils/`、`src/circleseeker/modules/` 纳入门禁范围。

如确实希望在 `pyproject.toml` 中对核心模块更严格，可用 overrides（注意：若全仓运行 mypy，非核心模块仍会报错，因此更推荐“按目录跑”）：

```toml
[[tool.mypy.overrides]]
module = [
  "circleseeker.cli.*",
  "circleseeker.config",
  "circleseeker.core.*",
]
disallow_untyped_defs = true
disallow_incomplete_defs = true

# 可进一步加严（建议后续阶段再开）：
# disallow_untyped_calls = true
# disallow_any_generics = true
```

### 1.3 优先修复范围（建议）

| 优先级 | 范围 | 原因 |
|---|---|---|
| P0 | `circleseeker/cli/*`, `circleseeker/config.py`, `circleseeker/core/*` | 用户入口 + 流程编排，类型收益最大 |
| P1 | `circleseeker/utils/dependency_checker.py` 等“早失败”模块 | 能更早暴露环境/依赖问题 |
| P2 | `circleseeker/modules/*` | 业务逻辑复杂，适合后续逐步清理 |

### 1.4 原错误热点（作为 backlog 参考）

> 说明：下表是计划创建时的参考统计；请以 `mypy` 当前输出为准。

| 文件路径 | 参考错误数（计划创建时） | 主要问题类型 |
|----------|---------------------------|--------------|
| `src/circleseeker/utils/dependency_checker.py` | 13 | 类型不兼容、参数类型 |
| `src/circleseeker/modules/ecc_unify.py` | 11 | 索引操作、对象不可索引 |
| `src/circleseeker/modules/ecc_summary.py` | 5 | 变量类型注解缺失 |
| `src/circleseeker/modules/cecc_build.py` | 4 | Any 返回类型 |
| `src/circleseeker/modules/um_classify.py` | 4 | Any 返回类型 |
| `src/circleseeker/modules/ecc_dedup.py` | 3 | Any 返回类型 |
| `src/circleseeker/modules/read_filter.py` | 3 | 可选类型处理 |
| `src/circleseeker/external/minimap2_align.py` | 2 | 方法签名不兼容 |
| `src/circleseeker/external/cd_hit.py` | 2 | 方法签名不兼容 |
| 其他文件 | ~26 | 混合问题 |

### 1.5 修复策略（更可落地）

1. **按需补齐类型存根**
   - 例如：`types-PyYAML`、`pandas-stubs`（项目当前未使用 `requests`，无需 `types-requests`）
2. **核心模块先写清“输入/输出”类型**
   - CLI option → `PipelineOptions` → pipeline 执行入口（优先保证主链路的类型一致）
3. **对动态结构使用 TypedDict/Protocol/类型守卫**
   - 避免用大量 `Any` 直接“压平”错误
4. **阶段性使用 `cast()`/类型窄化**
   - 作为迁移手段，不追求一步到位完美类型

### 1.6 验证命令（分两类）

```bash
# 作为“门禁”的核心路径检查（建议在 CI 跑）
mypy src/circleseeker/cli src/circleseeker/config.py src/circleseeker/core --show-error-codes

# 作为“指标/报表”的全仓扫描（本地/手工）
mypy src/circleseeker --show-error-codes
```

---

## 任务 2：补充集成测试（以回归价值为导向）

### 2.0 小任务拆分（建议逐步完成）

- [x] 2A. 将 `test_cresil_integration.py` 逐步收敛为纯 pytest 测试（去脚本式 main/print）
- [x] 2B. 新增配置优先级测试（CLI > config > default）
- [x] 2C. 新增 checkpoint/resume/force 行为测试
- [x] 2D. 新增 dry-run 输出/步骤结构测试（不跑外部工具）
- [x] 2E. （可选）补充 `@pytest.mark.slow` 的端到端测试并对外部工具做 skip

### 2.1 当前状态

- 已有：`tests/integration/test_cli_integration.py`（覆盖版本、help、debug gate、生成配置、show-steps、init-config 等）
- 已有：`tests/integration/test_cresil_integration.py`（更偏功能脚本式测试）

### 2.2 目标结构（精简 + 可维护）

```
tests/integration/
├── test_cli_integration.py              # 已有
├── test_config_priority.py              # 新增：CLI/config/default 优先级
├── test_checkpoint_resume.py            # 新增：checkpoint/resume/force 行为
├── test_pipeline_dry_run.py             # 新增：dry-run 输出结构/步骤列表
└── test_cresil_integration.py           # 已有（建议逐步改成纯 pytest 风格、去掉 print/main）
```

### 2.3 新增用例重点（建议）

- 配置优先级：CLI 覆盖 config、缺省值回退、缺参数报错路径
- checkpoint：创建/读取/`--force` 清理、异常后 resume 行为
- dry-run：不执行外部工具，仍能产出预期的“步骤/计划”输出与目录结构
- `validate`/`show-checkpoint`：确保 CLI 子命令行为稳定（尤其是 `--debug` gate）

### 2.4 外部工具依赖策略

- 默认集成测试不依赖 minimap2/cd-hit 等外部二进制：优先 mock 或 `skipif`。
- 真正端到端（跑全 pipeline）的测试用 `@pytest.mark.slow`（CI 默认不跑，release 前手动跑）。

### 2.5 运行命令

```bash
pytest tests/integration/ -v -m integration
pytest tests/ -v
pytest tests/integration/ -v -m "integration and slow"
```

---

## 任务 3：版本号与一致性（根据变更决定发布节奏）

### 3.0 小任务拆分（建议逐步完成）

- [x] 3A. 明确本轮发布策略：直接发布 `0.10.0`（里程碑）
- [x] 3B. 新增版本一致性单元测试（兼容 Python 3.9+）
- [x] 3C. 同步更新版本号（`__version__.py`/`pyproject.toml`/文档）并打 tag（可选）

### 3.1 版本策略建议

- 本轮发布：`0.10.0`（里程碑）
- 若本轮主要是 **内部质量改进**（类型/测试/文档），且 CLI/输出不变：优先发布 `0.9.16` 这类 patch。
- 若有 **用户可感知变更**（CLI 行为/参数、输出结构、默认流程变化等）：再升 `0.10.0`。

### 3.2 需修改的文件（保持一致）

- `src/circleseeker/__version__.py`
- `pyproject.toml`
- 文档页（如 CLI reference 顶部版本号）

### 3.3 版本一致性检查（示例）

> 注意：项目最低支持 Python 3.9；为了避免引入 TOML 解析依赖，可以只解析 `[project]` 段落的 `version` 行。

```python
import re
from pathlib import Path

def read_project_version(pyproject_path: Path) -> str:
    in_project = False
    for raw_line in pyproject_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line == "[project]":
            in_project = True
            continue
        if in_project and line.startswith("[") and line.endswith("]"):
            break
        if in_project:
            m = re.match(r'^version\s*=\s*"([^"]+)"\s*$', line)
            if m:
                return m.group(1)
    raise AssertionError("Missing [project].version")

def test_version_consistency():
    from circleseeker.__version__ import __version__
    assert __version__ == read_project_version(Path("pyproject.toml"))
```

### 3.4 Git 标签（可选）

保留原流程（按最终决定的版本号打 tag）。

---

## 任务 4：Sphinx 文档站点（先能构建，再重构迁移）

### 4.0 小任务拆分（建议逐步完成）

- [x] 4A. 补齐 Markdown 构建依赖（`myst-parser`，放在 `.[docs]`）
- [x] 4B. 增加 Sphinx 最小骨架（`docs/source/conf.py`/`index.rst`/`Makefile`）
- [x] 4C. 跑通 `cd docs && make html`（在干净环境可复现）
- [x] 4D. （可选）第二阶段再做文档迁移与信息架构重排

### 4.1 当前状态

- `docs/` 已有多份 Markdown 文档
- `pyproject.toml` 已有 docs 依赖（Sphinx/RTD theme/sphinx-click），但若要直接渲染 Markdown，建议加入 `myst-parser`

### 4.2 阶段化落地

**阶段 A：最小可用（MVP）**
- 新增 `docs/source/conf.py`、`docs/source/index.rst`、`docs/Makefile`
- 通过 `myst-parser` 直接把现有 Markdown 纳入构建（不先移动文件）
- 目标：`cd docs && make html` 稳定通过

**阶段 B：信息架构与迁移（可选，等 MVP 稳定后再做）**
- 再把文档搬到 `docs/source/user_guide/`、`developer_guide/` 等目录
- 再补齐 API autodoc、CLI 文档（sphinx-click）等

### 4.3 文档依赖建议

- 推荐优先维护在 `pyproject.toml` 的 `project.optional-dependencies.docs`（单一来源），并补上：
  - `myst-parser>=0.18.0`

### 4.4 构建命令

```bash
# 安装文档依赖（建议：使用 pyproject extras，避免重复维护 requirements.txt）
pip install -e ".[docs]"

cd docs && make html
open docs/_build/html/index.html
```

---

## 验证清单（调整后的验收口径）

- [x] `pytest tests/ -v` → all passed（允许 external toolchain 相关用例 skip）
- [x] `pytest tests/integration/ -v -m integration` → all passed
- [x] `mypy src/circleseeker/cli src/circleseeker/config.py src/circleseeker/core` → 0 errors
- [ ] `CircleSeeker --version` → 输出与 `pyproject.toml` 一致
- [x] `cd docs && make html` → Build successful

---

## 风险与回滚（更新版）

| 风险 | 等级 | 缓解/回滚 |
|------|------|-----------|
| mypy 规则过严导致大面积返工 | 中 | 采用“核心门禁 + 分阶段扩大范围”，避免全仓 strict 一步到位 |
| 集成测试依赖外部工具导致 CI 不稳定 | 中 | 默认 mock/skip；e2e 标记为 slow/optional |
| 文档大迁移引起链接断裂 | 低 | 先 MVP 可构建，迁移放到第二阶段并做索引/重定向 |

---

## 附录：关键文件路径

```
src/circleseeker/
├── __version__.py                   # 版本号
├── cli/
│   ├── __init__.py
│   ├── main.py                      # CLI 主入口
│   ├── pipeline.py                  # Pipeline 执行逻辑
│   └── commands/
│       ├── run.py
│       ├── checkpoint.py
│       ├── config.py
│       └── validate.py
├── config.py                        # 配置类
├── core/
│   └── pipeline.py                  # Pipeline 核心
└── utils/
    └── dependency_checker.py        # 依赖检查（建议优先纳入 mypy 门禁）

tests/
├── conftest.py                      # 全局 fixtures
├── integration/
│   ├── test_cli_integration.py      # 已有
│   └── test_cresil_integration.py   # 已有
└── unit/
    └── ...

docs/
├── Implementation_Plan.md           # 本文档
├── CLI_Reference.md                 # 现有文档
├── Pipeline_Modules.md              # 现有文档
└── images/logo.png                  # Logo

pyproject.toml                       # 项目配置（mypy, 版本号, 依赖）
```
