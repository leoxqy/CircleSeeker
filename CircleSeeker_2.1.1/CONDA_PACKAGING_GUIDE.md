# CircleSeeker Conda打包指南

## 📦 快速开始

### 1. 本地构建测试

```bash
# 安装必要工具
conda install -c conda-forge conda-build conda-verify anaconda-client

# 方法一：使用现有脚本（推荐）
cd conda-recipe
bash test_local.sh

# 方法二：手动构建
conda build conda-recipe -c conda-forge -c bioconda --python=3.9
```

### 2. 安装测试

```bash
# 创建测试环境
conda create -n circleseeker_test python=3.9
conda activate circleseeker_test

# 安装本地构建的包
conda install /path/to/built/circleseeker-2.1.1-py39_0.tar.bz2 -c conda-forge -c bioconda

# 验证安装
circleseeker --help
python -c "import circleseeker; print(circleseeker.__version__)"
```

## 🚀 发布流程

### 选项A：发布到Anaconda Cloud（个人频道）

```bash
# 登录Anaconda
anaconda login

# 上传包
anaconda upload /path/to/built/circleseeker-2.1.1-py39_0.tar.bz2

# 用户安装命令
conda install -c your-username -c conda-forge -c bioconda circleseeker
```

### 选项B：发布到Bioconda（推荐）

1. **Fork bioconda-recipes仓库**
   ```bash
   git clone https://github.com/YOUR_USERNAME/bioconda-recipes.git
   cd bioconda-recipes
   ```

2. **创建新分支**
   ```bash
   git checkout -b add-circleseeker
   ```

3. **添加recipe**
   ```bash
   mkdir -p recipes/circleseeker
   cp /path/to/your/conda-recipe/meta.yaml recipes/circleseeker/
   ```

4. **修改meta.yaml的source部分**
   ```yaml
   source:
     url: https://github.com/circleseeker/circleseeker2/archive/refs/tags/v{{ version }}.tar.gz
     sha256: YOUR_SHA256_HASH  # 使用 sha256sum 计算
   ```

5. **测试recipe**
   ```bash
   # 使用bioconda-utils测试
   pip install bioconda-utils
   bioconda-utils lint --packages circleseeker
   bioconda-utils build --packages circleseeker
   ```

6. **提交PR**
   ```bash
   git add recipes/circleseeker
   git commit -m "Add circleseeker 2.1.1"
   git push origin add-circleseeker
   ```

   然后在GitHub上创建Pull Request到bioconda/bioconda-recipes

## 📝 检查清单

### 构建前检查
- [x] 版本号一致（pyproject.toml、meta.yaml、__version__.py）
- [x] 所有依赖在conda-forge或bioconda中可用
- [x] 删除所有临时文件和缓存
- [x] 更新.gitignore

### 测试项目
- [ ] 命令行工具能正常运行
- [ ] Python包能正常导入
- [ ] 所有外部工具依赖能正确安装
- [ ] 在Linux和macOS上测试通过

## 🔧 常见问题

### 1. 构建失败：找不到依赖
```bash
# 确保添加了正确的channels
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
```

### 2. cyrcular包问题
cyrcular已在bioconda中可用，无需从GitHub安装：
```bash
conda install -c bioconda cyrcular
```

### 3. 版本更新
更新版本时，需要同时修改：
- `pyproject.toml` 中的 version
- `conda-recipe/meta.yaml` 中的 version
- `src/circleseeker/__version__.py` 中的 __version__

## 📚 参考资源

- [Bioconda贡献指南](https://bioconda.github.io/contributor/index.html)
- [Conda-build文档](https://docs.conda.io/projects/conda-build/en/latest/)
- [CircleSeeker GitHub](https://github.com/circleseeker/circleseeker2)

## ⚡ 快速命令汇总

```bash
# 本地构建
cd conda-recipe && bash test_local.sh

# 手动构建
conda build conda-recipe -c conda-forge -c bioconda

# 查看构建输出路径
conda build conda-recipe --output

# 清理构建缓存
conda build purge
```

---
生成时间：2025-09-25
CircleSeeker版本：2.1.1