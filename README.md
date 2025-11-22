# BioProject-build2025.11.01-pro-
# BioProject 专业版
![版本](https://img.shields.io/badge/build-2025.11.01-blue)
![Python](https://img.shields.io/badge/Python-3.8%2B-green)
![许可证](https://img.shields.io/badge/license-科研专用-orange)

一款 **工业级生物序列分析平台**，集成科研级算法与高性能架构，支持DNA/RNA序列的全方位分析，适用于科研机构、高校及生物实验室。

## 🚀 核心功能
基于 build2025.11.01 版本，包含以下工业级特性：
- **序列基础分析**：长度、GC含量、分子量、碱基组成、序列复杂度计算
- **科研级核心工具**：
  - ORF预测（支持原核/真核生物、多遗传密码表）
  - 引物设计（精准Tm值计算、二级结构风险评估）
  - Smith-Waterman局部比对（长序列启发式优化）
  - 多序列比对（MAFFT包装器，支持错误恢复）
  - 进化树构建（NJ法+UPGMA回退机制）
- **批量处理**：多文件并行分析、进度跟踪、错误恢复
- **高级可视化**：GC含量趋势图、ORF分布图、引物特性散点图
- **数据持久化**：分析历史存储、序列缓存、批量任务断点续传
- **工业级体验**：多线程优化、内存监控、自适应界面、快捷键支持

## 📋 系统要求
### 基础环境
- Python 3.8 及以上
- 推荐配置：多核CPU + 8GB+ RAM（批量处理长序列需16GB+）
- 操作系统：Windows 10+/macOS 12+/Linux（Ubuntu 20.04+/CentOS 8+）

### 必需依赖包
见 `requirements.txt`，包含 Biopython、PyQt5、NumPy 等核心库

### 可选工具（增强功能）
- MAFFT：多序列比对核心工具（需添加到系统PATH）
- RNAfold：二级结构检测（用于引物风险评估）

## 🔧 安装步骤
### 1. 克隆仓库
```bash
git clone https://github.com/你的用户名/BioProject.git
cd BioProject
