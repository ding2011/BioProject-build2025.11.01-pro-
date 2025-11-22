# BioProject-build2025.11.01-pro-
![版本](https://img.shields.io/badge/build-2025.11.01-blue)
![Python](https://img.shields.io/badge/Python-3.8%2B-green)
![许可证](https://img.shields.io/badge/license-科研专用-orange)

一款 **工业级生物序列分析平台**，集成科研级算法与高性能架构，支持DNA/RNA序列的全方位分析，适用于科研机构、高校及生物实验室。


![图标](1.ico)

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
## 📸 软件界面
![主界面](2.png)

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
git clone https://github.com/顶011/BioProject.git
cd BioProject
2. 安装依赖
bash
# 创建虚拟环境（推荐）
python -m venv bioproject-env
# Windows激活：
bioproject-env\Scripts\activate
# macOS/Linux激活：
source bioproject-env/bin/activate

# 安装依赖包
pip install -r requirements.txt
3. 安装可选工具（增强功能）
MAFFT：
官网下载：https://mafft.cbrc.jp/alignment/software/，安装后添加到系统环境变量
RNAfold：
官网下载：https://www.tbi.univie.ac.at/RNA/RNAfold.html，安装后添加到系统环境变量
🎯 快速使用
1. 启动软件
bash
python bioproject_main.py  # 假设主脚本文件名改为 bioproject_main.py
2. 核心操作流程
输入序列：直接粘贴 FASTA 格式 / 纯序列，或通过「📥 导入 FASTA」加载文件
选择工具：
基础分析：点击「🧬 碱基统计」获取序列核心信息
ORF 预测：调整最小长度 / 遗传密码表，点击「🔍 ORF 预测」
引物设计：设置长度（18-25bp）/Tm 值（55-65℃），点击「🧪 引物设计」
查看结果：右侧面板切换「分析结果」「可视化」标签页
导出数据：支持序列、图表、分析历史导出为 FASTA/CSV/PNG 格式
3. 批量处理
点击「批量处理」分组框，添加多个 FASTA 文件
选择分析类型（如「批量碱基统计」）
点击「开始批量分析」，查看进度与汇总结果
⌨️ 常用快捷键
F5：快速分析（碱基统计 + ORF 预测 + 引物设计）
F6：单独 ORF 预测
F7：单独引物设计
Ctrl+O：导入序列
Ctrl+S：导出结果
F1：查看帮助文档
F12：关于软件
📂 目录结构
❌ 常见问题排查
1. 依赖安装失败
问题：PyQt5 安装报错 → 解决方案：pip install PyQt5==5.15.0 --no-cache-dir
问题：Biopython 版本不兼容 → 确保安装 1.79 + 版本：pip install biopython>=1.79
2. 软件启动无响应
检查 Python 版本是否≥3.8
关闭其他占用大量内存的程序（建议预留 4GB 以上内存）
3. 批量分析报错
检查序列文件格式（仅支持 FASTA 格式，无特殊字符）
长序列建议分批处理（每批≤20 个文件）
📄 许可证说明
本软件仅供 科研和学习使用，禁止商业用途。使用前请确认：
不得用于商业产品开发
引用时需注明软件名称及版本：BioProject build2025.11.01 专业版
禁止二次分发或修改后重新打包
📞 技术支持
如有问题请联系：ding20110114@outlook.com




注意说明：
引物设计功能还未完善
