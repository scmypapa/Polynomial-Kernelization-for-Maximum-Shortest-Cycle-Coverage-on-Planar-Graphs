# Polynomial Kernelization for Maximum Shortest-Cycle Coverage on Planar Graphs

[论文 PDF](./Polynomial-Kernelization-for-Maximum-Shortest-Cycle-Coverage-on-Planar-Graphs.pdf) · [论文 LaTeX 源文件](./manuscript.tex) · [审稿修订版归档](https://github.com/scmypapa/Polynomial-Kernelization-for-Maximum-Shortest-Cycle-Coverage-on-Planar-Graphs/tree/reproducibility-v2)

稿件研究正整数加权平面图上的最大最短环覆盖问题。给定 $k$ 和阈值 $T$，判定能否选出 $k$ 条不同的全局最短简单环，使其边并集总权重至少为 $T$。

## 复现论文第 8 节和附录 A

在 Python 3.11 环境中，从仓库根目录运行以下命令。

    python -m pip install -r requirements.txt
    python -B verify_shortest_cycle_coverage_followup.py
    python -B verify_shortest_cycle_polynomial_kernel.py

- 第一份脚本在 120 个模型上对照穷举，检查按行计算的覆盖得分、区间动态规划状态和标记割保留的最优覆盖值。
- 第二份脚本在 107 个模型上检查边带、两阶段传播、分量追踪和同侧重叠。定向模型包含三个自配对区域，以及同一原区间选中多个块的情形。脚本核对了 12,288 次替换后的固定记账，其中 17 次为严格上界。

附录 A 说明测试流程、计数口径和适用范围。主定理的通用核化与位长界限由正文证明给出。

## 重新生成论文和图示

图 1–3 的 TikZ 源文件在 [figures](./figures) 中。安装 Tectonic 后，在仓库根目录运行以下命令生成论文 PDF。

    tectonic manuscript.tex
