# Polynomial Kernelization for Maximum Shortest-Cycle Coverage on Planar Graphs

[论文 PDF](./Polynomial-Kernelization-for-Maximum-Shortest-Cycle-Coverage-on-Planar-Graphs.pdf)

稿件研究正整数加权平面图上的最大最短环覆盖问题。给定 $k$ 和阈值 $T$，判定能否选出 $k$ 条不同的全局最短简单环，使其边并集总权重至少为 $T$。

## 复现论文第 8 节

在 Python 3.11 环境下，从仓库根目录运行以下命令。

```sh
python -m pip install -r requirements.txt
python -B verify_shortest_cycle_coverage_followup.py
python -B verify_shortest_cycle_polynomial_kernel.py
```

- `verify_shortest_cycle_coverage_followup.py` 在 120 个模型上对照穷举，检查按行计算的覆盖得分、区间动态规划状态和标记割保留的最优覆盖值。
- `verify_shortest_cycle_polynomial_kernel.py` 在 104 个模型上检查边带构造、两阶段传播、交互分量追踪和同侧重叠计算。它通过同目录的第一份脚本构造实例。

两份脚本验证论文所列有限模型及局部机制。通用核化的正确性与完整编码长度由正文证明给出。
