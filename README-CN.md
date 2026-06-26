# Mini Engineering Transport Phenomena（迷你工程输运现象）

**从零开始、零依赖的 C 语言实现**，涵盖大学级输运现象——动量、热量和质量传递的统一框架。每个模块对应 MIT（及其他顶尖大学）的一门或多门课程，将教科书中的公式转化为可运行的 C 代码，实现理论与实践的桥接。

## 子模块

| 子模块 | 主题 | 参考课程 |
|--------|--------|-------------|
| [mini-computational-transport-cfd-basics](mini-computational-transport-cfd-basics/) | 有限差分/有限体积离散、迎风/幂律/混合格式、对流-扩散求解器、SIMPLE 算法、顶盖驱动方腔流 | MIT 2.29, MIT 16.90 |
| [mini-convective-diffusion-reaction](mini-convective-diffusion-reaction/) | Peclet/Damköhler/Thiele 数、Fick 定律、Stefan-Maxwell 扩散、催化剂效率因子、Weisz-Prater 判据 | MIT 10.50, MIT 10.37 |
| [mini-dimensionless-numbers-re-nu-pr-sc](mini-dimensionless-numbers-re-nu-pr-sc/) | Buckingham Pi 定理、Re/Pr/Sc/Nu 等无量纲数、典型几何体对流传热关联式 | MIT 10.50, MIT 2.51 |
| [mini-interphase-transport-boundary-condition](mini-interphase-transport-boundary-condition/) | 无滑移/滑移边界条件、相间跳跃条件（质量/动量/能量）、传热传质边界条件、表面张力 | MIT 10.50, Stanford ME 351 |
| [mini-momentum-heat-mass-analogy](mini-momentum-heat-mass-analogy/) | Reynolds 类比（St = f/2）、Chilton-Colburn 类比、速度/温度/浓度边界层、无量纲群耦合 | MIT 10.50, MIT 2.51 |
| [mini-transport-in-micro-nano-channel](mini-transport-in-micro-nano-channel/) | Knudsen 流区（连续到自由分子）、滑移模型、电渗、电泳、流动电势、双电层模型 | MIT 2.25, Stanford ME 457 |

## 设计理念

- **零外部依赖** — 纯 C（C99/C11），仅使用 `libc` 和 `libm`
- **模块自包含** — 每个目录自带 `Makefile`、`include/`、`src/`、`examples/`、`demos/`、`tests/`
- **理论到代码的映射** — 每个模块包含 `docs/` 目录，内有课程对齐说明
- **实用演示程序** — SIMPLE 求解器、效率因子计算器、Nu 关联式浏览器、跳跃条件验证器、类比绘图、Knudsen 流区分类器

## 构建方式

每个模块相互独立。进入模块目录后运行：

```bash
cd mini-computational-transport-cfd-basics
make all    # 构建全部
make test   # 运行测试
```

需要 **GCC** 和 **GNU Make**。

## 项目结构

```
mini-eng-transport-phenomena/
├── mini-computational-transport-cfd-basics/       # 有限差分/体积法、对流-扩散求解器、SIMPLE
├── mini-convective-diffusion-reaction/            # CDR 输运、Fick 定律、催化剂效率
├── mini-dimensionless-numbers-re-nu-pr-sc/        # Buckingham Pi、Re/Nu/Pr/Sc、Nu 关联式
├── mini-interphase-transport-boundary-condition/  # 边界条件、相间跳跃平衡
├── mini-momentum-heat-mass-analogy/               # Reynolds/Chilton-Colburn 类比、边界层理论
└── mini-transport-in-micro-nano-channel/          # Knudsen 流区、电动现象、EDL 模型
```

## 许可证

MIT
