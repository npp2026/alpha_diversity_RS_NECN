# v1.6.9 方法口径更新

## Submission cleanup

本次交付在不改变已确认分析算法的前提下完成投稿清理：补充了 MS/SI 对照记录，统一了运行依赖清单，移除未使用的 `randomForest` 核心依赖，并保留固定15变量、原始 RF OOF、年度 OOB-QM 和 QRF 老龄验证设计的代码约束。投稿附件进一步精简：删除根目录 `tests/`、`validation/`、`environment/` 和 `docs/`，将输入输出、运行顺序、方法映射和 MS/SI 对照合并到根 `README.md`。上游断点引擎内部的 `validation/` 目录及其合成数据仍保留。方法对照记录还标出了 SI Table S1 中需要在文稿侧修正的“一次筛选”表述。

本次按作者确认更新上传的 `FEM_Supplementary_Code_v1.6.9_analysis_only.zip`，版本号仍为 1.6.9，交付文件名增加 `updated`，用于与原包区分。

| 环节 | 更新后的口径 |
|---|---|
| 生产模型 | 两个响应继续使用原有固定15变量；变量顺序和 RF 超参数保留。 |
| VSURF | 仍用于模块02的折内验证筛选；不读取为生产输入，不自动替换生产变量清单。 |
| Figure2及来源统计 | 场景 R²、来源移除、边际贡献、替代效应及 Bootstrap 使用原始 RF 外层测试折 OOF 预测。 |
| 年度制图 | 保留训练 OOB 拟合的 QM 校正；年度预测脚本及共享 QM 实现保留。 |

## 实际修改

- `run_nested_MS.R` 的场景对象改为 `predictions=rf`，保留 `raw_predictions`，增加 `prediction_type="RF_raw"`。QM 存入单独的 `calibrated_predictions`，仍保留 RF_raw/RF_QM/XGBoost_raw 的独立诊断导出。
- `R/oof_contracts.R` 增加原始 RF 选择入口；三个后处理脚本均先明确选择原始 RF，再检查样本配对和计算统计量。没有原始预测时停止，不回退到 QM。
- 统计表增加 `Prediction_Type=RF_raw`；原输出文件名、R²/贡献公式、共享样本及折、默认50 km/2000次 Bootstrap、BH 分组规则保留。
- 生产训练增加固定15变量数量检查和流程分离元数据；固定变量清单与训练计算保留。删除说明中尚待作者确认生产名单的旧表述。
- 同步根 README、模块内链接、完整性检查器、文件清单和 SHA256；根测试记录、发布验证记录和未锁定环境模板不再随投稿附件分发。

## 已有结果如何重算

原 v1.6.9 场景对象已包含 `raw_predictions` 时，可以直接重算后处理，无需为了切换预测口径重新训练模型。将 `all_results.rds` 复制到独立目录，再从更新包根目录执行：

```bash
mkdir -p /absolute/results/raw_rf_postprocess
cp /absolute/results/old_nested/all_results.rds /absolute/results/raw_rf_postprocess/
bash 02_model_scenario_validation/run_postprocess_submission.sh /absolute/results/raw_rf_postprocess
```

旧 RDS 不被写回；原 QM 统计表应单独保留。只有 QM 或没有明确标识的预测不能转换为原始 OOF，需取得原始预测或重新运行模块02。仅含原始 `predictions` 的对象需明确设置 `prediction_type="RF_raw"`；只有确认其确实来自原始 RF OOF 时才能添加该标识。

## 验证范围

当前投稿包保留 `scripts/validate_release.py` 的文件/链接/语法完整性检查、`scripts/check_environment.R` 的依赖检查，以及 `scripts/run_MS_validation.R` 调用的模块级 R self-check。本环境没有 R/Rscript，也没有原始研究数据，因此 R 语法解析、R self-check、模型训练/年度制图及论文数值复现未运行。Python/Bash/JavaScript 检查与文件一致性检查不能代替 R 执行验证。

在具备 R 的环境中运行：

```bash
python3 scripts/validate_release.py --require-r
Rscript scripts/check_environment.R
Rscript scripts/run_MS_validation.R . ../FEM_validation_v169_updated
```
