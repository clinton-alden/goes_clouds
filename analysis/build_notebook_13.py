"""Build Notebook 13 as a comparison-only notebook."""
from pathlib import Path
import nbformat as nbf

nb = nbf.v4.new_notebook()
nb["cells"] = [
    nbf.v4.new_markdown_cell("""# 13) Mammoth/CUES binary time-series comparison

This notebook performs **comparison only**. All RGB selection, spatial averaging, ERA5-Land matching, notebook-5b TSI-rule classification, CUES cleaning, clear-sky calculation, and $k_t$ classification are completed externally by `build_mammoth_cues_binary_timeseries.py`.

Inputs are two finished binary time series:

- RGB: `0=clear`, `1=cloud`, using every available raw Mammoth RGB timestamp.
- CUES shortwave: `1=cloud` for $k_t<0.55$ and `0=clear` for $k_t>0.85$; the interval from 0.55 through 0.85 is excluded as ambiguous."""),
    nbf.v4.new_markdown_cell("## 1. Load the two independently processed binary time series"),
    nbf.v4.new_code_cell("""from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (confusion_matrix, classification_report,
                             balanced_accuracy_score, cohen_kappa_score, f1_score)

ROOT = Path.cwd()
ANALYSIS = ROOT if (ROOT / 'output_13_mammoth_domain_cloud_fraction').exists() else ROOT / 'analysis'
OUT = ANALYSIS / 'output_13_mammoth_domain_cloud_fraction'
RGB_BINARY = OUT / 'mammoth_cues_rgb_binary_all_available.csv'
SW_BINARY = OUT / 'mammoth_cues_sw_binary_all_available.csv'

rgb = pd.read_csv(RGB_BINARY, parse_dates=['time'])
sw = pd.read_csv(SW_BINARY, parse_dates=['time'])
print(f'RGB binary rows: {len(rgb):,} ({rgb.time.min()} to {rgb.time.max()})')
print(f'SW binary rows:  {len(sw):,} ({sw.time.min()} to {sw.time.max()})')
display(rgb.head(), sw.head())"""),
    nbf.v4.new_markdown_cell("""## 2. Match the binary time series and restrict to 17Z–23Z

The series are matched by exact 5-minute timestamp. Only UTC hours 17 through 23 inclusive are retained to remove lower-sun-angle scenes."""),
    nbf.v4.new_code_cell("""compare = rgb[['time','rgb_cloud_binary']].merge(
    sw[['time','sw_cloud_binary','k_t','cos_sza']], on='time', how='inner', validate='one_to_one'
)
compare = compare.loc[compare.time.dt.hour.between(17, 23)].reset_index(drop=True)
compare['agree'] = compare.rgb_cloud_binary.eq(compare.sw_cloud_binary)
print(f'Exact matched timestamps in 17Z–23Z window: {len(compare):,}')
print(f'Agreement: {compare.agree.mean():.3%}')
display(compare.head())"""),
    nbf.v4.new_markdown_cell("## 3. Confusion matrix and binary skill metrics"),
    nbf.v4.new_code_cell("""truth, prediction = compare.sw_cloud_binary, compare.rgb_cloud_binary
cm = confusion_matrix(truth, prediction, labels=[0,1])
report = classification_report(truth, prediction, labels=[0,1], target_names=['clear','cloud'], output_dict=True, zero_division=0)
metrics = pd.DataFrame([{
    'n': len(compare), 'accuracy': (truth == prediction).mean(),
    'balanced_accuracy': balanced_accuracy_score(truth, prediction),
    'cohen_kappa': cohen_kappa_score(truth, prediction),
    'macro_f1': f1_score(truth, prediction, average='macro', zero_division=0),
    'weighted_f1': f1_score(truth, prediction, average='weighted', zero_division=0),
    'clear_precision': report['clear']['precision'], 'clear_recall': report['clear']['recall'],
    'clear_f1': report['clear']['f1-score'], 'cloud_precision': report['cloud']['precision'],
    'cloud_recall': report['cloud']['recall'], 'cloud_f1': report['cloud']['f1-score'],
}])
display(metrics.T)
fig, ax = plt.subplots(figsize=(5,4))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=False,
            xticklabels=['RGB clear','RGB cloud'], yticklabels=['SW clear','SW cloud'], ax=ax)
ax.set_xlabel('Updated-rule RGB binary'); ax.set_ylabel('$k_t$ shortwave binary')
ax.set_title('Mammoth/CUES binary comparison')
fig.tight_layout(); fig.savefig(OUT/'mammoth_cues_binary_confusion_matrix.png', dpi=200)
plt.show()"""),
    nbf.v4.new_markdown_cell("""## 4. Confusion matrices by UTC hour

Only hours containing unambiguous matched observations are shown. “Overall F1” in each title is macro F1 across the clear and cloud classes for that hour."""),
    nbf.v4.new_code_cell("""import numpy as np

compare['utc_hour'] = compare.time.dt.hour
hours = sorted(compare.utc_hour.unique())
ncols = 4
nrows = int(np.ceil(len(hours) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3.5*nrows))
axes = np.atleast_1d(axes).ravel()
hourly_metrics = []
for ax, hour in zip(axes, hours):
    group = compare.loc[compare.utc_hour == hour]
    group_cm = confusion_matrix(group.sw_cloud_binary, group.rgb_cloud_binary, labels=[0,1])
    group_f1 = f1_score(group.sw_cloud_binary, group.rgb_cloud_binary, average='macro', zero_division=0)
    hourly_metrics.append({'utc_hour': hour, 'n': len(group), 'macro_f1': group_f1})
    sns.heatmap(group_cm, annot=True, fmt='d', cmap='Blues', cbar=False, ax=ax,
                xticklabels=['clear','cloud'], yticklabels=['clear','cloud'])
    ax.set_title(f'{hour:02d} UTC | overall F1={group_f1:.3f} | n={len(group):,}')
    ax.set_xlabel('RGB'); ax.set_ylabel('$k_t$')
for ax in axes[len(hours):]: ax.axis('off')
fig.suptitle('Mammoth/CUES confusion matrices by UTC hour', fontsize=15, y=1.01)
fig.tight_layout(); fig.savefig(OUT/'mammoth_cues_confusion_matrices_by_utc_hour.png', dpi=200, bbox_inches='tight')
plt.show()
hourly_metrics = pd.DataFrame(hourly_metrics)
display(hourly_metrics)"""),
    nbf.v4.new_markdown_cell("""## 5. Confusion matrices by month

Each subplot uses all years for that calendar month. “Overall F1” is the month-specific macro F1 across clear and cloud."""),
    nbf.v4.new_code_cell("""compare['month'] = compare.time.dt.month
month_names = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
fig, axes = plt.subplots(3, 4, figsize=(16, 11))
monthly_metrics = []
for month, ax in enumerate(axes.ravel(), start=1):
    group = compare.loc[compare.month == month]
    group_cm = confusion_matrix(group.sw_cloud_binary, group.rgb_cloud_binary, labels=[0,1])
    group_f1 = f1_score(group.sw_cloud_binary, group.rgb_cloud_binary, average='macro', zero_division=0)
    monthly_metrics.append({'month': month, 'month_name': month_names[month-1], 'n': len(group), 'macro_f1': group_f1})
    sns.heatmap(group_cm, annot=True, fmt='d', cmap='Blues', cbar=False, ax=ax,
                xticklabels=['clear','cloud'], yticklabels=['clear','cloud'])
    ax.set_title(f'{month_names[month-1]} | overall F1={group_f1:.3f} | n={len(group):,}')
    ax.set_xlabel('RGB'); ax.set_ylabel('$k_t$')
fig.suptitle('Mammoth/CUES confusion matrices by calendar month', fontsize=15, y=1.01)
fig.tight_layout(); fig.savefig(OUT/'mammoth_cues_confusion_matrices_by_month.png', dpi=200, bbox_inches='tight')
plt.show()
monthly_metrics = pd.DataFrame(monthly_metrics)
display(monthly_metrics)"""),
    nbf.v4.new_markdown_cell("## 6. Diagnostic comparison time series"),
    nbf.v4.new_code_cell("""week_start = pd.Timestamp('2024-07-01 15:00')
week = compare.loc[compare.time.between(week_start, week_start + pd.Timedelta(days=7))]
fig, ax = plt.subplots(figsize=(14,4))
ax.step(week.time, week.sw_cloud_binary, where='mid', label='$k_t$ shortwave binary', lw=1.6)
ax.step(week.time, week.rgb_cloud_binary, where='mid', label='updated-rule RGB binary', lw=1.1, alpha=.8)
ax.set_yticks([0,1], ['clear','cloud']); ax.set_xlabel('UTC'); ax.grid(alpha=.25); ax.legend()
fig.tight_layout(); fig.savefig(OUT/'mammoth_cues_binary_week_20240701.png', dpi=200)
plt.show()"""),
    nbf.v4.new_markdown_cell("## 7. Export the matched comparison and metrics"),
    nbf.v4.new_code_cell("""compare.to_csv(OUT/'mammoth_cues_binary_comparison_all_available.csv', index=False)
metrics.to_csv(OUT/'mammoth_cues_binary_metrics_all_available.csv', index=False)
hourly_metrics.to_csv(OUT/'mammoth_cues_binary_metrics_by_utc_hour.csv', index=False)
monthly_metrics.to_csv(OUT/'mammoth_cues_binary_metrics_by_month.csv', index=False)
print('Saved comparison outputs to', OUT.resolve())"""),
]
nb["metadata"]["kernelspec"] = {"display_name":"Python 3","language":"python","name":"python3"}
nb["metadata"]["language_info"] = {"name":"python","version":"3"}
nbf.write(nb, Path("13_mammoth_domain_cloud_fraction_vs_sw_ratio.ipynb"))
