# Molecular SHAP Visualization Tool

A production-ready Python tool for generating publication-quality visualizations of molecular structures with atoms highlighted according to their SHAP importance values across multiple model predictions.

## 📋 Overview

This toolkit creates two complementary visualization types:

1. **Individual Model Visualizations** – One image per model showing atoms colored by their SHAP importance values
2. **Consensus Highlights** – A unified view highlighting atoms important in the top 20% for any model, with labels indicating which models prioritize each atom

Perfect for:
- Explaining ensemble model predictions on molecular structures
- Identifying critical atoms for molecular property prediction
- Chemical interpretability and explainability studies
- Publication-quality molecular visualizations

## 🚀 Quick Start

### Installation

```bash
pip install rdkit matplotlib numpy
```

### Basic Usage

```python
from shap_visualizer import visualize_shap_on_molecule
import numpy as np

# Prepare data
y_pred = [0.85, 0.92, 0.78]
shap_vectors = np.array([
    [0.1, 0.8, 0.3, 0.05, 0.9],
    [0.05, 0.7, 0.4, 0.1, 0.85],
    [0.2, 0.6, 0.9, 0.02, 0.7]
])
smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin

# Generate visualizations
output = visualize_shap_on_molecule(y_pred, shap_vectors, smiles)

# Access results
print(output["model_0"])      # → ./shap_plots/model_0_shap.png
print(output["consensus"])    # → ./shap_plots/consensus_highlights.png
```

## 📁 Files

### Core Module
- **shap_visualizer.py** – Main visualization module
  - `MolecularSHAPVisualizer` class with full configuration
  - `visualize_shap_on_molecule()` entry point function
  - Type hints and comprehensive docstrings
  - ~410 lines of production-ready code

### Examples & Documentation
- **example_usage.py** – Complete examples and reference
  - 3 detailed usage examples (Aspirin, Caffeine, complex molecule)
  - ML pipeline integration patterns
  - API reference documentation
  - Troubleshooting guide

## 🎨 Visualization Details

### Individual Model Visualizations (model_i_shap.png)

Each model gets its own visualization:
- **Atoms** colored on gradient: Blue (low importance) → Red (high importance)
- **Color bar** showing SHAP value range for that model
- **Title** displaying: "Model i | Predicted Y = {value}"
- **Bonds** shown in black
- **Atom symbols** centered with black text

### Consensus Highlights (consensus_highlights.png)

A unified view across all models:
- **Red atoms** – In top 20% SHAP importance for at least one model
- **Gray atoms** – Not in top 20% for any model
- **Labels** indicating which models highlighted each atom
  - Format: "M0, M2, M4" (1-3 models)
  - Format: "M0-2, M4" (many models, abbreviated)
- **Title** explaining the visualization

## 🔧 API Reference

### visualize_shap_on_molecule()

```python
visualize_shap_on_molecule(
    y_pred: List[float],
    shap_vectors: np.ndarray,
    smiles_string: str,
    output_dir: str = "./shap_plots",
    highlight_percentile: float = 0.8
) -> Dict[str, str]
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `y_pred` | List[float] | Model predictions (length = n_models) |
| `shap_vectors` | np.ndarray | SHAP values, shape (n_models, n_atoms) |
| `smiles_string` | str | SMILES representation of molecule |
| `output_dir` | str | Output directory (default: "./shap_plots") |
| `highlight_percentile` | float | Percentile threshold (0-1, default: 0.8 for top 20%) |

**Returns:**

Dictionary with file paths:
```python
{
    "model_0": "./shap_plots/model_0_shap.png",
    "model_1": "./shap_plots/model_1_shap.png",
    "model_2": "./shap_plots/model_2_shap.png",
    "consensus": "./shap_plots/consensus_highlights.png"
}
```

**Raises:**
- `ValueError` – Invalid SMILES, mismatched array shapes, or 1D shap_vectors

### MolecularSHAPVisualizer Class

For advanced usage and custom configurations:

```python
from shap_visualizer import MolecularSHAPVisualizer

visualizer = MolecularSHAPVisualizer(
    highlight_percentile=0.8,  # Top 20%
    dpi=300,                   # Resolution
    figure_size=(10, 8),       # Figure dimensions
    font_size=12               # Base font size
)

output = visualizer.visualize(y_pred, shap_vectors, smiles)
```

## 📊 Input Format Examples

### Example 1: Small molecule, 3 models
```python
y_pred = [0.85, 0.92, 0.78]
shap_vectors = np.array([
    [0.1, 0.8, 0.3, 0.05, 0.9],      # 5 atoms, Model 0
    [0.05, 0.7, 0.4, 0.1, 0.85],     # 5 atoms, Model 1
    [0.2, 0.6, 0.9, 0.02, 0.7]       # 5 atoms, Model 2
])
```

### Example 2: Ensemble from trained models
```python
import shap

# Train ensemble
models = [model_1, model_2, model_3]

# Get predictions
y_pred = [m.predict(X)[0] for m in models]

# Compute SHAP values
explainers = [shap.TreeExplainer(m) for m in models]
shap_vectors = np.array([e.shap_values(X)[0] for e in explainers])

# Visualize
visualize_shap_on_molecule(y_pred, shap_vectors, smiles)
```

## 🔍 Understanding the Visualizations

### Color Mapping
- **Blue atoms** = Low SHAP importance (less critical for prediction)
- **Yellow atoms** = Medium importance
- **Red atoms** = High importance (critical for prediction)

Per-model normalization ensures clear distinction even when SHAP ranges differ.

### Consensus Highlights
The consensus view identifies atoms important in **any** model:
- Useful for ensemble agreement analysis
- Shows model disagreement (e.g., M0 and M2 highlight different atoms)
- 80th percentile (top 20%) is default but adjustable

### Percentile Interpretation
| Percentile | Top N% | 10 atoms | 100 atoms |
|-----------|--------|----------|-----------|
| 0.8 | 20% | 2 atoms | 20 atoms |
| 0.9 | 10% | 1 atom | 10 atoms |
| 0.95 | 5% | <1 atom | 5 atoms |

## ⚙️ Configuration & Customization

### Changing Output Resolution
```python
# High-quality publication figures
visualizer = MolecularSHAPVisualizer(dpi=600)

# Draft quality
visualizer = MolecularSHAPVisualizer(dpi=150)
```

### Adjusting Highlight Threshold
```python
# Top 10% only
visualize_shap_on_molecule(..., highlight_percentile=0.9)

# Top 30%
visualize_shap_on_molecule(..., highlight_percentile=0.7)
```

### Custom Figure Size
```python
visualizer = MolecularSHAPVisualizer(figure_size=(14, 10))
```

## 🧪 Testing & Validation

Run the example script to see all features:
```bash
python example_usage.py
```

This displays:
- 3 worked examples with different molecule sizes
- Expected outputs and statistics
- API reference
- Integration patterns
- Troubleshooting guide

## 🔬 Scientific Background

### SHAP Values
SHAP (SHapley Additive exPlanations) values provide:
- Feature-level importance scores
- Consistent, theoretically-grounded explanations
- Per-instance predictions

### Multi-model Consensus
Ensemble predictions often benefit from consensus analysis:
- Identifies robust important features
- Reveals model disagreement
- Better generalization insights

## 📝 Output Examples

### model_0_shap.png
```
┌─────────────────────────────────────┐
│  Model 0 | Predicted Y = 0.850      │
│                                     │
│     Blue atoms: low importance      │
│     Red atoms: high importance      │
│                                     │
│  [Color scale: 0.05 ─── 0.90]       │
└─────────────────────────────────────┘
```

### consensus_highlights.png
```
┌─────────────────────────────────────┐
│ Consensus Highlights - Top 20% SHAP │
│                                     │
│  ● Red: highlighted by model(s)     │
│  ● Gray: not highlighted            │
│  Labels: M0,M2 means models 0 & 2   │
└─────────────────────────────────────┘
```

## 🛠️ Troubleshooting

### Invalid SMILES Error
```python
# Check SMILES validity
from rdkit import Chem
mol = Chem.MolFromSmiles("your_smiles")
assert mol is not None, "Invalid SMILES"
```

### Shape Mismatch
```python
# Ensure compatibility
assert shap_vectors.ndim == 2
assert shap_vectors.shape[0] == len(y_pred)
print(shap_vectors.shape)  # Should be (n_models, n_atoms)
```

### Overlapping Labels
Labels automatically abbreviate for clarity:
- 1-3 models: "M0, M1, M2"
- 4+ models: "M0-3, M5" (compressed format)

## 📚 Dependencies

| Package | Purpose | Version |
|---------|---------|---------|
| rdkit | Molecular structure handling | Latest |
| matplotlib | Visualization | 3.0+ |
| numpy | Numerical operations | 1.19+ |
| pathlib | File operations | Built-in |

## 🎓 Integration with Your Workflow

### Step 1: Train Models
```python
model_1 = your_model.fit(X_train, y_train)
model_2 = your_model.fit(X_train, y_train)
```

### Step 2: Get Predictions
```python
y_pred = [model_1.predict(X_mol), model_2.predict(X_mol)]
```

### Step 3: Compute SHAP
```python
explainers = [shap.Explainer(m) for m in [model_1, model_2]]
shap_vectors = np.array([e(X_mol) for e in explainers])
```

### Step 4: Visualize
```python
visualize_shap_on_molecule(y_pred, shap_vectors, smiles)
```

## 📄 License & Citation

Production-ready code following scientific computing best practices.

## 🤝 Contributing

For improvements, integration patterns, or issues:
- Check docstrings in source code
- Review example_usage.py for patterns
- Ensure type hints are maintained
- Test with molecules of varying sizes

## 📞 Support

Refer to:
1. **Docstrings** in `shap_visualizer.py` for API details
2. **example_usage.py** for worked examples
3. **This README** for overview and integration patterns

---

**Created for molecular machine learning explainability research**

Last updated: December 2025
