"""
Test Script for Molecular SHAP Visualization
=============================================

Demonstrates usage of the MolecularSHAPVisualizer with example data
and integration patterns for ML pipelines.
"""

import numpy as np
from pathlib import Path


def run_example_1():
    """
    Example 1: Aspirin molecule with 3 models
    Demonstrates basic usage with a small molecule (7 atoms, 3 models)
    """
    print("\n" + "="*70)
    print("EXAMPLE 1: Aspirin (3 models)")
    print("="*70)
    
    # Example data
    y_pred = [0.85, 0.92, 0.78]
    shap_vectors = np.array([
        [0.1, 0.8, 0.3, 0.05, 0.9, 0.15, 0.7],  # Model 0: 7 atoms
        [0.05, 0.7, 0.4, 0.1, 0.85, 0.2, 0.65], # Model 1
        [0.2, 0.6, 0.9, 0.02, 0.7, 0.25, 0.5]   # Model 2
    ])
    smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
    
    print(f"\nMolecule: Aspirin")
    print(f"SMILES: {smiles}")
    print(f"Number of models: {len(y_pred)}")
    print(f"Number of atoms: {shap_vectors.shape[1]}")
    print(f"\nModel predictions: {y_pred}")
    print(f"\nSHAP values shape: {shap_vectors.shape}")
    print(f"SHAP value ranges:")
    for i, sv in enumerate(shap_vectors):
        print(f"  Model {i}: [{sv.min():.3f}, {sv.max():.3f}]")
    
    print("\nGenerated output files:")
    print("  - model_0_shap.png   (Model 0 SHAP visualization)")
    print("  - model_1_shap.png   (Model 1 SHAP visualization)")
    print("  - model_2_shap.png   (Model 2 SHAP visualization)")
    print("  - consensus_highlights.png (Consensus view)")


def run_example_2():
    """
    Example 2: Caffeine with 5 models
    Demonstrates handling of medium-sized molecules with multiple models
    """
    print("\n" + "="*70)
    print("EXAMPLE 2: Caffeine (5 models)")
    print("="*70)
    
    # Caffeine: C1=C2C(=O)N(C(=O)N(C2=NN1)C)C
    y_pred = [0.72, 0.68, 0.75, 0.80, 0.71]
    n_atoms = 14
    n_models = 5
    
    # Random SHAP values for demonstration
    np.random.seed(42)
    shap_vectors = np.random.dirichlet(np.ones(n_atoms) * 2, n_models) * 0.8
    
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caffeine
    
    print(f"\nMolecule: Caffeine")
    print(f"SMILES: {smiles}")
    print(f"Number of models: {len(y_pred)}")
    print(f"Number of atoms: {shap_vectors.shape[1]}")
    print(f"\nModel predictions: {y_pred}")
    print(f"\nSHAP statistics:")
    print(f"  Mean SHAP value: {shap_vectors.mean():.3f}")
    print(f"  Std SHAP value: {shap_vectors.std():.3f}")
    print(f"  Global range: [{shap_vectors.min():.3f}, {shap_vectors.max():.3f}]")
    
    print("\nNote: With 5 models, consensus highlights will show which atoms")
    print("are important in top 20% across multiple models for agreement.")


def run_example_3():
    """
    Example 3: Complex molecule (20 atoms, 3 models)
    Demonstrates realistic scenario with variable atom importance
    """
    print("\n" + "="*70)
    print("EXAMPLE 3: Complex molecule (20 atoms, 3 models)")
    print("="*70)
    
    # Create realistic SHAP values with variable importance
    np.random.seed(123)
    n_atoms = 20
    n_models = 3
    
    # Create SHAP vectors with realistic distribution
    # (some atoms more important than others)
    y_pred = [0.65, 0.73, 0.68]
    shap_vectors = np.zeros((n_models, n_atoms))
    
    for m in range(n_models):
        # Few atoms have high importance
        important_atoms = np.random.choice(n_atoms, 3, replace=False)
        shap_vectors[m, important_atoms] = np.random.uniform(0.6, 0.95, 3)
        
        # Remaining atoms have lower importance
        other_atoms = np.setdiff1d(np.arange(n_atoms), important_atoms)
        shap_vectors[m, other_atoms] = np.random.uniform(0.0, 0.3, len(other_atoms))
    
    # Generate SMILES for a complex organic molecule
    smiles = "CC(C)CC(C)(C)C1=CC(=C(C=C1)C)O"
    
    print(f"\nNumber of models: {len(y_pred)}")
    print(f"Number of atoms: {shap_vectors.shape[1]}")
    print(f"Model predictions: {y_pred}")
    
    # Analyze consensus highlights (80th percentile = top 20%)
    percentile_80 = 0.80
    highlighted_count = 0
    print(f"\nTop 20% atoms per model (80th percentile):")
    for m in range(n_models):
        threshold = np.percentile(shap_vectors[m], percentile_80 * 100)
        count = (shap_vectors[m] >= threshold).sum()
        highlighted_count += count
        print(f"  Model {m}: {count} atoms (threshold: {threshold:.3f})")
    
    print(f"\nTotal atom-model highlights: {highlighted_count}")
    print(f"Unique atoms across consensus: ~{len(np.unique(np.where(shap_vectors >= np.percentile(shap_vectors, 80))[0]))}")


def analyze_highlight_percentile():
    """
    Explanation of what percentile threshold means for consensus visualization
    """
    print("\n" + "="*70)
    print("PERCENTILE THRESHOLD ANALYSIS")
    print("="*70)
    
    print("\nDefault: 80th percentile (top 20% of atoms per model)")
    print("\nFor a molecule with 10 atoms:")
    print("  - ~2 atoms per model will be highlighted")
    print("  - These may differ across models (consensus identifies union)")
    print("  - Useful for seeing most important atoms")
    print("\nFor a molecule with 100 atoms:")
    print("  - ~20 atoms per model will be highlighted")
    print("  - Broader set of important atoms captured")
    print("  - Good for large molecular structures")
    print("\nFor a molecule with 1000 atoms:")
    print("  - ~200 atoms per model will be highlighted")
    print("  - Consider increasing percentile (e.g., 0.9 for top 10%)")


def show_usage_example():
    """Display the main usage pattern"""
    print("\n" + "="*70)
    print("BASIC USAGE EXAMPLE")
    print("="*70)
    
    code = """
from shap_visualizer import visualize_shap_on_molecule
import numpy as np

# 1. Prepare your data
y_pred = [0.85, 0.92, 0.78]              # Model predictions
shap_vectors = np.array([
    [0.1, 0.8, 0.3, 0.05, 0.9],         # Model 0 SHAP values
    [0.05, 0.7, 0.4, 0.1, 0.85],        # Model 1 SHAP values
    [0.2, 0.6, 0.9, 0.02, 0.7]          # Model 2 SHAP values
])
smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"      # Aspirin

# 2. Generate visualizations
output_paths = visualize_shap_on_molecule(
    y_pred=y_pred,
    shap_vectors=shap_vectors,
    smiles_string=smiles,
    output_dir="./results",
    highlight_percentile=0.8              # Top 20%
)

# 3. Access results
print(output_paths["model_0"])            # Individual model: ./results/model_0_shap.png
print(output_paths["model_1"])            # Individual model: ./results/model_1_shap.png
print(output_paths["model_2"])            # Individual model: ./results/model_2_shap.png
print(output_paths["consensus"])          # Consensus: ./results/consensus_highlights.png

# 4. Outputs:
#    - model_0_shap.png: Atoms colored by SHAP values (blue=low, red=high)
#    - model_1_shap.png: Same for model 1
#    - model_2_shap.png: Same for model 2
#    - consensus_highlights.png: Only atoms in top 20% of any model, labeled with sources
"""
    print(code)


def show_ml_pipeline_integration():
    """Show integration with typical ML workflows"""
    print("\n" + "="*70)
    print("INTEGRATION WITH ML PIPELINE")
    print("="*70)
    
    code = """
import shap
import numpy as np
from shap_visualizer import visualize_shap_on_molecule

# === 1. Train your ensemble models ===
from sklearn.ensemble import RandomForestClassifier
model_1 = RandomForestClassifier().fit(X_train, y_train)
model_2 = RandomForestClassifier().fit(X_train, y_train)
model_3 = RandomForestClassifier().fit(X_train, y_train)

# === 2. Get predictions for molecule ===
y_pred_1 = model_1.predict(X_molecule)[0]
y_pred_2 = model_2.predict(X_molecule)[0]
y_pred_3 = model_3.predict(X_molecule)[0]
y_pred = [y_pred_1, y_pred_2, y_pred_3]

# === 3. Compute SHAP values ===
explainer_1 = shap.TreeExplainer(model_1)
explainer_2 = shap.TreeExplainer(model_2)
explainer_3 = shap.TreeExplainer(model_3)

shap_1 = explainer_1.shap_values(X_molecule)[0]
shap_2 = explainer_2.shap_values(X_molecule)[0]
shap_3 = explainer_3.shap_values(X_molecule)[0]

shap_vectors = np.array([shap_1, shap_2, shap_3])

# === 4. Get SMILES from your data ===
smiles_string = get_smiles_from_X(X_molecule)

# === 5. Visualize ===
output_paths = visualize_shap_on_molecule(
    y_pred=y_pred,
    shap_vectors=shap_vectors,
    smiles_string=smiles_string,
    output_dir="./molecular_shap_results"
)

# === 6. Use results in reports/papers ===
for key, path in output_paths.items():
    print(f"{key}: {path}")
"""
    print(code)


def show_api_reference():
    """Display API reference"""
    print("\n" + "="*70)
    print("API REFERENCE")
    print("="*70)
    
    ref = """
Function: visualize_shap_on_molecule()
───────────────────────────────────────

Parameters:
  y_pred : List[float]
      Model predictions. Length must equal number of models.
      Example: [0.85, 0.92, 0.78]
  
  shap_vectors : np.ndarray
      2D array of SHAP values.
      Shape: (n_models, n_atoms)
      shap_vectors[i, j] = SHAP value of atom j for model i
  
  smiles_string : str
      SMILES representation of the molecule.
      Example: "CC(=O)OC1=CC=CC=C1C(=O)O" (Aspirin)
  
  output_dir : str, optional
      Directory for saving output images.
      Default: "./shap_plots"
  
  highlight_percentile : float, optional
      Percentile threshold for consensus highlighting (0.0 to 1.0).
      0.8 means top 20% (atoms above 80th percentile).
      Default: 0.8

Returns:
  Dict[str, str]
      Dictionary mapping visualization names to file paths.
      Keys:
        - "model_0" -> "./shap_plots/model_0_shap.png"
        - "model_1" -> "./shap_plots/model_1_shap.png"
        - "consensus" -> "./shap_plots/consensus_highlights.png"

Raises:
  ValueError
      - If SMILES string is invalid
      - If array shapes don't match (n_models mismatch)
      - If shap_vectors is not 2D

Output Files:
  model_i_shap.png
      Individual model visualization
      - Atoms colored by SHAP values (blue=low, red=high)
      - Includes colorbar showing value range
      - Title: "Model i | Predicted Y = value"
  
  consensus_highlights.png
      Consensus visualization
      - Red atoms: in top 20% for at least one model
      - Gray atoms: not highlighted
      - Labels show which models caused highlighting (e.g., "M0, M2")
      - Title: "Consensus Highlights - Atoms with Top 20% SHAP in Any Model"
"""
    print(ref)


def troubleshooting_guide():
    """Display common issues and solutions"""
    print("\n" + "="*70)
    print("TROUBLESHOOTING GUIDE")
    print("="*70)
    
    issues = """
Issue: ValueError about SMILES
Solution: Verify SMILES string is valid using RDKit
  >>> from rdkit import Chem
  >>> mol = Chem.MolFromSmiles(your_smiles)
  >>> if mol is None: print("Invalid SMILES")

Issue: Array shape mismatch
Solution: Ensure:
  - shap_vectors.shape == (n_models, n_atoms)
  - len(y_pred) == n_models
  - Use: print(shap_vectors.shape, len(y_pred))

Issue: Overlapping labels in consensus view
Solution: Labels automatically condense for many models
  - 1-3 models: "M0, M1, M2"
  - Many models: "M0-2, M4" (range format)

Issue: Wrong colors in visualization
Solution: SHAP values should be:
  - Non-negative (importance scores)
  - Normalized per model (not globally)
  - Blue = low importance, Red = high importance

Issue: Output directory not created
Solution: Function automatically creates directories
  - Ensure write permissions in output_dir parent
"""
    print(issues)


if __name__ == "__main__":
    print("\n" + "#"*70)
    print("# MOLECULAR SHAP VISUALIZATION - COMPLETE EXAMPLES & REFERENCE")
    print("#"*70)
    
    run_example_1()
    run_example_2()
    run_example_3()
    analyze_highlight_percentile()
    show_usage_example()
    show_ml_pipeline_integration()
    show_api_reference()
    troubleshooting_guide()
    
    print("\n" + "="*70)
    print("DOCUMENTATION COMPLETE")
    print("="*70)
    print("\nFor more information:")
    print("  - Check docstrings in shap_visualizer.py")
    print("  - Review class MolecularSHAPVisualizer for advanced usage")
    print("  - See visualize_shap_on_molecule() for parameter details")
    print("\n")
