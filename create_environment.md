# Setup Mattersim & MatterGen Environment

## Quick Setup

### 1. Create Environment
```bash
cd /projects/mmi/Ridwan/Github/superconductor
conda env create -f environment.yml
```

### 2. Activate Environment
```bash
conda activate mattersim
```

### 3. Install MatterGen
```bash
pip install mattergen==1.0.3
```

### 4. Verify Installation
```bash
python -c "import mattergen; import mattersim; print('✓ Ready to use mattergen and mattersim')"
```

---

## Verify PyTorch
```bash
python -c "import torch; print(f'PyTorch: {torch.__version__}'); print(f'CUDA: {torch.cuda.is_available()}')"
```

## Environment Details
- **Name:** mattersim
- **Python:** 3.10.18
- **PyTorch:** 2.2.1
- **Includes:** Mattersim 1.1.2, E3NN, PyMatGen, Lightning, Hydra, and all dependencies

## Deactivate
```bash
conda deactivate
```
