# Electride Generation with MatterGen

This directory contains workflows for generating both binary and ternary electride crystal structures using MatterGen's Crystal Structure Prediction (CSP) mode.

## Table of Contents

- [Model Fine-Tuning](#model-fine-tuning)
- [Binary Electride Generation](#binary-electride-generation)
  - [Overview](#binary-overview)
  - [Workflow Steps](#binary-workflow-steps)
  - [File Descriptions](#binary-file-descriptions)
  - [Configuration](#binary-configuration)
- [Ternary Electride Generation](#ternary-electride-generation)
  - [Overview](#ternary-overview)
  - [Prerequisites](#ternary-prerequisites)
  - [Workflow Steps](#ternary-workflow-steps)
  - [Advanced Customization](#advanced-customization)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

---

# Model Fine-Tuning

Before generating structures, you may want to fine-tune the MatterGen base model in CSP (Crystal Structure Prediction) mode for better performance on electride compositions.

## Fine-Tuning MatterGen in CSP Mode

The `mattergen_csp.sh` script trains a MatterGen model specifically for CSP tasks using the Materials Project mp_20 dataset.

### Submit Training Job

```bash
# Make script executable
chmod +x mattergen_csp.sh

# Submit training job (requires GPU cluster)
sbatch mattergen_csp.sh
```

### Training Configuration

**Key parameters:**
- **Config**: CSP mode (pure composition-conditioned generation)
- **Dataset**: mp_20 (Materials Project)
- **Properties**: None (pure CSP without property conditioning)
- **Strategy**: DDP with 2 GPUs
- **Max epochs**: 2000
- **Memory**: 64GB RAM
- **Time limit**: 72 hours

### Monitor Training

```bash
# Check training progress
tail -f mattergen_csp_<jobid>.out

# View error logs if needed
tail -f mattergen_csp_<jobid>.err
```

### Using the Fine-Tuned Model

After training completes, the script will display the checkpoint path. Update this path in your generation scripts:

**For binary electrides:**
```bash
# In binary_electride/generate_binary_csp.sh
CSP_MODEL="$HOME/SOFT/mattergen_test/outputs/singlerun/2025-10-16/18-55-08"
```

**For ternary electrides:**
```bash
# In ternary_electride/generate_ternary_csp.sh
CSP_MODEL="$HOME/SOFT/mattergen_test/outputs/singlerun/2025-10-16/18-55-08"
```

### When to Use Fine-Tuned vs Base Model

- **Fine-tuned model**: Better for binary electrides and specific composition families
- **Base model**: Sufficient for exploratory ternary electride generation
- **Recommendation**: Fine-tune if generating >1000 structures or targeting specific chemistry

---

# Binary Electride Generation

## Binary Overview

Binary electrides follow the formula **A_l C_n** where:
- **A**: Electropositive metals (Group I, II, or III elements)
- **C**: Electronegative non-metals (Group V, VI, or VII elements)

The workflow searches for compositions with excess electrons that could form electride structures, then generates candidate crystal structures using MatterGen.

## Binary Workflow Steps

### 1. Search for Binary Electride Compositions

```bash
cd binary_electride
python search_binary_electrides.py
```

This creates:
- `binary_electride_compositions.json` - Full composition data with excess electrons
- `binary_electride_compositions.txt` - Simple list of formulas

**Key parameters** (in `main()` function):
- `max_atoms=15` - Maximum atoms per unit cell
- `excess_electron_range=(0.1, 2.0)` - Excess electrons for electride formation
- `max_compositions=-1` - Generate all compositions (-1 means no limit)

### 2. Generate Crystal Structures

#### Option A: Single Job (Small Scale)

```bash
sbatch generate_binary_csp.sh
```

**Key parameters** in the script:
- `STRUCTURES_PER_ATOM=2.0` - Number of structures per atom per cell per composition
- `MAX_COMPOSITIONS=-1` - Process all compositions (-1 means all)
- `START_INDEX=0` - Starting composition index

**Note on structure counts**: Each composition is expanded into multiple supercells (up to max 20 atoms). The script automatically distributes structures evenly across all supercell sizes and calculates optimal batch sizes.

#### Option B: Parallel Jobs (Large Scale)

For faster processing of many compositions, use the parallel job submission script:

```bash
cd binary_electride
./submit_parallel_jobs.sh
```

**What it does:**
- Automatically detects total compositions from `binary_electride_compositions.json`
- Splits them into batches (default: 400 compositions per job)
- Submits multiple independent SLURM jobs that run in parallel
- Each job processes a different subset of compositions

**Configuration:**

Edit `submit_parallel_jobs.sh` to adjust batch size:

```bash
COMPOSITIONS_PER_JOB=400  # Increase for fewer larger jobs, decrease for more smaller jobs
```

**Monitoring progress:**

```bash
# Check running jobs
squeue -u $USER

# Count completed compositions
ls -d ../results/binary_csp_electrides/*_structures | wc -l

# View job outputs
tail -f gen_bin_ele_csp_batch*.out
```

**Resume capability**: The generation script automatically skips compositions that already have generated structures, allowing you to resubmit if a job times out or fails.

### 3. Review Results

Generated structures are saved to `../results/binary_csp_electrides/`:

```
../results/binary_csp_electrides/
├── Li1N1_structures/
│   ├── generated_crystals_cif.zip
│   └── generated_crystals.extxyz
├── Li2N1_structures/
│   ├── generated_crystals_cif.zip
│   └── generated_crystals.extxyz
├── generation_statistics.json
└── generation_summary.json
```

## Binary File Descriptions

### Search Scripts
- `search_binary_electrides.py` - Composition search based on valence electron counting

### Generation Scripts
- `generate_structures_batch.py` - Main generation script with supercell expansion
- `generate_binary_csp.sh` - SLURM job script for structure generation
- `submit_parallel_jobs.sh` - Helper to submit multiple parallel jobs
- `summarize_generation.py` - Creates summary statistics after generation

### Output Files
- `binary_electride_compositions.json` - Compositions with metadata
- `binary_electride_compositions.txt` - Simple formula list
- `generation_statistics.json` - Per-composition generation stats
- `generation_summary.json` - Overall generation summary
- `failed_compositions.txt` - List of failed compositions (if any)

## Binary Configuration

### CSP Model Configuration

The binary generation script uses a fine-tuned CSP model for better performance. See the [Model Fine-Tuning](#model-fine-tuning) section for training details.

After fine-tuning with `mattergen_csp.sh`, update the model path in `binary_electride/generate_binary_csp.sh`:

```bash
CSP_MODEL="$HOME/SOFT/mattergen_test/outputs/singlerun/2025-10-16/18-55-08"
```

Replace the path with your trained model's checkpoint location.

### Supercell Expansion

For each composition, the script automatically generates multiple supercells up to 20 atoms. For example, `Li1N1` generates:
- Li1N1 (2 atoms)
- Li2N2 (4 atoms)
- Li3N3 (6 atoms)
- ...
- Li10N10 (20 atoms)

---

# Ternary Electride Generation

## Ternary Overview

This workflow generates crystal structures for potential ternary electride compositions using MatterGen's **Crystal Structure Prediction (CSP) mode** with the base model (no fine-tuning required).

**Strategy:**
1. Search for compositions with excess valence electrons (potential electrides)
2. Generate structures per composition using MatterGen CSP mode

## Ternary Prerequisites

- ASE and Pymatgen installed
- MatterGen base model accessible
- GPU access on cluster

## Ternary Workflow Steps

### Step 1: Search for Ternary Electride Compositions

#### On Local Machine (or Cluster):

```bash
cd ternary_electride
python search_ternary_electrides.py
```

**What it does:**
- Searches through combinations of: `A_l B_m C_n`
  - **A**: Group I/II/III metals + lanthanides (Li, Na, K, Ca, Sr, Ba, La, etc.)
  - **B**: Group III/IV semi-metals (Al, Si, Ge, etc.)
  - **C**: Group V/VI/VII non-metals (N, O, F, P, S, etc.)
- Calculates excess valence electrons: `excess = val(A)×l + val(B)×m - val(C)×n`
- Keeps compositions with 0.1-4.0 excess electrons (electride range)
- Limits to max 20 atoms per composition
- Finds up to 1000 valid compositions

**Output files:**
```
ternary_electride_compositions.json  # Full data with metadata
ternary_electride_compositions.txt   # Just formulas for easy viewing
```

**Example compositions found:**
```
Li3AlN2     | 2.0 e⁻ | 6 atoms
Na3AlN2     | 2.0 e⁻ | 6 atoms  
Ca2SiO3     | 1.0 e⁻ | 6 atoms
K3GaN2      | 2.0 e⁻ | 6 atoms
```

#### Customize the Search:

Edit `search_ternary_electrides.py` to adjust:

```python
compositions = search_ternary_electrides(
    max_atoms=20,                      # Change max atoms per unit cell
    excess_electron_range=(0.1, 4.0),  # Change e⁻ range
    max_compositions=1000              # Change number of compositions
)
```

### Step 2: Transfer to Cluster

```bash
# Transfer composition file and generation script
scp ternary_electride_compositions.json HPC_HOST:~/SOFT/mattergen_test/
scp generate_ternary_electrides.sh HPC_HOST:~/SOFT/mattergen_test/
```

### Step 3: Generate Structures with MatterGen CSP Mode

#### Option A: Single Job (Small to Medium Scale)

```bash
ssh HPC_HOST
cd ~/SOFT/mattergen_test/ternary_electride

# Make script executable
chmod +x generate_ternary_csp.sh

# Submit generation job
sbatch generate_ternary_csp.sh

# Monitor progress
tail -f logs/mattergen_ternary_<jobid>.out
```

#### Option B: Parallel Jobs (Large Scale)

For faster processing of many compositions:

```bash
cd ternary_electride
./submit_parallel_jobs.sh
```

**What it does:**
- Automatically detects total compositions from `ternary_electride_compositions.json`
- Splits them into batches (default: 400 compositions per job)
- Submits multiple independent SLURM jobs that run in parallel
- Each job processes a different subset of compositions

**Configuration:**

Edit `submit_parallel_jobs.sh` to adjust batch size:

```bash
COMPOSITIONS_PER_JOB=400  # Increase for fewer larger jobs, decrease for more smaller jobs
```

**Monitoring progress:**

```bash
# Check running jobs
squeue -u $USER

# Count completed compositions
ls -d ../results/ternary_electrides/*_structures | wc -l

# View job outputs
tail -f gen_ter_ele_csp_batch*.out
```

#### What Happens:

For each composition (e.g., `Li3AlN2`):
1. MatterGen runs in **CSP mode** (conditioned on exact composition)
2. Generates structures adjustable via `STRUCTURES_PER_ATOM`
3. Saves structures as CIF files in `results/ternary_electrides/Li3AlN2_structures/`

**Key Parameters in Script:**

```bash
MODEL="mattergen_base"           # Base model (no fine-tuning)
                                 # Or use fine-tuned model path from mattergen_csp.sh
STRUCTURES_PER_ATOM=2.0          # Target structures per atom
MAX_COMPOSITIONS=-1              # Process all compositions (-1 = all)
```

**Note**: You can use either the base model or a fine-tuned model (see [Model Fine-Tuning](#model-fine-tuning)). The fine-tuned model may provide better results but requires additional training time.

**Note on Structure Counts:**
- Each composition is expanded to multiple supercells (e.g., Li1B1N1, Li2B2N2, ..., Li6B6N6)
- Structures are distributed as evenly as possible across all supercells
- Batch size is automatically calculated per supercell to minimize waste
- Handles remainders intelligently:
  - If not evenly divisible, some supercells get 1 extra structure
  - Difference between any two supercells is at most 1 structure
  
**Examples:**
- 60 structures, 6 supercells: Each gets 10 structures (even split)
- 65 structures, 6 supercells: 5 get 11 structures, 1 gets 10 structures
- 100 structures, 7 supercells: 2 get 15 structures, 5 get 14 structures

The algorithm prefers batch sizes that divide evenly (1, 2, 4, 5, 8, 10, 16, 20, 32, 64)

#### CSP Mode vs Property Conditioning:

| Mode | Input | Output | Use Case |
|------|-------|--------|----------|
| **CSP** | Chemical formula (e.g., Li3AlN2) | Structures with exact composition | Specific composition search |
| **Property** | Property value (e.g., is_electride=1) | Structures with that property | Exploration/discovery |

**We use CSP mode here** because we have specific target compositions from our search.

### Step 4: Review Results

```bash
# Check summary
cat results/ternary_electrides/generation_summary.json

# List generated compositions
ls results/ternary_electrides/*/

# Count structures
find results/ternary_electrides -name "*.cif" | wc -l
```

**Expected output structure:**
```
results/ternary_electrides/
├── Li3AlN2_structures/
│   ├── structure_000.cif
│   ├── structure_001.cif
│   ├── ...
│   └── structure_019.cif (20 total)
├── Na3AlN2_structures/
│   ├── structure_000.cif
│   ├── ...
├── generation_summary.json
└── failed_compositions.txt (if any failed)
```

## Advanced Customization

### Customize Element Selection

Edit `search_ternary_electrides.py` to try different element combinations:

```python
# Example 1: Focus on alkali metals only
GROUP_A = ['Li', 'Na', 'K', 'Rb', 'Cs']
GROUP_B = ['Al', 'Si']
GROUP_C = ['N', 'O']

# Example 2: Include rare earths more prominently
GROUP_A = ['La', 'Ce', 'Pr', 'Nd', 'Y']
GROUP_B = ['Si', 'Ge', 'Sn']
GROUP_C = ['N', 'O']

# Example 3: Search binary compositions (set GROUP_C to single element)
GROUP_A = ['Ca', 'Sr', 'Ba']
GROUP_B = ['Al']
GROUP_C = ['N']  # Only nitrides
```

---

# Troubleshooting

## Binary Electride Issues

### GPU Memory Issues
If you encounter OOM errors, reduce `STRUCTURES_PER_ATOM` in the generation script.

### Job Timeout
Use `submit_parallel_jobs.sh` to split the work across multiple jobs, each processing fewer compositions.

### Resume Generation
Simply resubmit the job - it will automatically skip compositions that already have generated structures.

## Ternary Electride Issues

### Problem: Too many compositions (>10,000)

**Solution:** Reduce search space:
```python
compositions = search_ternary_electrides(
    max_atoms=15,  # Reduce from 20
    max_compositions=500  # Cap at 500
)
```

### Problem: Generation is too slow

**Solution:** Use parallel job submission (see Option B in generation steps):

For binary electrides:
```bash
cd binary_electride
./submit_parallel_jobs.sh
```

For ternary electrides:
```bash
cd ternary_electride
./submit_parallel_jobs.sh
```

This automatically splits compositions into batches and runs multiple jobs in parallel. Adjust `COMPOSITIONS_PER_JOB` in the script to control batch size.

### Problem: Out of GPU memory

**Solution:** The batch size is auto-calculated. If you get OOM errors, reduce structures per composition:
```bash
STRUCTURES_PER_ATOM=1.0  # Reduce from 2.0
# This will automatically use smaller batch sizes per supercell
```

### Problem: Some compositions fail to generate

**Solution:** Check `failed_compositions.txt` and retry:
```bash
# Failed compositions are saved automatically
cat results/[binary_csp_electrides|ternary_electrides]/failed_compositions.txt

# Manually retry specific compositions
mattergen-generate results/retry_Li3AlN2 \
    --model_path=mattergen_base \
    --sampling_config_name=csp \
    --target_compositions='[{"Li": 3, "Al": 1, "N": 2}]' \
    --batch_size=20 \
    --num_batches=1
```

---

# Workflow Summary

```
┌─────────────────────────────────────────────────────────────┐
│ 1. SEARCH COMPOSITIONS (Local/Cluster)                      │
│    Binary:  python search_binary_electrides.py              │
│    Ternary: python search_ternary_electrides.py             │
│    → composition files with excess electron candidates      │
└───────────────────┬─────────────────────────────────────────┘
                    │
┌───────────────────▼─────────────────────────────────────────┐
│ 2. GENERATE STRUCTURES (Cluster - GPU)                      │
│    Binary:  sbatch generate_binary_csp.sh                   │
│    Ternary: sbatch generate_ternary_electrides.sh           │
│    → Multiple structures per composition                    │
│    → Results saved to results/ directory                    │
└─────────────────────────────────────────────────────────────┘
```

## Expected Results

**From binary compositions:**
- Generate structures based on STRUCTURES_PER_ATOM setting
- Review `generation_summary.json` for statistics
- Structures saved as CIF and extxyz formats

**From ~1000 ternary compositions:**
- Generate ~20,000 structures (20 per composition)
- Structures organized by composition in results directory

**Estimated Timeline:**
- Composition search: ~1 minute
- Structure generation: 24-72 hours (depends on number of compositions and GPU availability)

---

# Citation

If you use this workflow, cite:
- **MatterGen:** Zeni et al., Nature (2025) - https://doi.org/10.1038/s41586-025-08628-5

---

# Summary

This workflow enables **high-throughput computational discovery** of binary and ternary electride materials by:
1. Systematic composition space search based on valence electron counting
2. AI-driven structure generation using MatterGen CSP mode
