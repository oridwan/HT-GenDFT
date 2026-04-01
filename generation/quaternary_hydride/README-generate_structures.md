## Structure Generation with MatterGen CSP

This README describes how to use `generate_structures_batch.py` to generate crystal structures for multiple compositions using MatterGen in CSP mode.

### Purpose

- Automates the generation of crystal structures for a batch of compositions.
- Supports both fine-tuned MatterGen checkpoints and pretrained models.
- Expands each composition to all valid supercells up to a maximum atom count, allowing exploration of different cell sizes.

### Workflow

1. **Input:**
   - Takes a JSON file containing a list of compositions (e.g., from binary electride search).
2. **Supercell Expansion:**
   - For each composition, generates all possible supercells (multiples of the base formula) that do not exceed the maximum atom count.
   - Example: For Li1B3P2 (6 atoms), supercells are Li1B3P2 (6 atoms), Li2B6P4 (12 atoms), Li3B9P6 (18 atoms).
3. **Structure Generation:**
   - Calls `mattergen-generate` for each supercell, distributing the requested number of structures proportionally to atom count.
   - Handles batching to avoid GPU out-of-memory errors.
4. **Output:**
   - Combines generated CIF and extxyz files for all supercells into a single output directory per composition.
   - Saves summary statistics and error logs.

### Usage

Run the script from the command line:

```bash
python generate_structures_batch.py \
  --compositions binary_electride_compositions.json \
  --output-dir ./output_structures \
  --model mattergen_base \
  --n-structures 20
```

#### Key Options

- `--compositions` / `-c`: Path to the JSON file with compositions.
- `--output-dir` / `-o`: Base directory for output structures.
- `--model`: MatterGen model checkpoint path or pretrained name (e.g., `mattergen_base`).
- `--n-structures` / `-n`: Number of structures to generate per composition (default: 20).
- `--structures-per-atom`: Proportional mode, generates a set number of structures per atom across all supercells.
- `--max-atoms`: Maximum atoms per cell for supercell expansion (default: 20).
- `--max-batch-size`: Maximum structures per GPU batch (default: 100).
- `--timeout`: Timeout per batch in seconds (default: 1800).
- `--max-compositions` / `-m`: Maximum number of compositions to process.
- `--start-index`: Start index in the compositions list.
- `--skip-existing`: Skip compositions that already have generated structures (resume mode).

### Example: Supercell Expansion

For a composition like Li1B3P2 (6 atoms):

- 1x: Li1B3P2 (6 atoms)
- 2x: Li2B6P4 (12 atoms)
- 3x: Li3B9P6 (18 atoms)

Supercells larger than 20 atoms are excluded if `max_atoms=20`.

### Output

- For each composition, the script creates a directory containing:
  - Combined CIF files (`generated_crystals_cif.zip`)
  - Combined extxyz file (`generated_crystals.extxyz`)
  - Error logs if generation fails
  - Generation statistics (`generation_statistics.json`)

### Error Handling

- If structure generation fails for any composition, details are logged and temporary directories are preserved for debugging.

### Notes

- The script is designed to be robust for large batch jobs and can resume from previous runs.
- Adjust batch size and timeout as needed for your hardware.

---
For more details, see the comments and documentation within `generate_structures_batch.py`.