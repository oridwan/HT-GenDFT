## Binary Electride Composition Search


This directory contains the script `search_binary_electrides.py`, which is used to identify potential binary electride compositions based on valence electron counting. The script generates a list of candidate compositions suitable for further structure generation with MatterGen.


**Example output formula and calculation:**

   Li13Sb3

This represents a binary electride candidate with 13 lithium atoms and 3 antimony atoms in the formula unit.

**Excess electron calculation:**

- Lithium (Li) valence electrons: 1
- Antimony (Sb) valence electrons: 3 (as defined in the script)

Excess electrons = (valence of Li) × 13 − |valence of Sb| × 3

Excess electrons = 1 × 13 − 3 × 3 = 13 − 9 = 4

So, Li13Sb3 has 4 excess electrons, making it a potential binary electride candidate.

### What the Script Does

- **Purpose:**
  - Searches for binary compounds of the form $A_lC_n$, where $A$ is an electropositive metal (Group I, II, or III) and $C$ is an electronegative non-metal (Group V, VI, or VII).
  - Identifies compositions with a positive excess of valence electrons, which is a key indicator for potential electrides.

- **How it Works:**
  1. **Element Selection:**
     - Uses predefined lists of electropositive metals and electronegative non-metals.
  2. **Composition Generation:**
     - Iterates over possible stoichiometries $l$ and $n$ (number of atoms of $A$ and $C$) within a maximum total atom constraint.
     - Only considers compositions in their simplest integer ratio (i.e., $gcd(l, n) = 1$).

     **Example:**

     Suppose we are considering $A = \text{Li}$ and $C = \text{Sb}$.

     - The script tries combinations of $l$ and $n$ (e.g., $l = 13$, $n = 3$).
     - It checks if the greatest common divisor (gcd) of $l$ and $n$ is 1. For $l = 13$, $n = 3$, $gcd(13, 3) = 1$, so this is already the simplest ratio.
     - If $l = 12$, $n = 6$, then $gcd(12, 6) = 6$, so the simplest ratio would be $l_p = 2$, $n_p = 1$ (by dividing both by 6).
     - Only compositions where $l$ and $n$ are coprime (gcd = 1) are kept, ensuring the formula is in its simplest form (e.g., Li13Sb3, not Li12Sb6).
  3. **Excess Electron Calculation:**
     - For each composition, calculates the excess valence electrons using:
       $$\text{excess} = \text{valence}(A) \times l - |\text{valence}(C)| \times n$$
     - Selects those with excess electrons in a user-defined range (default: 0.1 to 4.0).
  4. **Output:**
     - Saves the valid compositions to a JSON file (`binary_electride_compositions.json`) and a plain text file listing the formulas.
     - Prints summary statistics and example compositions.

### Usage

Run the script directly:

```bash
python search_binary_electrides.py
```

This will generate the output files in the current directory. You can adjust parameters such as the maximum number of atoms, the range of excess electrons, or the maximum number of compositions by modifying the script or adapting the `main()` function.

### Next Steps

1. Review the generated `binary_electride_compositions.json` file.
2. Transfer the file to your cluster or desired location for further processing.
3. Use the generated compositions as input for structure generation workflows (e.g., with MatterGen).

---
For more details, see the comments and documentation within `search_binary_electrides.py`.