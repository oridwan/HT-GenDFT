# SuperConductorFlow Workflow Diagram (Current Repository)

```text
Composition Space
  |
  v
MatterGen Crystal Structure Prediction
  |
  v
Candidate Structures
  |
  v
Prescreen Stability with MatterSim
(Energy Above Convex Hull)
  |
  v
Initial Density Functional Theory Workflow
  |
  +--> Structure Relaxation
        |
        v
      Self-Consistent Static Calculation
        |
        v
      Partial Charge Density
        |
        v
      Electron Localization Function
  |
  v
Electride Analysis and Filtering
  |
  v
Refined Relaxation
(Tight Density Functional Theory, Three-Step Adaptive Relaxation)
  |
  +--> Optional Branch 1: Refined Stability Postcheck
  |
  +--> Optional Branch 2: Electronic Workflow
  |       Self-Consistent Field
  |       -> Partial Charge Density
  |       -> Electron Localization Function
  |       -> Electronic Band Structure
  |       -> Projected Density of States
  |
  +--> Optional Branch 3: Phonon Workflow
  |
  v
Final Electride Set
  |
  v
Final Reports, Plots, and Assets
```

## Notes
- Prescreen threshold is configurable (`--hull-threshold`), default in current script is `0.1 eV/atom`.
- Phonon and detailed electronic workflows run on refined structures with `RELAX_DONE` state (relaxation completed).
- Stability postcheck is optional and can be combined with electronic/phonon outputs.
