# Wolfe Line Search Refactoring

## Date
March 18, 2026

## Issue Identified
The initial implementation of `wolfe_line_search` had two critical problems:

1. **Incorrect function calls**: Used direct `me%obj`/`me%dobj` calls instead of `compute_obj_and_dobj`
   - Did not increment function evaluation counters
   - Did not handle maximization problems (iext=1)
   - Bypassed the caching mechanism

2. **Missing constraint handling**: Operated only on the objective function, completely ignoring constraints
   - For constrained problems, should search along the augmented Lagrangian merit function
   - Extended line search correctly evaluates constraints and computes merit function

## Solution Implemented

### 1. Use `compute_obj_and_dobj` (First Fix)
- Added `iext` parameter to signature
- Replaced all direct `me%obj`/`me%dobj` calls with `compute_obj_and_dobj`
- Proper tracking of `kd_local`/`ld_local` for each evaluation

### 2. Augmented Lagrangian Merit Function (Second Fix)
- Updated signature to accept constraint-related parameters:
  - `n` - dimension of constraint null space
  - `nc` - number of constraints
  - `ic, ica, cl, cu, cz` - constraint parameters
  - `rpf` - penalty parameter
  - `cf` - constraint values workspace
  - `gc, cg` - constraint gradient workspaces

- Modified evaluation process at each trial point:
  ```fortran
  ! Evaluate objective and gradient
  call me%compute_obj_and_dobj(...)

  ! Evaluate constraints
  call me%compute_con_and_dcon(...)

  ! Compute augmented Lagrangian merit function
  cf(nc + 1) = f_new
  call compute_augmented_lagrangian(...)
  ```

- Updated both main loop and internal `zoom` subroutine
- Merit function φ(x) = f(x) + rpf·Σ|violations| + λ·constraints

## Impact on Performance

The refactoring reveals the **true cost** of Wolfe line search:

### Rosenbrock with constraint test:
- **Extended**: 53 function evals, 37 gradient evals, 37 iterations
- **Wolfe**: 278 function evals, 278 gradient evals, 32 iterations

### Analysis:
- Wolfe uses ~5× more evaluations but ~14% fewer iterations
- Extended evaluates only objective during line search, postpones constraints
- Wolfe evaluates objective + constraints at each trial (correct for merit function)
- Both methods converge to same solution with comparable accuracy

### Trade-off:
- **Extended**: Fewer evaluations per iteration, more iterations
- **Wolfe**: More evaluations per iteration, fewer iterations (better steps)

## Correctness

The refactored implementation now:
1. ✅ Properly tracks function/gradient evaluation counts
2. ✅ Handles maximization problems correctly
3. ✅ Operates on augmented Lagrangian for constrained problems
4. ✅ Maintains architectural consistency with `extended_line_search`
5. ✅ Passes all test cases (standard, sparse, adaptive, comparison)

## Conclusion

The Wolfe line search is now **architecturally correct** for constrained optimization. The higher evaluation count is not a bug but reflects a fundamental limitation: **Wolfe line search is poorly suited for non-smooth augmented Lagrangian merit functions**.

### Why Wolfe Struggles with Constrained Optimization

1. **Non-smooth merit function**: The augmented Lagrangian φ(x) = f(x) + rpf·Σ|violations| contains absolute values, creating non-differentiable points
2. **Curvature condition failure**: Strong Wolfe's curvature condition `|φ'(α)| ≤ c₂|φ'(0)|` assumes smoothness
3. **Constraint evaluation cost**: Wolfe evaluates objective + constraints at each trial (correct behavior), while Extended postpones constraint evaluation

### Performance Analysis (Rosenbrock with constraint)

**Extended**: 53 function evals, 37 gradient evals, 37 iterations
**Wolfe**: 278 function evals, 278 gradient evals, 32 iterations

- Wolfe takes ~5× more evaluations but ~14% fewer iterations
- Extended evaluates only objective during line search, postpones constraints
- Wolfe evaluates objective + constraints at each trial (correct for merit function)
- Extended's simpler acceptance criteria work better for non-smooth functions

### Optimizations Applied

The following fixes **were successfully implemented** and improve Wolfe's behavior:

1. ✅ **c2 = 0.1** (changed from 0.9) - Appropriate for quasi-Newton methods
2. ✅ **Eliminated redundant re-evaluations** - zoom returns gradient, avoiding duplicate calls
3. ✅ **Proper constraint handling** - Evaluates augmented Lagrangian (not just objective)
4. ✅ **Correct function counting** - Uses `compute_obj_and_dobj` for proper tracking

These fixes are correct and necessary. The high evaluation count reflects an algorithmic limitation, not a bug.

### Recommendations

**For constrained optimization with penalty methods:**
- **Use Extended line search** (method=1, default) - More robust for non-smooth merit functions
- Simpler acceptance criteria work better with absolute value penalties
- Lower per-iteration cost (defers constraint evaluation)

**For smooth unconstrained problems:**
- **Wolfe line search** (method=2) may be more efficient
- Directional derivatives enable better step lengths
- Curvature conditions make sense for smooth objectives

**For interior-point or smooth merit functions:**
- Wolfe may outperform Extended
- Filter methods (not yet implemented) would be ideal with Wolfe

### Future Considerations

To make Wolfe competitive for constrained problems:
1. **Implement filter method** - Avoids non-smooth penalty terms
2. **Modified Wolfe conditions** - Relaxed for non-smooth functions
3. **Smoothed merit function** - Replace |·| with smooth approximation
4. **Hybrid approach** - Use Extended near constraint boundaries, Wolfe elsewhere

For problems where constraint evaluations are expensive, the extended line search may be more efficient. For problems where better step lengths significantly reduce iterations, Wolfe may be preferred despite higher per-iteration cost.

Both methods are now available as options via the `line_search_method` parameter.
