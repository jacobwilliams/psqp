# Adaptive Penalty Parameter Updates

## Overview

Implemented automatic adaptive penalty parameter (`rpf`) updates in PSQP, eliminating the need for manual tuning and preventing convergence failures due to poorly chosen penalty coefficients.

**Priority**: #1 improvement from ALGORITHM_ANALYSIS.md  
**Effort**: Low-Medium  
**Impact**: High - Improves robustness and eliminates user tuning

---

## What Was Implemented

### 1. Dual Adaptive Strategy

The penalty parameter `rpf` (used in the augmented Lagrangian merit function) is now automatically adjusted during optimization using two complementary mechanisms:

#### Strategy A: Constraint Stagnation Detection
- Monitors constraint violation progress between major iterations
- If `cmax > threshold * cmax_previous` for multiple iterations:
  - Increment stagnation counter
  - After `rpf_stagnation_limit` stagnations, increase `rpf` by `rpf_increase_factor`
  - This handles cases where constraints aren't improving

#### Strategy B: Large Lagrange Multiplier Detection
- Monitors magnitude of active Lagrange multipliers
- If `max|λ| > 100 * rpf`:
  - Immediately increase `rpf` to `max|λ| / 10`
  - This prevents penalty parameter from being orders of magnitude too small

### 2. Configuration Parameters

Added public tuning parameters to `psqp_class`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `rpf_adaptive` | `.true.` | Enable/disable adaptive updates |
| `rpf_min` | `1.0e-6` | Minimum penalty parameter value |
| `rpf_max` | `1.0e6` | Maximum penalty parameter value |
| `rpf_increase_factor` | `10.0` | Multiplier for stagnation-triggered increases |
| `rpf_stagnation_threshold` | `0.9` | Threshold for detecting constraint stagnation |
| `rpf_stagnation_limit` | `3` | Iterations of stagnation before increase |

### 3. Internal State Tracking

Added private member variables for algorithm state:
- `cmax_previous`: Constraint violation from previous major iteration
- `rpf_stagnation_count`: Counter for consecutive stagnations

---

## Changes to Code

### Modified Files

1. **src/psqp_module.f90**
   - Added 6 public configuration parameters to `psqp_class`
   - Added 2 private state tracking variables
   - Implemented adaptive update logic in main optimization loop (after line 627)
   - Print messages when `rpf` is increased (when `iprnt > 1`)

### Test Files

2. **test/test_adaptive_penalty.f90** (new)
   - Demonstrates adaptive mechanism with deliberately poor initial `rpf = 1.0e-6`
   - Problem: constrained quadratic minimization
   - Shows automatic increase to `rpf ≈ 0.18` (factor of 180,000x)
   - Verifies successful convergence despite poor initial parameter

---

## How to Use

### Default Behavior (Recommended)

Adaptive mode is enabled by default. Just set an initial guess for `rpf`:

```fortran
rpar(5) = 1.0e-4_wp  ! Initial penalty parameter (will adapt automatically)

call solver%psqpn(nf, nb, nc, x, ix, xl, xu, cf, ic, cl, cu, &
                  ipar, rpar, f, gmax, cmax, iprnt, iterm, &
                  obj_func, dobj_func, con_func, dcon_func)
```

### Custom Configuration

Adjust adaptation parameters if needed:

```fortran
type(psqp_class) :: solver

! Enable aggressive adaptation
solver%rpf_increase_factor = 100.0_wp  ! Increase faster
solver%rpf_stagnation_limit = 1        ! Trigger after 1 stagnation

! Or disable adaptation (use fixed rpf)
solver%rpf_adaptive = .false.
```

---

## Verification

### Existing Tests

All existing tests pass with adaptive mode enabled:
- `test_sparse_optimization`: Converges without triggering adaptation (well-conditioned)
- `test_sparse_nonlinear`: Converges without triggering adaptation (smooth progress)

### New Test

`test_adaptive_penalty`:
- Starts with `rpf = 1.0e-6` (intentionally poor)
- Adaptive mechanism increases to `rpf = 0.18` on iteration 1
- Converges successfully in 3 iterations
- Constraint violation: `cmax = 0.0`
- **Result**: ✓ PASSED

---

## Algorithm Details

### Update Location

Adaptive updates occur in the main optimization loop:

```fortran
! Line 627: After computing new penalty parameters for constraints
call compute_new_penalty_parameters(nf, n, nc, ica, cz, cp)

! NEW: Adaptive rpf update logic
if (me%rpf_adaptive .and. me%nit > 0) then
   ! Check stagnation
   if (cmax > me%rpf_stagnation_threshold * me%cmax_previous) then
      me%rpf_stagnation_count = me%rpf_stagnation_count + 1
      if (me%rpf_stagnation_count >= me%rpf_stagnation_limit) then
         rpf = min(rpf * me%rpf_increase_factor, me%rpf_max)
         me%rpf_stagnation_count = 0
      end if
   else
      me%rpf_stagnation_count = 0
   end if
   
   ! Check large multipliers
   if (nf - n > 0) then
      rp = maxval(abs(cz(1:nf-n)))
      if (rp > 100.0_wp * rpf) then
         rpf = min(rp / 10.0_wp, me%rpf_max)
      end if
   end if
end if
me%cmax_previous = cmax
```

### Augmented Lagrangian Merit Function

The penalty parameter `rpf` scales constraint violations:

```
φ(x) = f(x) + rpf * Σ|constraint_violations| - Σ(λ_i * active_constraint_residuals)
```

Adaptive updates ensure `rpf` is:
- Large enough to enforce constraints
- Not so large that it causes numerical issues
- Adjusted during optimization rather than requiring user expertise

---

## Benefits

✅ **No Manual Tuning**: Users no longer need to guess initial `rpf` values  
✅ **Robust Convergence**: Prevents failures from poor penalty parameters  
✅ **Automatic Scaling**: Adapts to problem conditioning during solve  
✅ **Transparent**: Print messages show when/why `rpf` increases (with `iprnt > 1`)  
✅ **Backward Compatible**: Existing code works unchanged, gains automatic adaptation

---

## References

- Nocedal & Wright, *Numerical Optimization*, 2nd ed., Chapter 18 (Penalty and Augmented Lagrangian Methods)
- SNOPT Technical Reports on automatic penalty parameter selection
- Fletcher & Powell (1963) - Original penalty/barrier method theory

---

## Future Enhancements

Possible extensions (not currently implemented):

1. **Decrease rpf**: Currently only increases; could decrease when constraints easily satisfied
2. **Per-constraint penalties**: Adapt `cp(i)` array more aggressively (currently only rpf adapts)
3. **Heuristic refinement**: Tune thresholds (100x, 0.9, etc.) based on empirical testing
4. **Filter method**: Alternative globalization strategy (see improvement #5 in analysis)

---

## Status

✅ **Implemented**: January 2026  
✅ **Tested**: All existing tests pass, new test demonstrates functionality  
✅ **Documented**: This file + inline comments in code  
✅ **Default**: Enabled by default for all users
