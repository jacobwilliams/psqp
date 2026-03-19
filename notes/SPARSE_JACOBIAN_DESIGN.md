# Sparse Jacobian Support - Design Document

## Overview
Add sparse Jacobian support to PSQP while maintaining backward compatibility with dense mode.

## Current Architecture

### Constraint Jacobian Storage (Dense Mode)
- Stored in `cg(nf*nc)` - full dense matrix
- Column-major format: constraint `kc` gradient at indices `(kc-1)*nf+1` to `kc*nf`
- User callback `dcon(me, nf, kc, x, gc)` computes one constraint gradient at a time
- Memory: O(nf × nc)

### Key Operations Using Jacobian
1. **Matrix-vector products**: `g(j) = dot_product(cg((l-1)*nf+1:l*nf), s)`
2. **Constraint evaluation**: Computing `A*s` where A is Jacobian
3. **QP subproblem**: Uses Jacobian columns in dual range-space method
4. **Active set updates**: Accessing individual constraint gradients

### Files & Routines Requiring Updates

#### Core Routines in `psqp_module.f90`:
- `compute_con_and_dcon` (line 661) - Builds/retrieves Jacobian
- `dual_range_space_quad_prog` (line 792) - QP solver, accesses `cg`
- `ops_after_constr_deletion` (line 1386) - Uses constraint vectors
- `psqp` main routine - Allocates and passes `cg` array

---

## Design for Sparse Support

### 1. Sparse Matrix Data Structure

Create new `psqp_sparse_module.f90`:

```fortran
module psqp_sparse_module
    use psqp_kind_module, only: wp => psqp_wp

    type :: sparse_matrix_csr
        !! Compressed Sparse Row format for constraint Jacobian
        integer :: nrows = 0              ! number of constraints (nc)
        integer :: ncols = 0              ! number of variables (nf)
        integer :: nnz = 0                ! number of non-zeros
        integer, allocatable :: row_ptr(:)     ! (nc+1) row start indices
        integer, allocatable :: col_ind(:)     ! (nnz) column indices
        real(wp), allocatable :: values(:)     ! (nnz) non-zero values
    contains
        procedure :: initialize => sparse_init
        procedure :: destroy => sparse_destroy
        procedure :: matvec => sparse_matvec       ! y = A*x
        procedure :: matvec_t => sparse_matvec_t   ! y = A^T*x
        procedure :: get_row => sparse_get_row     ! extract row
        procedure :: get_col => sparse_get_col     ! extract column
        procedure :: set_col => sparse_set_col     ! set column values
    end type sparse_matrix_csr

    type :: sparse_pattern
        !! Sparsity pattern (structure only, no values)
        integer :: nrows = 0
        integer :: ncols = 0
        integer :: nnz = 0
        integer, allocatable :: row_ptr(:)
        integer, allocatable :: col_ind(:)
    end type sparse_pattern

end module psqp_sparse_module
```

**Why CSR format?**
- Efficient row access (needed for constraint-by-constraint operations)
- Good for matrix-vector products
- Standard format, easy to convert from COO
- CSC could also work (would need column access instead)

### 2. API Extensions to `psqp_class`

```fortran
type :: psqp_class
    ! ... existing members ...

    ! Sparse mode support
    logical :: use_sparse = .false.  !! Use sparse Jacobian mode
    type(sparse_pattern), allocatable :: jac_pattern  !! User-provided sparsity
    type(sparse_matrix_csr), allocatable :: jac_sparse  !! Sparse Jacobian

    ! New callback for sparse mode
    procedure(sparse_jac_func), pointer :: sparse_jac => null()

contains
    ! ... existing procedures ...

    ! New procedures
    procedure :: set_jacobian_pattern
    procedure :: compute_sparse_jacobian
    procedure :: get_constraint_gradient_sparse
end type psqp_class
```

### 3. New User Interfaces

#### Option 1: User provides sparsity pattern
```fortran
! User tells solver which elements are nonzero
call solver%set_jacobian_pattern(row_indices, col_indices, nnz)

! Then provides sparse callback
abstract interface
    subroutine sparse_jac_func(me, nf, nc, x, row_ptr, col_ind, values)
        !! Compute full sparse Jacobian at once
        class(psqp_class), intent(inout) :: me
        integer, intent(in) :: nf, nc
        real(wp), intent(in) :: x(nf)
        integer, intent(in) :: row_ptr(nc+1)   ! sparsity structure
        integer, intent(in) :: col_ind(:)      ! column indices
        real(wp), intent(out) :: values(:)     ! fill in non-zero values
    end subroutine sparse_jac_func
end interface
```

#### Option 2: Auto-detect sparsity from dense callback
```fortran
! User provides dense gradient callback as before
! Solver detects zeros and builds sparse structure automatically
call solver%detect_sparsity(x0, tol=1.0e-14_wp)
```

#### Option 3: Sparse gradient callback per constraint (hybrid)
```fortran
! User provides indices of non-zeros for each constraint
abstract interface
    subroutine sparse_dcon_func(me, nf, kc, x, indices, values, nnz)
        class(psqp_class), intent(inout) :: me
        integer, intent(in) :: nf, kc
        real(wp), intent(in) :: x(nf)
        integer, intent(out) :: indices(:)  ! indices of non-zeros
        real(wp), intent(out) :: values(:)  ! non-zero values
        integer, intent(out) :: nnz         ! number of non-zeros
    end subroutine sparse_dcon_func
end interface
```

**Recommendation**: Start with Option 1 (explicit pattern) as it's most efficient.

### 4. Internal Modifications Required

#### File: `psqp_module.f90`

##### A. `psqpn` subroutine (line 147)
```fortran
! Current: allocates ra((nf+nc+8)*nf + 3*nc + 1) for dense Jacobian
! Need: conditional allocation based on sparse flag
! Dense: lcg = nf*nc
! Sparse: lcg = 0 (use separate sparse_matrix_csr object)
```

##### B. `compute_con_and_dcon` (line 661)
```fortran
! Current: fills cg array column by column
! Need: branch on sparse mode
if (me%use_sparse) then
    call compute_sparse_jacobian(me, nf, nc, x, ...)
else
    ! existing dense code
end if
```

##### C. `dual_range_space_quad_prog` (line 792)
All Jacobian accesses via `cg((inew-1)*nf+1)` need abstraction:
```fortran
! Line 865: call mxvdir(nf, cz(j), cg((kc-1)*nf+1), s, s)
! Replace with:
call get_constraint_gradient(me, kc, temp_gc)
call mxvdir(nf, cz(j), temp_gc, s, s)

! Or for sparse:
if (me%use_sparse) then
    ! sparse-specific operation
else
    ! dense operation
end if
```

##### D. Key locations needing abstraction:
- **Line 700**: `call mxvcop(nf, cg((kc-1)*nf+1), gc)` → get column kc
- **Line 703**: `call mxvcop(nf, gc, cg((kc-1)*nf+1))` → set column kc
- **Line 865**: `call mxvdir(nf, cz(j), cg((kc-1)*nf+1), s, s)` → axpy with column
- **Line 1071**: `gmax = mxvdot(nf, cg(jcg), s)` → dot product with column
- **Line 1084**: `g(j) = mxvdot(nf, cg((l-1)*nf+1), s)` → dot product
- **Line 2015**: `call mxvdir(nf, -cz(j), cg((l-1)*nf+1), g, g)` → axpy

#### File: `psqp_matrix_module.f90`
Add sparse matrix operations:
- `sparse_matvec` - sparse matrix-vector multiply
- `sparse_get_column` - extract column from CSR (requires transpose or CSC)
- `sparse_get_row` - extract row from CSR

### 5. Helper Routines Needed

```fortran
! In psqp_class:

subroutine get_constraint_gradient(me, kc, gc)
    !! Abstraction layer - get gradient of constraint kc
    class(psqp_class), intent(in) :: me
    integer, intent(in) :: kc
    real(wp), intent(out) :: gc(:)

    if (me%use_sparse) then
        ! Extract row kc from sparse Jacobian (CSR format - efficient)
        call me%jac_sparse%get_row(kc, gc)
    else
        ! Dense: copy from cg array
        gc = cg((kc-1)*nf+1:kc*nf)
    end if
end subroutine

real(wp) function constraint_gradient_dot(me, kc, v)
    !! Dot product of constraint kc gradient with vector v
    class(psqp_class), intent(in) :: me
    integer, intent(in) :: kc
    real(wp), intent(in) :: v(:)

    if (me%use_sparse) then
        ! Sparse row-vector dot product
        constraint_gradient_dot = me%jac_sparse%row_dot(kc, v)
    else
        ! Dense
        constraint_gradient_dot = mxvdot(nf, cg((kc-1)*nf+1), v)
    end if
end function
```

### 6. Storage Format Trade-offs

**CSR (Compressed Sparse Row)** - RECOMMENDED
- ✅ Efficient constraint-by-constraint access (row access)
- ✅ Good for matrix-vector products A*v
- ✅ Natural for constraint Jacobian (each row = one constraint)
- ❌ Column access requires search or transpose

**CSC (Compressed Sparse Column)**
- ✅ Efficient variable-by-variable access
- ✅ Natural for some operations
- ❌ Row access (constraint gradient) requires search
- ❌ Less natural for this application

**COO (Coordinate Format)**
- ✅ Easy to construct
- ✅ Easy to convert to CSR/CSC
- ❌ Inefficient for operations
- 💡 Use as intermediate format

**Hybrid: CSR + Column cache**
- Store in CSR for constraint access
- Cache frequently-used columns in dense format
- Best performance but more complex

### 7. Backward Compatibility Strategy

```fortran
! Default behavior: dense mode (no user changes needed)
type(psqp_class) :: solver
call solver%psqpn(...)  ! works as before

! Opt-in to sparse mode:
type(psqp_class) :: solver
call solver%set_jacobian_pattern(rows, cols, nnz)
solver%sparse_jac => my_sparse_jac_func
call solver%psqpn(...)  ! automatically uses sparse mode
```

### 8. Testing Strategy

1. **Validation**: Run all existing dense tests in sparse mode
   - Convert dense Jacobian to sparse
   - Verify identical results

2. **Sparse-specific tests**:
   - Large sparse problems (CUTEst)
   - Verify memory reduction
   - Benchmark performance

3. **Edge cases**:
   - Very sparse (1% fill)
   - Moderately sparse (20% fill)
   - Nearly dense (80% fill) - should use dense
   - Diagonal constraints
   - Block-structured Jacobians

---

## Implementation Phases

### Phase 1: Infrastructure (1-2 days)
- [ ] Create `psqp_sparse_module.f90`
- [ ] Implement `sparse_matrix_csr` type
- [ ] Implement sparse matrix operations
- [ ] Add unit tests for sparse matrix ops

### Phase 2: API Extensions (1 day)
- [ ] Add sparse mode flag to `psqp_class`
- [ ] Add `set_jacobian_pattern` method
- [ ] Define `sparse_jac_func` interface
- [ ] Add helper methods (get_constraint_gradient, etc.)

### Phase 3: Core Algorithm Updates (2-3 days)
- [ ] Abstract Jacobian access in `dual_range_space_quad_prog`
- [ ] Update `compute_con_and_dcon` for sparse mode
- [ ] Modify memory allocation in `psqpn`
- [ ] Update all Jacobian accesses to use abstraction

### Phase 4: Testing & Validation (1-2 days)
- [ ] Create sparse test problems
- [ ] Validate against dense mode
- [ ] Performance benchmarking
- [ ] Documentation

### Phase 5: Optimizations (optional)
- [ ] Auto-detect sparsity
- [ ] Hybrid sparse/dense for nearly dense
- [ ] Column caching for frequently accessed
- [ ] Parallel sparse operations

---

## Example Usage

```fortran
program sparse_example
    use psqp_module
    use psqp_sparse_module

    type(psqp_class) :: solver
    integer, parameter :: nf = 1000, nc = 500

    ! Define sparsity: each constraint depends on ~10 variables
    integer :: rows(5000), cols(5000)
    ! ... fill in sparsity pattern ...

    ! Configure sparse mode
    call solver%set_jacobian_pattern(rows, cols, 5000)

    ! Solve (uses sparse Jacobian automatically)
    call solver%psqpn(nf, nb, nc, x, ix, xl, xu, cf, ic, cl, cu, &
                      ipar, rpar, f, cmax, gmax, iprnt, iterm, &
                      obj, dobj, con, dcon_sparse)

contains

    subroutine dcon_sparse(me, nf, nc, x, row_ptr, col_ind, values)
        ! User fills in sparse Jacobian values
        class(psqp_class), intent(inout) :: me
        integer, intent(in) :: nf, nc
        real(wp), intent(in) :: x(nf)
        integer, intent(in) :: row_ptr(nc+1), col_ind(:)
        real(wp), intent(out) :: values(:)

        ! Fill values based on col_ind structure
        ! ...
    end subroutine

end program
```

---

## Memory Usage Comparison

**Dense mode**: `nf × nc × 8 bytes`
- Example: 1000 vars × 500 constraints = 4 MB

**Sparse mode**: `nnz × 12 bytes` (value + index)
- Example: 5000 non-zeros = 60 KB
- **Savings: 98.5%**

For very large problems (nf=10000, nc=5000, 1% fill):
- Dense: 400 MB
- Sparse: 6 MB
- **Savings: 98.5%**

---

## Performance Considerations

**When sparse is faster:**
- Fill ratio < 20%
- Large problems (nf > 500, nc > 200)
- Structured sparsity (banded, block)

**When dense may be competitive:**
- Fill ratio > 50%
- Small problems (overhead dominates)
- Dense BLAS available

**Recommendation**: Provide auto-detection based on fill ratio
```fortran
if (nnz > 0.3 * nf * nc) then
    ! Use dense mode even if sparse pattern provided
    use_sparse = .false.
end if
```

---

## Open Questions

1. **CSR vs CSC**: Start with CSR (constraint-centric), add CSC later if needed?
2. **Column access**: Cache columns or compute on-demand?
3. **User interface**: Require explicit pattern or auto-detect?
4. **Parallel**: Add OpenMP for sparse matvec? (future enhancement)
5. **External libraries**: Link to SPARSKIT/SuiteSparse or implement from scratch?

**Recommendations**:
1. CSR format (constraint-centric is natural)
2. Compute columns on-demand initially
3. Explicit pattern required (performance-oriented users)
4. Keep simple for initial implementation
5. Pure Fortran initially (no external dependencies)
