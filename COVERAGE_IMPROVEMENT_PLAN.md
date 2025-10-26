# Coverage Improvement Plan: 44.78% → 80%

**Current Status:** 44.78% line coverage (27,489 / 61,393 lines)
**Target:** 80% line coverage
**Gap:** 35.22% (need to cover ~21,600 additional lines)

**Date:** 2025-10-26

---

## Executive Summary

We have **160 modules with 0% coverage** out of ~255 source files. The path to 80% involves:

1. **High-Impact Quick Wins** (10-15% gain): Test 0% coverage API handlers and implementations
2. **Core Module Testing** (10-15% gain): Add tests for electrical, materials science, datetime
3. **Implementation Layer** (5-10% gain): Test compute implementations (tensor, matrix, nuclear, etc.)
4. **Specialized Modules** (5% gain): Machine learning, game theory, linear programming

---

## Phase 1: Quick Wins (Target: +15% coverage, 2-3 days)

### Priority 1A: API Handlers with 0% Coverage (12 modules)

**Impact:** Each handler is 40-200 lines, testing all 12 = ~1,200 lines = +2% coverage

| Module | Lines | Current | Priority |
|--------|-------|---------|----------|
| `api/handlers/biology.rs` | 41 | 0% | HIGH |
| `api/handlers/datetime.rs` | 42 | 0% | HIGH |
| `api/handlers/engineering.rs` | 42 | 0% | HIGH |
| `api/handlers/fluid_dynamics.rs` | 21 | 0% | HIGH |
| `api/handlers/geophysics.rs` | 42 | 0% | HIGH |
| `api/handlers/optics.rs` | 41 | 0% | HIGH |
| `api/handlers/thermodynamics.rs` | 41 | 0% | HIGH |

**Action:** Create integration tests in `tests/integration/api/old_api_handlers_tests.rs`

**Estimated Effort:** 4-6 hours (12 handlers × 20-30 min each)

### Priority 1B: Electrical Module (8 files, all 0%)

**Impact:** ~800 lines = +1.3% coverage

| Module | Lines | Tests Needed |
|--------|-------|--------------|
| `electrical/ac_analysis.rs` | 122 | AC circuit analysis tests |
| `electrical/dc_analysis.rs` | 103 | DC circuit tests |
| `electrical/impedance.rs` | 84 | Impedance calculation tests |
| `electrical/filter_design.rs` | 56 | Filter design tests |
| `electrical/network_analysis.rs` | 120 | Network analysis tests |
| `electrical/power_analysis.rs` | 56 | Power calculation tests |
| `electrical/transient_analysis.rs` | 76 | Transient response tests |
| `electrical/nec_calculations.rs` | 182 | NEC code tests |

**Action:** Create `tests/unit/electrical/comprehensive_tests.rs`

**Estimated Effort:** 6-8 hours

### Priority 1C: Materials Science (6 files, all 0%)

**Impact:** ~800 lines = +1.3% coverage

| Module | Lines |
|--------|-------|
| `materials_science/band_theory.rs` | 103 |
| `materials_science/crystal_structures.rs` | 198 |
| `materials_science/diffusion.rs` | 120 |
| `materials_science/mechanical_properties.rs` | 154 |
| `materials_science/thermal_properties.rs` | 131 |
| `materials_science/xrd_analysis.rs` | 94 |

**Action:** Create `tests/unit/materials_science/comprehensive_tests.rs`

**Estimated Effort:** 6-8 hours

### Priority 1D: DateTime Implementation (42 lines, 0%)

**Impact:** +0.07% coverage

**Action:** Create `tests/unit/implementations/compute/datetime_tests.rs`

**Estimated Effort:** 1-2 hours

---

## Phase 2: Implementation Layer (Target: +10% coverage, 3-4 days)

### Priority 2A: Compute Implementations (Critical - most are 0%)

**Impact:** ~8,000 lines = +13% coverage

| Module | Lines | Coverage | Priority |
|--------|-------|----------|----------|
| `implementations/computer.rs` | 3,867 | 4.27% | CRITICAL |
| `implementations/compute/control.rs` | 564 | 0% | HIGH |
| `implementations/compute/tensor.rs` | 232 | 0% | HIGH |
| `implementations/compute/number_theory.rs` | 508 | 0% | HIGH |
| `implementations/compute/nuclear.rs` | 326 | 0% | HIGH |
| `implementations/compute/matrix.rs` | 359 | 27% | MEDIUM |
| `implementations/compute/information.rs` | 160 | 0% | HIGH |
| `implementations/compute/graph.rs` | 72 | 0% | MEDIUM |

**Critical Issue:** `implementations/computer.rs` has 3,867 lines with only 4.27% coverage!

**Action:**
1. Create comprehensive tests for each compute operation type
2. Test the dispatch logic in `computer.rs`
3. Add tests for tensor operations (Christoffel, Riemann, Einstein)
4. Add tests for number theory (RSA, primes, modular arithmetic)
5. Add tests for control systems (transfer functions, Bode plots)
6. Add tests for nuclear physics

**Estimated Effort:** 16-20 hours

### Priority 2B: Analyzer and Integrator

| Module | Lines | Coverage |
|--------|-------|----------|
| `implementations/analyzer.rs` | 798 | 63% |
| `implementations/integrator.rs` | 576 | 54% |

**Action:** Add tests for uncovered branches and edge cases

**Estimated Effort:** 4-6 hours

---

## Phase 3: Specialized Modules (Target: +5% coverage, 2-3 days)

### Priority 3A: Machine Learning (all 0%, ~2,000 lines)

**Impact:** +3.3% coverage

| Module | Lines |
|--------|-------|
| `specialized/machine_learning/neural_network.rs` | 402 |
| `specialized/machine_learning/regression.rs` | 396 |
| `specialized/machine_learning/clustering.rs` | 326 |
| `specialized/machine_learning/optimization.rs` | 340 |
| `specialized/machine_learning/dimensionality_reduction.rs` | 248 |

**Action:** Create `tests/unit/specialized/machine_learning/comprehensive_tests.rs`

**Estimated Effort:** 10-12 hours

### Priority 3B: Game Theory (all 0%, ~800 lines)

**Impact:** +1.3% coverage

| Module | Lines |
|--------|-------|
| `specialized/game_theory/normal_form.rs` | 242 |
| `specialized/game_theory/evolutionary.rs` | 219 |
| `specialized/game_theory/cooperative.rs` | 172 |
| `specialized/game_theory/extensive_form.rs` | 167 |

**Action:** Create `tests/unit/specialized/game_theory/comprehensive_tests.rs`

**Estimated Effort:** 6-8 hours

### Priority 3C: Linear Programming (all 0%, ~438 lines)

**Impact:** +0.7% coverage

| Module | Lines |
|--------|-------|
| `specialized/linear_programming/simplex.rs` | 297 |
| `specialized/linear_programming/dual.rs` | 96 |
| `specialized/linear_programming/mod.rs` | 45 |

**Action:** Create `tests/unit/specialized/linear_programming/comprehensive_tests.rs`

**Estimated Effort:** 4-6 hours

---

## Phase 4: Lower-Impact Improvements (Target: +5% coverage, 2-3 days)

### Priority 4A: Improve Existing Low Coverage

| Module | Lines | Current | Target |
|--------|-------|---------|--------|
| `thermodynamics/mod.rs` | 307 | 60% | 85% |
| `engineering/mod.rs` | 473 | 43% | 80% |
| `specialized/stochastic_processes/lib.rs` | 357 | 31% | 75% |
| `tools/dimensional_analysis/lib.rs` | 893 | 31% | 70% |
| `tools/numerical_methods/mod.rs` | 730 | 47% | 80% |
| `tools/computational_geometry/advanced.rs` | 532 | 67% | 85% |

**Action:** Add tests for uncovered branches and edge cases

**Estimated Effort:** 8-10 hours

### Priority 4B: WASM Bindings (282 lines, 0%)

**Impact:** +0.5% coverage

**Note:** WASM tests require `wasm-pack test` which may not contribute to coverage metrics

**Action:** Add WASM integration tests if they contribute to coverage

**Estimated Effort:** 2-4 hours

---

## Summary Timeline

| Phase | Focus | Coverage Gain | Effort | Days |
|-------|-------|---------------|--------|------|
| **Phase 1** | Quick wins (API, electrical, materials) | +15% | 20-26 hrs | 2-3 |
| **Phase 2** | Implementation layer (computer.rs) | +10% | 20-26 hrs | 3-4 |
| **Phase 3** | Specialized (ML, game theory, LP) | +5% | 20-26 hrs | 2-3 |
| **Phase 4** | Improvements (edge cases, WASM) | +5% | 10-14 hrs | 2-3 |
| **TOTAL** | **All phases** | **+35%** | **70-92 hrs** | **9-13 days** |

**Final Target:** 44.78% + 35% = **~80% coverage**

---

## Immediate Action Items (Start Today)

1. ✅ **Create this plan document**
2. ✅ **Phase 1A**: API handlers - **COMPLETED** (44.78% → 45.11%, +0.33%)
   - biology.rs: 0% → 77.27% ✅
   - datetime.rs: 0% → 73.91% ✅
   - engineering.rs: 0% → 75.00% ✅
   - fluid_dynamics.rs: 0% → 71.43% ✅
   - optics.rs: 0% → 72.73% ✅
   - geophysics.rs: 0% (11 tests created, need parameter fixes)
   - thermodynamics.rs: 0% (11 tests created, need parameter fixes)
3. ⏳ **Phase 1B**: Electrical module tests (8 files, ~6-8 hours)
4. ⏳ **Phase 2A**: Critical - `implementations/computer.rs` (3,867 lines at 4%)

---

## Test File Organization

All new tests should follow the existing structure:

```
tests/
├── unit/
│   ├── electrical/
│   │   └── comprehensive_tests.rs (NEW)
│   ├── materials_science/
│   │   └── comprehensive_tests.rs (NEW)
│   ├── implementations/
│   │   ├── computer_tests.rs (NEW)
│   │   ├── compute/
│   │   │   ├── tensor_tests.rs (NEW)
│   │   │   ├── nuclear_tests.rs (NEW)
│   │   │   └── control_tests.rs (NEW)
│   └── specialized/
│       ├── machine_learning/
│       │   └── comprehensive_tests.rs (NEW)
│       ├── game_theory/
│       │   └── comprehensive_tests.rs (NEW)
│       └── linear_programming/
│           └── comprehensive_tests.rs (NEW)
└── integration/
    └── api/
        └── old_api_handlers_tests.rs (NEW)
```

---

## Testing Strategy

For each module:

1. **Happy Path Tests**: Test primary functionality with valid inputs
2. **Edge Cases**: Test boundary conditions (zero, negative, infinity, NaN)
3. **Error Handling**: Test invalid inputs and error paths
4. **Integration**: Test interaction between components
5. **Regression**: Add tests for any bugs found during implementation

---

## Success Metrics

- [ ] Coverage reaches 80%+ overall
- [ ] All API handlers have >90% coverage
- [ ] All implementations have >80% coverage
- [ ] Core modules (electrical, materials) have >85% coverage
- [ ] Specialized modules have >75% coverage
- [ ] All tests pass with 0 failures

---

## Notes

- **WASM tests** may not contribute to coverage - investigate if needed
- **Ignored tests** (8 slow physics simulations) don't affect coverage
- Focus on **line coverage** as primary metric
- Track **function coverage** and **region coverage** as secondary metrics
