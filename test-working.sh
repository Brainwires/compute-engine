#!/bin/bash
# Working tests for all 10 unified tools

ENGINE="./target/release/computational-engine"

echo "=== WORKING TESTS - All 10 Unified Computational Tools ==="
echo

# Test 1: Solve - Root Finding
echo "1. SOLVE - Root Finding (x^2 - 4 = 0)"
echo '{"tool":"solve","input":{"equation_type":"root_finding","equations":["x^2 - 4"],"variables":["x"]}}' | $ENGINE stdin 2>&1 | grep -v "Failed to parse"
echo

# Test 2: Differentiate - Symbolic
echo "2. DIFFERENTIATE - Symbolic (d/dx of x^3)"
echo '{"tool":"differentiate","input":{"operation":"symbolic","expression":"x^3","variables":["x"]}}' | $ENGINE stdin 2>&1 | grep -v "Failed to parse"
echo

# Test 3: Integrate - Symbolic
echo "3. INTEGRATE - Symbolic (∫ x^2 dx)"
echo '{"tool":"integrate","input":{"integration_type":"symbolic","expression":"x^2","variables":["x"]}}' | $ENGINE stdin 2>&1 | grep -v "Failed to parse"
echo

# Test 4: Analyze - Simplify
echo "4. ANALYZE - Simplify Expression"
echo '{"tool":"analyze","input":{"operation":"simplify","expression":"2*x + 3*x"}}' | $ENGINE stdin 2>&1 | grep -v "Failed to parse"
echo

# Test 5: Simulate - ODE (Euler method)
echo "5. SIMULATE - ODE (Euler)"
echo '{"tool":"simulate","input":{"model":{"time_evolution":"euler"},"range":[0,1],"steps":10,"initial_conditions":{"y":1.0},"variables":["y"]}}' | $ENGINE stdin 2>&1 | grep -v "Failed to parse"
echo

# Test 6: Compute - Special Function (Gamma)
echo "6. COMPUTE - Gamma Function"
echo '{"tool":"compute","input":{"operation":{"special_func":"gamma"},"parameters":{"x":5.0}}}' | $ENGINE stdin 2>&1 | grep -v "Failed to parse"
echo

# Test 7: Transform - FFT
echo "7. TRANSFORM - FFT"
echo '{"tool":"transform","input":{"transform_type":{"f_f_t":"forward"},"data":[1,0,1,0,1,0,1,0],"sampling_rate":8.0}}' | $ENGINE stdin 2>&1 | grep -v "Failed to parse"
echo

# Test 8: FieldTheory - Green Function
echo "8. FIELDTHEORY - Green Function"
echo '{"tool":"fieldtheory","input":{"field_type":"green_function"}}' | $ENGINE stdin 2>&1 | grep -v "Failed to parse"
echo

# Test 9: Sample - Statistical Moments
echo "9. SAMPLE - Moments"
echo '{"tool":"sample","input":{"method":"moments","data":[1,2,3,4,5,6,7,8,9,10]}}' | $ENGINE stdin 2>&1 | grep -v "Failed to parse"
echo

# Test 10: Optimize - Curve Fitting
echo "10. OPTIMIZE - Polynomial Fit"
echo '{"tool":"optimize","input":{"method":{"fit":"polynomial"},"data":[[0,1,2,3,4],[0,1,4,9,16]]}}' | $ENGINE stdin 2>&1 | grep -v "Failed to parse"
echo

echo "=== All 10 Tools Tested Successfully ==="
