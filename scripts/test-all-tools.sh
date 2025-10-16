#!/bin/bash
# Test all 10 unified tools

ENGINE="./target/release/brainwires-compute-engine"

echo "=== Testing All 10 Unified Computational Tools ==="
echo

# Test 1: Solve
echo "1. SOLVE - Linear System"
echo '{"tool":"solve","input":{"equation_type":"linear_system","equations":["x + y = 5","x - y = 1"],"variables":["x","y"]}}' | $ENGINE stdin 2>&1
echo

# Test 2: Differentiate
echo "2. DIFFERENTIATE - Symbolic"
echo '{"tool":"differentiate","input":{"operation":"symbolic","expression":"x^2 + y^2","variables":["x"]}}' | $ENGINE stdin 2>&1
echo

# Test 3: Integrate
echo "3. INTEGRATE - Symbolic"
echo '{"tool":"integrate","input":{"integration_type":"symbolic","expression":"x^2","variable":"x"}}' | $ENGINE stdin 2>&1
echo

# Test 4: Analyze
echo "4. ANALYZE - Simplify"
echo '{"tool":"analyze","input":{"operation":"simplify","expression":"(x + 1)^2 - (x - 1)^2"}}' | $ENGINE stdin 2>&1
echo

# Test 5: Simulate
echo "5. SIMULATE - Brownian Motion"
echo '{"tool":"simulate","input":{"model":{"stochastic":"brownian_motion"},"range":[0,1],"steps":100,"parameters":{"drift":0.0,"volatility":1.0}}}' | $ENGINE stdin 2>&1
echo

# Test 6: Compute
echo "6. COMPUTE - Matrix Determinant"
echo '{"tool":"compute","input":{"operation":{"matrix":"determinant"},"data":[[1,2],[3,4]]}}' | $ENGINE stdin 2>&1
echo

# Test 7: Transform
echo "7. TRANSFORM - FFT"
echo '{"tool":"transform","input":{"transform_type":{"fft":"forward"},"data":[1,2,3,4,5,6,7,8]}}' | $ENGINE stdin 2>&1
echo

# Test 8: FieldTheory
echo "8. FIELDTHEORY - EM Antenna"
echo '{"tool":"fieldtheory","input":{"field_type":{"em":"antenna"},"frequency":2.4e9}}' | $ENGINE stdin 2>&1
echo

# Test 9: Sample
echo "9. SAMPLE - Moments"
echo '{"tool":"sample","input":{"method":"moments","data":[1,2,3,4,5,6,7,8,9,10]}}' | $ENGINE stdin 2>&1
echo

# Test 10: Optimize
echo "10. OPTIMIZE - Polynomial Fit"
echo '{"tool":"optimize","input":{"method":{"fit":"polynomial"},"data":[[0,1,2,3,4],[0,1,4,9,16]]}}' | $ENGINE stdin 2>&1
echo

echo "=== All 10 Tools Tested ==="
