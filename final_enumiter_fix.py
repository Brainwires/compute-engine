#!/usr/bin/env python3
"""
Final correct implementation:
1. Add EnumIter to ALL enums
2. Add Default + #[default] ONLY to unit-variant enums (no associated data)
3. Do NOT add Default to enums with variants that have associated data
"""

import re

with open('src/engine/equations.rs', 'r') as f:
    content = f.read()

# Step 1: Add strum import
if 'use strum::EnumIter;' not in content:
    content = content.replace(
        'use serde::{Deserialize, Deserializer, Serialize};',
        'use serde::{Deserialize, Deserializer, Serialize};\nuse strum::EnumIter;'
    )

lines = content.split('\n')
output = []
i = 0

while i < len(lines):
    line = lines[i]

    # Look for enum definitions
    if line.strip().startswith('pub enum '):
        enum_name = line.strip().split()[2].rstrip('{')

        # Find the derive line (go backwards)
        derive_idx = i - 1
        while derive_idx >= 0 and (lines[derive_idx].strip() == '' or
                                    lines[derive_idx].strip().startswith('#[serde') or
                                    lines[derive_idx].strip().startswith('#[strum')):
            derive_idx -= 1

        has_derive = derive_idx >= 0 and lines[derive_idx].strip().startswith('#[derive(')

        # Scan forward to check if enum has variants with associated data
        j = i + 1
        has_assoc_data = False
        is_unit_only = True
        while j < len(lines) and lines[j].strip() != '}':
            variant_line = lines[j].strip()
            if variant_line and not variant_line.startswith('//') and not variant_line.startswith('#['):
                # This is a variant
                if '(' in variant_line:
                    has_assoc_data = True
                    is_unit_only = False
                    break
            j += 1

        # Now output the enum with correct derives
        if has_derive:
            # Modify existing derive line
            derive_line = lines[derive_idx]
            if 'EnumIter' not in derive_line:
                derive_line = derive_line.rstrip(']\n')
                if is_unit_only and enum_name != 'SelectionCriteria':
                    derive_line += ', EnumIter, Default)]\n'
                else:
                    derive_line += ', EnumIter)]\n'

            # Output everything from derive line to enum line
            output.append(derive_line)
            for k in range(derive_idx + 1, i):
                output.append(lines[k])
            output.append(line)

            # For unit-only enums, add #[default] to first variant
            if is_unit_only and enum_name != 'SelectionCriteria':
                i += 1
                # Skip opening brace if on next line
                if i < len(lines) and lines[i].strip() == '{':
                    output.append(lines[i])
                    i += 1

                # Add #[default] before first variant
                if i < len(lines):
                    variant_line = lines[i]
                    if variant_line.strip() and not variant_line.strip().startswith('#['):
                        indent = len(variant_line) - len(variant_line.lstrip())
                        output.append(' ' * indent + '#[default]\n')
                    output.append(variant_line)
            i += 1
            continue
        else:
            # No derive line - add one
            if is_unit_only and enum_name != 'SelectionCriteria':
                output.append('#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]\n')
            else:
                output.append('#[derive(Debug, Clone, Serialize, Deserialize, EnumIter)]\n')

            # Add #[serde] / #[strum] lines
            for k in range(derive_idx + 1, i):
                output.append(lines[k])
            output.append(line)

            # For unit-only enums, add #[default] to first variant
            if is_unit_only and enum_name != 'SelectionCriteria':
                i += 1
                if i < len(lines) and lines[i].strip() == '{':
                    output.append(lines[i])
                    i += 1

                if i < len(lines):
                    variant_line = lines[i]
                    if variant_line.strip() and not variant_line.strip().startswith('#['):
                        indent = len(variant_line) - len(variant_line.lstrip())
                        output.append(' ' * indent + '#[default]\n')
                    output.append(variant_line)
            i += 1
            continue

    output.append(line)
    i += 1

# Write back
with open('src/engine/equations.rs', 'w') as f:
    f.write('\n'.join(output))

print("✓ Added EnumIter correctly to all enums with Default only on unit enums")
