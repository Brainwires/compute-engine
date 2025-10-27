#!/usr/bin/env python3
"""
Add #[default] attribute to first variant of each enum (except SelectionCriteria).
"""

import re

# Read the file
with open('src/engine/equations.rs', 'r') as f:
    lines = f.readlines()

output = []
i = 0
while i < len(lines):
    line = lines[i]
    output.append(line)

    # Check if this is an enum definition line (but not SelectionCriteria)
    if line.strip().startswith('pub enum ') and 'SelectionCriteria' not in line:
        # Next line should be opening brace or might be on same line
        i += 1
        if i < len(lines):
            output.append(lines[i])  # Add the { line

            # Next line should be first variant - add #[default] before it
            i += 1
            if i < len(lines):
                variant_line = lines[i]
                # Check if it's a variant (not already #[default])
                if variant_line.strip() and not variant_line.strip().startswith('#[') and not variant_line.strip().startswith('}'):
                    # Get indentation
                    indent = len(variant_line) - len(variant_line.lstrip())
                    # Add #[default] with same indentation
                    output.append(' ' * indent + '#[default]\n')
                output.append(variant_line)

    i += 1

# Write back
with open('src/engine/equations.rs', 'w') as f:
    f.writelines(output)

print("✓ Added #[default] attributes to first variant of each enum")
