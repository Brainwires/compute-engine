#!/usr/bin/env python3
"""
Properly add EnumIter and Default to ALL enums with #[default] on FIRST variant.
"""

import re

# Read the file
with open('src/engine/equations.rs', 'r') as f:
    content = f.read()

# Step 1: Add EnumIter import (if not already there)
if 'use strum::EnumIter;' not in content:
    content = content.replace(
        'use serde::{Deserialize, Deserializer, Serialize};',
        'use serde::{Deserialize, Deserializer, Serialize};\nuse strum::EnumIter;'
    )

# Step 2: Add EnumIter and Default to all enum derives (except SelectionCriteria)
lines = content.split('\n')
output = []
i = 0

while i < len(lines):
    line = lines[i]

    # Check if this is a derive line before an enum
    if line.strip().startswith('#[derive(') and i + 2 < len(lines):
        # Check next couple lines to see if it's an enum definition
        next_lines = '\n'.join(lines[i:i+5])

        if 'pub enum' in next_lines and 'SelectionCriteria' not in next_lines:
            # Add EnumIter and Default if not present
            if 'EnumIter' not in line:
                line = line.replace(')]', ', EnumIter, Default)]')

            output.append(line)

            # Skip to the line with "pub enum"
            i += 1
            while i < len(lines) and 'pub enum' not in lines[i]:
                output.append(lines[i])
                i += 1

            if i < len(lines):
                enum_line = lines[i]
                output.append(enum_line)  # pub enum Name {

                # Next line should be the opening brace or first variant
                i += 1
                if i < len(lines):
                    bracket_or_variant = lines[i]

                    if '{' in bracket_or_variant and '{' in enum_line:
                        # { was on same line as pub enum, current line is first variant
                        variant_line = bracket_or_variant
                    else:
                        # { is on its own line
                        output.append(bracket_or_variant)
                        i += 1
                        if i < len(lines):
                            variant_line = lines[i]
                        else:
                            continue

                    # Add #[default] before first variant if not already there
                    if variant_line.strip() and not variant_line.strip().startswith('#['):
                        indent = len(variant_line) - len(variant_line.lstrip())
                        output.append(' ' * indent + '#[default]\n')

                    output.append(variant_line)
        else:
            output.append(line)
    else:
        output.append(line)

    i += 1

# Write back
with open('src/engine/equations.rs', 'w') as f:
    f.write('\n'.join(output))

print("✓ Added EnumIter, Default derives and #[default] attributes properly")
