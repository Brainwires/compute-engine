#!/usr/bin/env python3
"""
Add EnumIter to ALL enums in equations.rs, regardless of current derive state.
"""

import re

# Read the file
with open('src/engine/equations.rs', 'r') as f:
    lines = f.readlines()

output = []
i = 0

while i < len(lines):
    line = lines[i]

    # Check if this line has "pub enum"
    if line.strip().startswith('pub enum '):
        enum_name = line.strip().split()[2].rstrip('{')

        # Skip SelectionCriteria (has manual Default impl)
        if enum_name == 'SelectionCriteria':
            output.append(line)
            i += 1
            continue

        # Look back to find/add #[derive(...)]
        # Check if previous non-empty line is #[derive(...)]
        j = i - 1
        while j >= 0 and (lines[j].strip() == '' or lines[j].strip().startswith('#[serde') or lines[j].strip().startswith('#[strum')):
            j -= 1

        if j >= 0 and lines[j].strip().startswith('#[derive('):
            # Has a derive line - make sure it has EnumIter
            derive_line_idx = j
            derive_line = lines[derive_line_idx].strip()

            if 'EnumIter' not in derive_line:
                # Add EnumIter to existing derives
                lines[derive_line_idx] = derive_line.replace(')]', ', EnumIter)]') + '\n'

            # Output all lines from derive to current
            for k in range(derive_line_idx, i):
                if k not in [x for x in range(len(output))]:  # Avoid duplicates
                    output.append(lines[k])
            output.append(line)
        else:
            # No derive line - add one
            indent = len(line) - len(line.lstrip())
            output.append(' ' * indent + '#[derive(Debug, Clone, Serialize, Deserialize, EnumIter)]\n')

            # Add back any #[serde] or #[strum] lines that came before
            for k in range(j + 1, i):
                output.append(lines[k])

            output.append(line)
        i += 1
    else:
        # Not an enum line, just add it
        # But skip if we already added it in the enum processing above
        if i == 0 or not lines[i-1].strip().startswith('pub enum'):
            output.append(line)
        i += 1

# Write back
with open('src/engine/equations.rs', 'w') as f:
    f.writelines(output)

print("✓ Added EnumIter to all enums")
