#!/usr/bin/env python3
"""
Automatically add EnumIter, Display, EnumString derives to all enums in equations.rs
and add strum documentation attributes.
"""

import re
import sys

# Read the file
with open('src/engine/equations.rs', 'r') as f:
    content = f.read()

# Pattern to find enum definitions
enum_pattern = r'#\[derive\(([^\]]+)\)\]\s+#\[serde\(rename_all = "snake_case"\)\]\s+pub enum (\w+) \{'

def add_strum_derives(match):
    derives = match.group(1)
    enum_name = match.group(2)

    # Skip if already has EnumIter
    if 'EnumIter' in derives:
        return match.group(0)

    # Add strum derives
    new_derives = f"{derives}, EnumIter, Display, EnumString"

    result = f'#[derive({new_derives})]\n'
    result += '#[serde(rename_all = "snake_case")]\n'
    result += '#[strum(serialize_all = "snake_case")]\n'
    result += f'pub enum {enum_name} {{'

    return result

# Replace all enum definitions
content = re.sub(enum_pattern, add_strum_derives, content)

# Write back
with open('src/engine/equations.rs', 'w') as f:
    f.write(content)

print("✓ Added strum derives to all enums")
