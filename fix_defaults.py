#!/usr/bin/env python3
"""
Add Default derive and #[default] attribute to first variant of each enum.
"""

import re

# Read the file
with open('src/engine/equations.rs', 'r') as f:
    content = f.read()

# Step 1: Add Default to all EnumIter derives (except SelectionCriteria)
content = re.sub(
    r'#\[derive\(Debug, Clone, Serialize, Deserialize, EnumIter\)\]\s*\n#\[serde\(rename_all = "snake_case"\)\]\s*\npub enum (?!SelectionCriteria)',
    r'#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]\n#[serde(rename_all = "snake_case")]\npub enum ',
    content
)

# Step 2: Add #[default] to first variant after each enum definition
# Pattern: find "pub enum Name {\n    FirstVariant"
def add_default_to_first_variant(match):
    enum_def = match.group(0)
    # Add #[default] before the first variant
    result = re.sub(
        r'(\{\s*\n\s+)(\w+)',  # Find first variant after opening brace
        r'\1#[default]\n    \2',
        enum_def,
        count=1  # Only replace first occurrence
    )
    return result

# Find all enum definitions and add #[default] to first variant
pattern = r'pub enum \w+ \{[^}]+?\n    \w+'
content = re.sub(pattern, add_default_to_first_variant, content, flags=re.MULTILINE)

# Write back
with open('src/engine/equations.rs', 'w') as f:
    f.write(content)

print("✓ Added Default derives and #[default] attributes")
