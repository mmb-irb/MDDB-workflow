#!/usr/bin/env python3

import re
import os
from pathlib import Path


def extract_yaml_documentation(yaml_file_path):
    """Extract documentation from YAML file comments."""
    with open(yaml_file_path, 'r') as f:
        content = f.read()

    lines = content.split('\n')

    # Extract file description from the header
    file_description = []
    i = 0
    # Skip document start marker
    while i < len(lines) and lines[i].strip() in ['---', '']:
        i += 1

    # Get file description comments at the top
    while i < len(lines) and lines[i].strip().startswith('#'):
        line = lines[i].strip('# \n')
        if line and not line.startswith('-'):  # Skip comment dividers
            file_description.append(line)
        i += 1
        if i < len(lines) and not lines[i].strip():  # Stop at blank line
            break

    # Process sections
    sections = []
    i = 0

    while i < len(lines):
        line = lines[i]

        # Detect section header pattern: #----...----
        if re.match(r'^#-+$', line.strip()):
            section_title = None
            section_overview = []
            section_content = []
            # The section title is always on the line immediately after the divider
            section_title = lines[i + 1].strip('# \n')
            # Skip to end of header block (including the title line)
            i += 3

            # Collect section overview (only general comments, not field-specific ones)
            while i < len(lines):
                line = lines[i]
                i += 1
                if not line.strip(): # Stop on blank line
                    break
                section_overview.append(line.strip('# \n'))

            # Collect section content until next section or end
            while i < len(lines):
                line = lines[i]
                if re.match(r'^#-+$', line.strip()):
                    # Found next section, break
                    break
                section_content.append(line)
                i += 1

            sections.append((section_title, section_overview, section_content))
        else:
            i += 1

    # Generate RST content
    rst_content = ".. _input_file_documentation: generated with generate_input_docs.py\n\n"
    rst_content += "Input File\n"
    rst_content += "==========================\n\n"
    rst_content += ".. note::\n"
    rst_content += "   This documentation is automatically generated from the YAML template file.\n\n"
    
    # Add file description
    if file_description:
        rst_content += "\n".join(file_description) + "\n\n"

    # Process each section
    for title, overview, content in sections:
        rst_content += f"{title.capitalize()}\n"
        rst_content += "-" * len(title) + "\n\n"

        # Add section overview
        if overview:
            rst_content += "\n".join(overview) + "\n\n"

        # Process fields in this section
        i = 0
        while i < len(content):
            # Skip blank lines
            while i < len(content) and not content[i].strip():
                i += 1

            if i >= len(content):
                break

            # Collect field description comments
            field_description = []
            while i < len(content) and content[i].strip().startswith('#'):
                comment = content[i].strip('# \n')
                field_description.append(comment)
                i += 1

            # Skip blank lines after comments
            while i < len(content) and not content[i].strip():
                i += 1

            if i >= len(content):
                break

            # Check for field definition
            line = content[i]
            if ':' in line and not line.strip().startswith('#'):
                field_name = line.split(':')[0].strip()
                field_value = line.split(
                    ':', 1)[1].strip() if ':' in line else ''

                # Add field to RST
                rst_content += f"``{field_name}``\n"
                rst_content += "~" * (len(field_name) + 4) + "\n\n"

                if field_description:
                    description_text = "\n".join(field_description)
                    # Replace multiple consecutive newlines with double newlines for RST paragraphs
                    description_text = re.sub(r'\n\s*\n', '\n\n', description_text)
                    # Format list items properly in the description text
                    # Insert a newline before list items for proper RST formatting
                    description_text = re.sub(r'([^\n])\n([ \t]*-[ \t]+)', r'\1\n\n\2', description_text)
                    
                    if "IMPORTANT:" in description_text:
                        description_text = description_text.replace(
                            "IMPORTANT:", "\n.. warning::\n").strip() + "\n\n"
                    elif "NOTE:" in description_text:
                        description_text = description_text.replace(
                            "NOTE:", "\n.. note::\n").strip() + "\n\n"
                    # Check if there's an example in the description
                    if "Example:" in description_text:
                        parts = description_text.split("Example:", 1)
                        rst_content += parts[0] + "\n\n"
                        example_text = parts[1].strip()
                        # Split the example into lines
                        example_lines = example_text.split("\n")
                        formatted_lines = []

                        # Process each line
                        for line in example_lines:
                            line = line.strip()
                            if line:
                                # Don't add tab to lines that already start with the field name in context
                                if field_name and line.startswith(field_name):
                                    formatted_lines.append(line)
                                else:
                                    formatted_lines.append("  " + line)

                        # Join the lines back with newlines
                        example_text = "\n".join(formatted_lines)
                        example_text
                        rst_content += "::\n\n\t" + example_text.replace("\n", "\n\t") + "\n\n"
                    else:
                        rst_content += description_text + "\n\n"
                
                if field_value:
                    rst_content += f"Default: ``{field_value}``\n\n"

            i += 1

    return rst_content


def main():
    # Define paths
    print("Generating input file documentation...")
    repo_root = Path('/home/rchaves/repo/workflow')
    yaml_template = repo_root / 'model_workflow' / \
        'resources' / 'inputs_file_template.yml'
    output_rst = repo_root / 'docs' / 'source' / 'input.rst'

    # Extract documentation
    rst_content = extract_yaml_documentation(yaml_template)

    # Write to RST file
    with open(output_rst, 'w') as f:
        f.write(rst_content)

    print(f"Documentation generated: {output_rst}")


if __name__ == "__main__":
    main()
