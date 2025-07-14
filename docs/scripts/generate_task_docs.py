#!/usr/bin/env python3

import re
import sys
import os
import inspect
from pathlib import Path
from model_workflow.mwf import project_requestables, md_requestables, DEPENDENCY_FLAGS, Task

# Add the repository root to the Python path so we can import modules
repo_root = Path('..').resolve()
sys.path.insert(0, str(repo_root))


def clean_docstring(docstring):
    """Clean up a docstring for documentation.

    Removes any "Args:", "Returns:", "Raises:" or "Notes:" sections from the docstring.
    """
    if not docstring:
        return ""
    
    # Remove leading/trailing whitespace and normalize internal whitespace
    lines = docstring.strip().split('\n')
    
    # Remove common leading whitespace from all lines after the first
    if len(lines) > 1:
        leading_spaces = min(len(line) - len(line.lstrip()) 
                            for line in lines[1:] if line.strip())
        for i in range(1, len(lines)):
            if lines[i]:
                lines[i] = lines[i][leading_spaces:] if leading_spaces <= len(lines[i]) else lines[i]
    
    # Remove Args, Returns, Raises sections
    new_lines = []
    skip_section = False
    for line in lines:
        if line.strip().startswith(('Args:', 'Returns:', 'Raises:', 'Notes:')):
            skip_section = True
            continue
        if not line.strip() and skip_section:
            skip_section = False
            continue
        if not skip_section:
            new_lines.append(line)

    s = ' '.join(new_lines)
    s = s[0].lower() + s[1:]
    return s


def get_github_link(func):
    """Get GitHub link for a function."""
    try:
        # Get the source file path
        source_file = inspect.getfile(func)
        # Get the line number
        lineno = inspect.getsourcelines(func)[1]
        # Convert to relative path from repo root
        rel_path = Path(source_file).relative_to(repo_root)
        return f"https://github.com/mmb-irb/MDDB-workflow/blob/master/{rel_path}#L{lineno}"
    except (OSError, ValueError):
        return None

def generate_task_docs():
    """Generate documentation for all tasks in requestables."""
    # Start building the RST content
    rst_content = ".. _task_documentation: generated with generate_task_docs.py\n\n"
    rst_content += "Workflow Tasks\n"
    rst_content += "==================\n\n"
    rst_content += "This page documents the available tasks in the MDDB Workflow system.\n"
    rst_content += "These tasks can be specified with the ``-i`` (include) or ``-e`` (exclude) flags.\n\n"
    

    # Document project-level tasks
    rst_content += "Project Tasks\n"
    rst_content += "---------------\n\n"
    rst_content += "These tasks are executed once per project:\n\n"
    for task_name in sorted(project_requestables):
        func = project_requestables[task_name]
        # Get the actual function if it's a Task object
        if isinstance(func, Task):
            func = func.func  
        rst_content += f"* ``{task_name}``"
        rst_content += f" `[source] <{get_github_link(func)}>`__"
        rst_content += f": {clean_docstring(func.__doc__)}"        
        rst_content += "\n\n"
    
    # Document MD-level tasks
    rst_content += "MD Tasks\n"
    rst_content += "-----------\n\n"
    rst_content += "These tasks are executed for each MD in the project:\n\n"
    
    # Categorize MD tasks
    md_file_tasks = []
    md_analysis_tasks = []
    
    for task_name in md_requestables:
        func = md_requestables[task_name]
        # Heuristic: file tasks are usually in input_files or processed_files
        if task_name in ['inpro', 'istructure', 'itrajectory', 'structure', 'trajectory']:
            md_file_tasks.append(task_name)
        else:
            md_analysis_tasks.append(task_name)
    
    # First document input and processed files
    if md_file_tasks:
        rst_content += "Files\n"
        rst_content += "~~~~~~~~\n\n"
        
        for task_name in sorted(md_file_tasks):
            func = md_requestables[task_name]
            if isinstance(func, Task):
                func = func.func  
            rst_content += f"* ``{task_name}``"
            rst_content += f" `[source] <{get_github_link(func)}>`__"
            rst_content += f": {clean_docstring(func.__doc__)}"
            rst_content += "\n\n"
    
    # Then document analyses
    if md_analysis_tasks:
        rst_content += "Analyses\n"
        rst_content += "~~~~~~~~~~~~~~\n\n"
        
        for task_name in sorted(md_analysis_tasks):
            func = md_requestables[task_name]
            if isinstance(func, Task):
                func = func.func    
            rst_content += f"* ``{task_name}``"
            rst_content += f" `[source] <{get_github_link(func)}>`__"
            rst_content += f": {clean_docstring(func.__doc__)}"
            rst_content += "\n\n"
    
    # Document dependency flag groups
    rst_content += "Task Groups\n"
    rst_content += "-------------\n\n"
    rst_content += "These are predefined groups of tasks that can be specified with a single flag.\n\n"
    
    for flag, tasks in DEPENDENCY_FLAGS.items():
        rst_content += f"* ``{flag}``: {', '.join([f'``{t}``' for t in tasks])}.\n\n"
    
    return rst_content

def main():
    print("Generating task documentation...")
    output_rst = repo_root / 'docs' / 'source' / 'tasks.rst'
    
    # Generate the documentation
    rst_content = generate_task_docs()
    
    # Write to RST file
    with open(output_rst, 'w') as f:
        f.write(rst_content)
    
    print(f"Task documentation generated: {output_rst}")

if __name__ == "__main__":
    main()
