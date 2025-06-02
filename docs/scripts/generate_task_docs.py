#!/usr/bin/env python3

import re
import sys
import os
import inspect
from pathlib import Path
from model_workflow.mwf import requestables, DEPENDENCY_FLAGS, Project, MD

# Add the repository root to the Python path so we can import modules
repo_root = Path('..').resolve()
sys.path.insert(0, str(repo_root))


def clean_docstring(docstring):
    """Clean up a docstring for documentation."""
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
    
    return '\n'.join(lines)

def get_task_docstring(task_name, func):
    """Get a nice docstring for a task function."""
    docstring = clean_docstring(inspect.getdoc(func))
    return docstring

def generate_task_docs():
    """Generate documentation for all tasks in requestables."""
    
    # Separate project and MD tasks for better organization
    global project_tasks, md_tasks
    
    # Use inspection to determine which class a function belongs to
    project_tasks = []
    md_tasks = []
    
    for task_name, func in requestables.items():
        # Check if the function is a method
        if inspect.ismethod(func):
            if func.__self__ is Project:
                project_tasks.append(task_name)
            elif func.__self__ is MD:
                md_tasks.append(task_name)
        else:
            # For functions, inspect the first parameter to determine where it belongs
            try:
                signature = inspect.signature(func)
                params = list(signature.parameters.values())
                if params and params[0].name == 'self':
                    # Check the function's qualname to determine its class
                    if 'Project' in func.__qualname__:
                        project_tasks.append(task_name)
                    elif 'MD' in func.__qualname__:
                        md_tasks.append(task_name)
            except (ValueError, TypeError):
                # If we can't determine, put it in project tasks as a fallback
                project_tasks.append(task_name)
    
    # Start building the RST content
    rst_content = ".. _task_documentation: generated with generate_task_docs.py\n\n"
    rst_content += "Workflow Tasks\n"
    rst_content += "==================\n\n"
    rst_content += "This page documents the available tasks in the Model Workflow system.\n"
    rst_content += "These tasks can be specified with the ``-i`` (include) or ``-e`` (exclude) flags.\n\n"
    

    # Document project-level tasks
    rst_content += "Project Tasks\n"
    rst_content += "---------------\n\n"
    rst_content += "These tasks are executed once per project:\n\n"
    
    for task_name in sorted(project_tasks):
        func = requestables[task_name]
        rst_content += f"``{task_name}``\n"
        rst_content += "~" * (len(task_name) + 4) + "\n\n"
        
        docstring = get_task_docstring(task_name, func)
        rst_content += f"{docstring}\n\n"
    
    # Document MD-level tasks
    rst_content += "MD Tasks\n"
    rst_content += "-----------\n\n"
    rst_content += "These tasks are executed for each MD in the project:\n\n"
    
    # Categorize MD tasks
    md_file_tasks = []
    md_analysis_tasks = []
    
    for task_name in md_tasks:
        func = requestables[task_name]
        # Heuristic: file tasks are usually in input_files or processed_files
        if task_name.startswith('i') or task_name in ['structure', 'trajectory']:
            md_file_tasks.append(task_name)
        else:
            md_analysis_tasks.append(task_name)
    
    # First document input and processed files
    if md_file_tasks:
        rst_content += "Files\n"
        rst_content += "~~~~~~~~\n\n"
        
        for task_name in sorted(md_file_tasks):
            func = requestables[task_name]
            rst_content += f"``{task_name}``\n"
            rst_content += "^" * (len(task_name) + 4) + "\n\n"
            
            docstring = get_task_docstring(task_name, func)
            rst_content += f"{docstring}\n\n"
    
    # Then document analyses
    if md_analysis_tasks:
        rst_content += "Analyses\n"
        rst_content += "~~~~~~~~~~~~~~\n\n"
        
        for task_name in sorted(md_analysis_tasks):
            func = requestables[task_name]
            rst_content += f"``{task_name}``\n"
            rst_content += "^" * (len(task_name) + 4) + "\n\n"
            
            docstring = get_task_docstring(task_name, func)
            rst_content += f"{docstring}\n\n"
    
    # Document dependency flag groups
    rst_content += "Task Groups\n"
    rst_content += "-------------\n\n"
    rst_content += "These are predefined groups of tasks that can be specified with a single flag.\n\n"
    
    for flag, tasks in DEPENDENCY_FLAGS.items():
        rst_content += f"``{flag}``\n"
        rst_content += "~" * (len(flag) + 4) + "\n\n"
        rst_content += f"Includes the following tasks: {', '.join([f'``{t}``' for t in tasks])}\n\n"
    
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
