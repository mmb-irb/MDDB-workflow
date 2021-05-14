from setuptools import setup, find_packages
from glob import glob

# Include all files in the resources directory
# Exclude the 'deprecated' directory
resources_directory = 'model_workflow/utils/resources'
deprecated_directory = resources_directory + '/deprecated'
resources_files = [ log for log in glob(resources_directory + '/*') if log != deprecated_directory ]

print(resources_files)

setup(
    name='model_workflow',
    version='0.0.1',
    python_requires='>=3.6.0',
    packages=find_packages(),
    install_requires=[
        'argparse'
    ],
    entry_points={
        'console_scripts': [
            'mwf = model_workflow.mwf:main',
        ]
    },
    data_files=[
        # Include all files in the resources directory but not recursively
        # This way, files in the 'dprecated' directory are ignored
        (resources_directory, resources_files)
    ]
)
