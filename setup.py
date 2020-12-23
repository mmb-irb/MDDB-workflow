from setuptools import setup, find_packages
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
    },)
