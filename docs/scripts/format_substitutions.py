from mddb_workflow.tools.check_inputs import (
    TRAJECTORY_SUPPORTED_FORMATS,
    TOPOLOGY_SUPPORTED_FORMATS,
    STRUCTURE_SUPPORTED_FORMATS
)


def setup(app):
    """Sphinx extension to add supported formats as substitutions."""
    # Add as RST prolog to make them available in all documents
    app.config.rst_prolog = f"""
.. |trajectory_formats| replace:: {', '.join(TRAJECTORY_SUPPORTED_FORMATS)}
.. |topology_formats| replace:: {', '.join(TOPOLOGY_SUPPORTED_FORMATS)}
.. |structure_formats| replace:: {', '.join(STRUCTURE_SUPPORTED_FORMATS)}
"""
    return {'version': '0.1'}
