# Import classes used for type hints
# DANI: Hay que importarlas siempre
# DANI: Sinó todos los scripts que usan e.g. Optional dan error porque no está importado
# DANI: Al cargar todos los módulos se producen errores de imports imposibles
# DANI: e.g. intentas importar structures, quien a su vez intenta importar los type hints

from pytraj import TrajectoryIterator
from typing import TYPE_CHECKING, Callable, Optional, Union, Generator
Coords = tuple[float, float, float]

if TYPE_CHECKING:
    from mddb_workflow.utils.structures import Structure, Residue, Atom
    from mddb_workflow.utils.register import Register
    from mddb_workflow.utils.cache import Cache
    from mddb_workflow.utils.file import File
    from mddb_workflow.utils.selections import Selection
    from mddb_workflow.mwf import Task, MD, Project
    from MDAnalysis import Universe
