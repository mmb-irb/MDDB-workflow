# Import classes used for type hints
# DANI: Hay que importarlas siempre
# DANI: Sinó todos los scripts que usan e.g. Optional dan error porque no está importado
# DANI: Al cargar todos los módulos se producen errores de imports imposibles
# DANI: e.g. intentas importar structures, quien a su vez intenta importar los type hints

from pytraj import TrajectoryIterator
from typing import TYPE_CHECKING, Callable, List, Set, Optional, Tuple, Union, Generator
Coords = Tuple[float, float, float]

if TYPE_CHECKING:
    from model_workflow.utils.structures import Structure, Residue, Atom
    from model_workflow.utils.register import Register
    from model_workflow.utils.cache import Cache
    from model_workflow.utils.file import File
    from model_workflow.utils.selections import Selection
    from model_workflow.mwf import Task, MD, Project
    from MDAnalysis import Universe
