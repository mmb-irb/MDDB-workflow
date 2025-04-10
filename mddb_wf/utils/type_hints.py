# Import classes used for type hints
# DANI: Hay que importarlas siempre
# DANI: Sinó todos los scripts que usan e.g. Optional dan error porque no está importado
# DANI: Al cargar todos los módulos se producen errores de imports imposibles
# DANI: e.g. intentas importar structures, quien a su vez intenta importar los type hints

from pytraj import TrajectoryIterator
from typing import TYPE_CHECKING, Callable, List, Optional, Tuple, Union

if TYPE_CHECKING:
    from mddb_wf.utils.structures import Structure, Residue, Atom
    from mddb_wf.utils.register import Register
    from mddb_wf.utils.file import File
    from mddb_wf.utils.selections import Selection
    from mddb_wf.mwf import MD
