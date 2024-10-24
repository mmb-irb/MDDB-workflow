from typing import TYPE_CHECKING 

if TYPE_CHECKING: # 
    from pytraj import TrajectoryIterator
    from model_workflow.utils.structures import Structure
    from model_workflow.utils.register import Register
    from model_workflow.utils.file import File