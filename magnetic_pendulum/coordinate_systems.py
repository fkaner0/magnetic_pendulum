from typing import Self  # NB: this only works for py > 3.11


class CoordinateSystem:
    def __init__(self, coords) -> None:
        self.coords = coords
        ## not finished
    
    def to(self, coord_system: Self) -> Self:
        """convert to another system - not sure if this is good struct"""
        if coord_system == type(self):
            return self
        else:
            new_system = coord_system.__name__
            function_name = f'_to_{new_system.lower()}'
            try:
                converter = getattr(self, function_name)
            except AttributeError:
                raise NotImplementedError
            return converter()
        

class Cartesian(CoordinateSystem):

    def _to_polar(self):
        pass

class Polar(CoordinateSystem):
    """NB: this polar form has theta defined as angle to NEGATIVE z axis"""

    def _to_cartesian(self):
        pass



