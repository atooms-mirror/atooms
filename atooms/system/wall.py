import numpy


class Wall:

    """
    Examples:
    ---------
    >>> w = Wall([1.0, 0.0, 0.0], [0,0,0])
    >>> w.distance(numpy.array([1.0, 0.0, 0.0]))
    >>> w.distance(numpy.array([[-2.5, 2.0, 1.0], [1.0, 0.0, 0.0]]).transpose())
    """

    def __init__(self, vector, point, species='A'):
        """A plane specified by `point` and a `vector` orthogonal to the plane"""
        vec = numpy.array(vector)
        self.normal = vec / numpy.sum(vec**2)**0.5
        self.point = numpy.array(point)
        self.species = species
        self._offset = numpy.dot(self.point, self.normal)

    def distance(self, point):
        """
        Return the vector distance between the `point` and the wall, where
        `point` is an array whose first dimension equals the spatial
        dimensionality.
        """
        scalar_distance = (numpy.dot(self.normal, point) - self._offset)
        if len(point.shape) == 1:
            return scalar_distance * self.normal
        elif len(point.shape) == 2:
            return scalar_distance * self.normal[:, None]
        else:
            raise ValueError('input vector has more than 2 dimensions')
