# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""
Base and total interaction class.

Actual interaction backends should implement this interface or
subclass Interaction by implementing the compute() method and
specifying the variables passed to it using the Interaction.variables
dictionary.
"""

import numpy


class InteractionBase(object):

    def __init__(self):
        self.variables = {'position': 'particle.position'}
        """
        A dictionary of variables needed to compute the interaction

        The keys must match the variable names in the
        `Interaction.compute()` interface, the fields must be
        canonicalizable variables accepted by `System.dump()`.

        It is possible to specify the required data type using the
        optionl colon syntax <property>[:<dtype>]. The dtype must a
        valid identifier for numpy array creation.
        """
        # TODO: order is not a good variable, we should expand the syntax using [:order]
        self.order = 'F'  # deprecated
        self.observable = ['energy', 'forces', 'virial', 'stress', 'hessian']
        for observable in self.observable:
            setattr(self, observable, None)

    def __add__(self, other):
        total = Interaction()
        for attr in self.observable:
            # Sanity check
            err = getattr(self, attr) is not None and getattr(other, attr) is None or \
                  getattr(self, attr) is None and getattr(other, attr) is not None
            assert not err, 'attribute {} not set in {} or {}'.format(attr, self, other)

            # Store the sum of the properties in the total interaction
            if getattr(self, attr) is not None and getattr(other, attr) is not None:
                setattr(total, attr, getattr(self, attr) + getattr(other, attr))
        return total

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def compute(self, observable, position):
        """
        Compute an `observable` from this interaction

        Subclasses must adapt the interface so as to match the keys
        specified by the `Interaction.variables` dict.
        """
        # Sanity checks
        assert observable in self.observable, \
            'unsupported observable %s'.format(observable)

        # Zeroing observables
        ndim, N = position.shape
        self.energy = None
        self.virial = None
        self.stress = None
        self.forces = None
        if observable == 'energy':
            self.energy = 0.0
        elif observable == 'forces' or observable is None:
            self.energy = 0.0
            self.virial = 0.0
            self.forces = numpy.zeros_like(position)
        elif observable == 'stress':
            self.energy = 0.0
            self.virial = 0.0
            self.forces = numpy.zeros_like(position)
            self.stress = numpy.zeros(ndim, ndim)
        elif observable == 'hessian':
            self.hessian = numpy.zeros((ndim, N, ndim, N))


class Interaction(InteractionBase):

    def __init__(self, *terms):
        InteractionBase.__init__(self)
        self.variables = {}
        self.term = []
        for term in terms:
            self.add(term)

    def add(self, term):
        assert set(self.observable) == set(term.observable), 'observables differ'
        self.term.append(term)
        self.variables.update(term.variables)

    def compute(self, observable, **kwargs):
        """
        Compute an `observable` from all terms of this interaction
        """
        if len(self.term) == 0:
            return

        for term in self.term:
            # Extract the relevant variables for this term
            term_kwargs = {}
            for key, value in term.variables.items():
                term_kwargs[key] = kwargs[key]
            term.compute(observable, **term_kwargs)

        # Sum all interaction terms
        total = sum(self.term)
        for attr in self.observable:
            setattr(self, attr, getattr(total, attr))
