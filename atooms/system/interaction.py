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
            if getattr(self, attr) is None and getattr(other, attr) is None:
                continue
            elif getattr(self, attr) is not None and getattr(other, attr) is not None:
                # Store the sum of the properties in the total interaction
                setattr(total, attr, getattr(self, attr) + getattr(other, attr))
            else:
                raise ValueError('attribute {} not set in {} or {}'.format(attr,
                                                                           self,
                                                                           other))
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
            'unsupported observable {}'.format(observable)

        # Zeroing observables
        ndim, N = position.shape
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
            for key in term.variables:
                term_kwargs[key] = kwargs[key]
            term.compute(observable, **term_kwargs)

        # Sum all interaction terms
        total = sum(self.term)
        for attr in self.observable:
            setattr(self, attr, getattr(total, attr))
