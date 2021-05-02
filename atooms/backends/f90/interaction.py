import numpy

from atooms.interaction import Interaction as _Interaction
from .helpers import _merge_source, _normalize_path
import f2py_jit

class Interaction(_Interaction):

    def __init__(self, model=None, potential='',
                 potential_parameters={}, cutoff='',
                 cutoff_parameters={}, interaction='interaction.f90',
                 helpers='helpers.f90', inline=True,
                 inline_safe=False, debug=False):
        """
        The interaction model can be defined in three ways:

        1) Passing a `model` string that matches any of the models
        defined in the atooms-models database (ex. "lennard_jones" or
        "gaussian_core")

        2) Passing a `model` dictionary with "potential" and "cutoff"
        keys and identical layout as the atooms-model database entries
        (ex. https://gitlab.info-ufr.univ-montp2.fr/atooms/models/blob/master/atooms/models/lennard_jones.json)

        3) Passing the `potential`, `cutoff`, `potential_parameters`
        and `cutoff_parameters` parameters. The `potential` and
        `cutoff` strings should match either

        a) potential and/or cutoff defined in atooms-models, such as
        "lennard_jones" or "cut_and_shift", or

        b) paths to Fortran 90 source codes that implement appropriate
        routines following the interfaces defined by the atooms-models
        package, see for instance
        https://gitlab.info-ufr.univ-montp2.fr/atooms/models/blob/master/atooms/models/lennard_jones.f90
        https://gitlab.info-ufr.univ-montp2.fr/atooms/models/blob/master/atooms/models/cut_shift.f90
        for an example of the interface the routines should implement.

        The parameters values are provided as dictionaries
        (`potential_parameters` and `cutoff_parameters``) matching the
        intent(in) variables entering the `setup()` routines of the
        fortran modules.
        """
        _Interaction.__init__(self, None)

        if model and not hasattr(model, 'get'):
            # This may be a string, so we look for the model in the
            # atooms-models database and replace the string with the dictionary
            try:
                from atooms.models import database
                model = database[model]
            except ImportError:
                raise ValueError('could not find model {}'.format(model))

        if model:
            # At this stage we expect a model dictionary
            assert len(model.get('potential')) == 1
            assert len(model.get('cutoff')) == 1
            potential = model.get('potential')[0].get('path')
            potential_parameters = model.get('potential')[0].get('parameters')
            cutoff = model.get('cutoff')[0].get('path')
            cutoff_parameters = model.get('cutoff')[0].get('parameters')

        self._module_path = None

        # Merge all sources into a unique source blob
        source = _merge_source(helpers, potential, cutoff, interaction)

        # Inline subroutines
        if inline:
            from f2py_jit.finline import inline_source
            # TODO: depending on f2py-jit version we can inline compute and smooth as well but this should be checked for bacward compatibility
            if inline_safe:
                source = inline_source(source, ignore='compute,smooth,tailor,forces')
            elif inline:
                source = inline_source(source, ignore='tailor,forces')

        # Compile and bundle the module with f2py
        if debug:
            extra_args = '--opt="-O3 -pg -fbounds-check -ffree-form -ffree-line-length-none"'
        else:
            extra_args = '--opt="-O3 -ffast-math -ffree-form -ffree-line-length-none"'

        # Build a unique module.
        # Every model with its own parameter combination corresponds to a unique module
        # and can be safely reused (up to changes in interaction / helpers)
        uid = f2py_jit.build_module(source,
                                    metadata={"interaction": interaction,
                                              "helpers": helpers,
                                              "potential": [potential, potential_parameters],
                                              "cutoff": [cutoff, cutoff_parameters]},
                                    extra_args=extra_args,
                                    db_file='.atooms_jit.json')

        # Setup potential and cutoff parameters
        _interaction = f2py_jit.import_module(uid)
        _interaction.potential.setup(**potential_parameters)
        _interaction.cutoff.setup(**cutoff_parameters)

        # Store module name (better not store the module itself, else we cannot deepcopy)
        self._uid = uid

        # Neighbors list
        self.neighbor_list = None

    def compute(self, observable, box, pos, ids):
        if self.neighbor_list is None:
            self._compute(observable, box, pos, ids)
        else:
            self._compute_with_neighbors(observable, box, pos, ids)

    def _compute(self, observable, box, pos, ids):
        _interaction = f2py_jit.import_module(self._uid)
        if observable in ['forces', 'energy']:
            if self.forces is None:
                self.forces = numpy.zeros_like(pos, order='F')
            self.energy, self.virial = _interaction.interaction.forces(box, pos, ids, self.forces)

        elif observable == 'gradw':
            if not hasattr(self, 'gradw'):
                self.gradw = numpy.zeros_like(pos, order='F')
            _interaction.interaction.gradw(box, pos, ids, self.gradw)

        elif observable == 'hessian':
            ndim, N = pos.shape
            if self.hessian is None:
                self.hessian = numpy.ndarray((ndim, N, ndim, N), order='F')
            _interaction.interaction.hessian(box, pos, ids, self.hessian)

    def _compute_with_neighbors(self, observable, box, pos, ids):
        f90 = f2py_jit.import_module(self._uid)
        if observable in ['forces', 'energy']:
            if self.forces is None:
                self.forces = numpy.zeros_like(pos, order='F')
            # TODO: rcut
            self.neighbor_list.adjust(box, pos, f90.cutoff.rcut_)
            self.neighbor_list.compute(box, pos, ids)
            self.energy, self.virial = f90.interaction_neighbors.forces(box, pos, ids,
                                                                        self.neighbor_list.neighbors,
                                                                        self.neighbor_list.number_of_neighbors,
                                                                        self.forces)

        elif observable == 'gradw':
            if not hasattr(self, 'gradw'):
                self.gradw = numpy.zeros_like(pos, order='F')
            _interaction.interaction.gradw(box, pos, ids, self.gradw)

        elif observable == 'hessian':
            ndim, N = pos.shape
            if self.hessian is None:
                self.hessian = numpy.ndarray((ndim, N, ndim, N), order='F')
            f90.interaction_neighbors.hessian(box, pos, ids,
                                              self.neighbor_list.neighbors,
                                              self.neighbor_list.number_of_neighbors,
                                              self.hessian)
