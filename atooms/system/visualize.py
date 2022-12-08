"""
Visualization functions for particles

Each show() method relies on a specific, optional visualization backend.
"""

# These functions used to belong to system.particle

_palette = [
    "#3c44aa",
    "#f0002f",
    "#f58f20",
    "#940650",
    "#f0e694",
]

def _hex_to_rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i + 2], 16) / 256 for i in (0, 2, 4))


def show_3dmol(particle, cell=None, radius=1.0, palette=None):
    """
    Visualize particles in cell using 3dmol http://3dmol.csb.pitt.edu/
    """
    from .particle import distinct_species
    import py3Dmol

    if palette is None:
        palette = ["#50514f", "#f25f5c", "#ffe066", "#247ba0", "#70c1b3",
                   "#0cce6b", "#c200fb", "#e2a0ff", "#6622cc", "#119822"]
    colors = {}
    for i, s in enumerate(distinct_species(particle)):
        colors[s] = palette[i]
    view = py3Dmol.view()
    view.setBackgroundColor('white')
    for p in particle:
        view.addSphere({'center': {'x': p.position[0], 'y': p.position[1], 'z': p.position[2]},
                        'radius': radius * p.radius, 'color': colors[p.species]})
    if cell is not None:
        view.addBox({'center': {'x': cell.center[0], 'y': cell.center[1], 'z': cell.center[2]},
                     'dimensions': {'w': cell.side[0], 'h': cell.side[1], 'd': cell.side[2]},
                     'wireframe': True, 'color': "#000000"})
    return view


def show_matplotlib(particle, cell, output_file=None, linewidth=3, alpha=0.3, show=False, now=False, outfile=None):
    """
    Make a snapshot of the `particle`s in the `cell` and save the
    image in `outfile`. The image is returned for further
    customization or visualization in jupyer notebooks.
    """
    import matplotlib.pyplot as plt
    from .particle import distinct_species
    if now:
        show = True
    if output_file is not None:
        outfile = output_file
    color_db = ['b', 'r', 'g', 'y']
    species = distinct_species(particle)
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.set_xlim((-cell.side[0]/2, cell.side[0]/2))
    ax.set_ylim((-cell.side[1]/2, cell.side[1]/2))
    for p in particle:
        c = plt.Circle(p.position[: 2], p.radius,
                       facecolor=color_db[species.index(p.species)],
                       edgecolor='black', alpha=alpha, linewidth=linewidth)
        ax.add_artist(c)
    if outfile is not None:
        fig.savefig(outfile, bbox_inches='tight')
    # TODO: show should be True by default
    if show:
        plt.show()
    return fig


# TODO: radius not working
def show_ovito(particle, cell, output_file=None, color='species',
               radius=0.35, viewport=None, callback=None, tmpdir=None,
               camera_dir=(0, 1, 0), camera_pos=(0, -10, 0),
               size=(640, 480), zoom=True, perspective=False,
               color_map='viridis', color_normalize=False, outfile=None):
    """
    Render image of particles in cell using ovito. The image is
    returned for visualization in jupyter notebooks.
    """
    import os
    from ovito.io import import_file  # pylint: disable=import-error
    from ovito.vis import Viewport, TachyonRenderer  # pylint: disable=no-name-in-module

    import tempfile
    from atooms.core.utils import mkdir

    if output_file is not None:
        outfile = output_file

    # Color coding
    color_attr = [getattr(p, color) for p in particle]
    is_discrete = isinstance(color_attr[0], (str, int))
    # Corresponding color system
    if is_discrete:
        # Discrete attribute
        color_db = sorted(list(set(color_attr)))
        color_db.sort()
        palette = [_hex_to_rgb(c) for c in _palette]
        colors = []
        for p in particle:
            colors.append(palette[color_db.index(getattr(p, color))])
    else:
        # Continuous attribute
        try:
            # Try with matplotlib
            import matplotlib
            import matplotlib.cm
            colormap = matplotlib.cm.get_cmap(color_map)
            if color_normalize:
                norm = matplotlib.colors.Normalize()
                norm.autoscale(color_attr)
                colors = colormap(norm(color_attr))
            else:
                colors = colormap(color_attr)
        except ImportError:
            # Fallback
            if color_normalize:
                c_min, c_max = min(color_attr), max(color_attr)
                if c_min != c_max:
                    colors = [[(c - c_min) / (c_max - c_min), 0.3,
                               (c_max - c) / (c_max - c_min)] for c in color_attr]
            else:
                colors = [[c, 0.3, c] for c in color_attr]

    # Make sure dirname exists
    if outfile is not None:
        mkdir(os.path.dirname(outfile))

    # Get a temporary file to write the sample
    fh = tempfile.NamedTemporaryFile('w', dir=tmpdir, suffix='.xyz', delete=False)
    tmp_file = fh.name

    # Self-contained EXYZ dump (it is not clean to use trajectories here)
    fh.write('{}\n'.format(len(particle)))
    fh.write('Properties=species:S:1:pos:R:3:radius:R:1:color:R:3 Lattice="{},0.,0.,0.,{},0.,0.,0.,{}"\n'.format(*cell.side))
    for i, p in enumerate(particle):
        c = colors[i]
        fh.write('{} {} {} {} {} {} {} {}\n'.format(p.species, *p.position, p.radius, *c))
    fh.close()

    # Ovito stuff. Can be customized by client code.
    pipeline = import_file(tmp_file)
    # Ovito seems to ignore the lattice info of exyz file
    # so we forcibly set the cell info here
    pipeline.source.data.cell_[0, 0] = cell.side[0]
    pipeline.source.data.cell_[1, 1] = cell.side[1]
    pipeline.source.data.cell_[2, 2] = cell.side[2]
    pipeline.source.data.cell_[:, 3] = -cell.side/2
    # Scale radius by default
    vis_element = pipeline.source.data.particles.vis
    vis_element.radius *= radius
    # Apply client code callback
    if callback:
        callback(pipeline)
    pipeline.add_to_scene()

    # Define viewport
    if viewport:
        vp = viewport
    else:
        if perspective:
            vp = Viewport(type=Viewport.Type.Perspective, camera_dir=camera_dir, camera_pos=camera_pos)
        else:
            vp = Viewport(type=Viewport.Type.Ortho, camera_dir=camera_dir, camera_pos=camera_pos)

    # Render
    if zoom:
        vp.zoom_all()
    if outfile is None:
        outfile = tmp_file + '.png'

    vp.render_image(filename=outfile,
                    size=size,
                    renderer=TachyonRenderer())

    # Scene is a singleton, so we must clear it
    pipeline.remove_from_scene()

    from atooms.core.utils import rmf
    rmf(tmp_file)

    # Try to display the image (e.g. in a jupyter notebook)
    try:
        from IPython.display import Image
        return Image(outfile)
    except ImportError:
        return outfile
