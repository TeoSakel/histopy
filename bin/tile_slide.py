"""
Aperio Slide Tiler.

Given a Aperio (svs) <slide> it partitions it into tiles and saves them in <outdir>.
A file coordinates.csv is also created in <outdir> with the coordinates of every tile
in the original slide.

Usage: tile_slide.py [options] <slide> <outdir>

Arguments:
    slide                           Slide to partition
    outdir                          Directory to write results

Options:
    -h --help                       Show this message and exit
    -o FORMAT --format=FORMAT       Image format to save tiles. Either jpeg or png [default: jpeg]
    -Q VALUE --quality=VALUE        For jpeg format, compression quality [default: 75]
    -s PIXEL --size=PIXEL           Tile size [default: 510]
    -e PIXEL --overlap=PIXEL        Number of extra pixels to add to each interior edge of a tile [default: 1]
    -L --limit-bounds               Renders only non-empty slide region
    -m LEVEL --magnification=LEVEL  Magnification level to use for tiling [default: 20.]
    -n --normalize                  Normalize staining (based on Macenko et al 2009)
    --heuristic=VALUE               Heuristics to estimate tissue coverage separated by comma.
                                    Possible values are: edge, dens.
    -b VALUE --beta=VALUE           Tissue coverage lower-cutoff separated by comma. Must be in [0,1)
"""
import os
import csv
import json
# from functools import partial  # to change parameters defaults

import numpy as np
from PIL import Image
from docopt import docopt
from histopy.slide_parsers import OpenSlide, DeepZoomTiler
from histopy.preprocessing import tissue_cover_edge, \
        tissue_cover_dens, tissue_cover_gray, normalize_staining


TILE_STATS = {
    'cover_edge': tissue_cover_edge,
    'cover_dens': tissue_cover_dens,
    'cover_gray': tissue_cover_gray
}


def _check_option(option, argv, type_str, type_fun):
    if option.startswith('-'):
        value = argv[option]
        msg = 'The option "{} {}" cannot be converted to {}'
    else:
        value = argv['--' + option]
        msg = 'The option "{}={}" cannot be converted to {}'
    msg = msg.format(option, value, type_str)
    try:
        ans = type_fun(value)
    except ValueError:
        raise ValueError(msg)
    return ans


def main(dz, outdir, criteria=tuple(), normalize=None, tile_fmt='jpeg', jpeg_qual=75):
    """Loop through and save tile images in outdir."""
    meta_file = os.path.join(outdir, "tile_metadata.csv")
    meta_fields = ['col', 'row',                         # tile coordinates in grid
                   'tile_width', 'tile_height',          # tile dimensions
                   'slide_level', 'slide_x', 'slide_y',  # tile coordinates in slide
                   'slide_width', 'slide_height']        # tile dimensions in slide
    meta_fields.extend(TILE_STATS.keys())

    if callable(normalize):
        norm = normalize
    else:
        norm = np.asarray  # does nothing

    with open(meta_file, 'w+', newline='') as f:
        # TODO: write json instead?
        meta = csv.writer(f)
        meta.writerow(meta_fields)
        for address, tile in dz.itertiles():
            col, row = address
            filename = "{x}_{y}.{ext}".format(x=col, y=row, ext=tile_fmt)
            filename = os.path.join(outdir, filename)
            # get tile meta
            px, py = tile.size  # size of tile
            (x, y), lvl, (w, h) = dz.get_tile_coordinates(address)
            tile_stats = {'col': col, 'row': row,
                          'tile_width': px, 'tile_height': py,
                          'slide_level': lvl, 'slide_x': x, 'slide_y': y,
                          'slide_width': w, 'slide_height': h}
            tile = np.asarray(tile)
            tile_stats.update({stat: estimator(tile)
                               for stat, estimator in TILE_STATS.items()})
            if any(tile_stats[crit] <= beta for crit, beta in criteria):
                continue
            tile = norm(tile)
            tile = Image.fromarray(tile)
            tile.save(filename, quality=jpeg_qual)
            meta.writerow(tile_stats[field] for field in meta_fields)


if __name__ == '__main__':
    argv = docopt(__doc__, version='0.1')
    # TODO: add `-f --force` option to skip overwritting files
    # TODO: add `-v --verbose` option to give some feedback
    # TODO: add more normalization methods (see preprocessing.py)

    outdir = argv['<outdir>']

    # Type-check tile options
    tile_fmt = argv['--format']
    assert tile_fmt in ('jpeg', 'jpg', 'png'), \
        'Image format must be jpeg or png. "{}" given'.format(tile_fmt)

    jpeg_qual = _check_option('quality',       argv, 'integer', int)
    tile_size = _check_option('size',          argv, 'integer', int)
    tile_over = _check_option('overlap',       argv, 'integer', int)
    tile_zoom = _check_option('magnification', argv, 'float',   float)
    tile_limit = argv['--limit-bounds']

    # Parse Normalization function
    if argv['--normalize']:
        normalize_fun = normalize_staining
    else:
        normalize_fun = None

    # Setup inclusion criteria
    # Parse Heuristic
    if argv['--heuristic'] is None:
        heuristics = ()
    else:
        heuristics = tuple('cover_' + h for h in argv['--heuristic'].split(','))
        for c in heuristics:
            if c not in TILE_STATS:
                raise ValueError('Uknown heuristic `{}` specified'.format(c))

    # Parse beta
    if argv['--beta'] is None:
        beta = ()
    else:
        beta = argv['--beta'].split(',')
        try:
            beta = list(map(float, beta))
        except ValueError as e:
            ValueError("Error while parsing beta parameter: " + e.args[0])
        if len(beta) == 1:
            beta *= len(heuristics)
        assert len(beta) == len(heuristics), \
            "Number of betas does not much number of heuristics"
        assert all(0 <= b < 1 for b in beta), \
            "Some betas are not in the [0, 1) interval"

    criteria = tuple(zip(heuristics, beta))  # combine

    # Read Silde and Tile
    slide = OpenSlide(argv['<slide>'])

    dz = DeepZoomTiler(slide, tile_size=tile_size, overlap=tile_over,
                       limit_bounds=tile_limit, magnification=tile_zoom)

    # write dz parameters to tiling_metadata.json
    os.makedirs(outdir, exist_ok=True)
    dz_meta = {'size':          dz.tile_size,
               'overlap':       dz.overlap,
               'bounds':        dz.limit_bounds,
               'grid':          dz.shape,
               'dimenions':     dz.dimensions,
               'magnification': dz.magnification,
               'normalize':     argv['--normalize'],
               'tile_format':   tile_fmt,
               'jpeg_qual':     jpeg_qual}
    dz_meta_path = os.path.join(outdir, 'tiling_params.json')
    with open(dz_meta_path, 'w') as fid:
        json.dump(dz_meta, fid)

    main(dz, outdir, criteria=criteria, normalize=normalize_fun, tile_fmt=tile_fmt, jpeg_qual=jpeg_qual)
