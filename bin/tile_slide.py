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
"""
import os
import csv
import json

from docopt import docopt
from histopy.slide_parsers import OpenSlide, DeepZoomTiler


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


def main(dz, outdir, tile_fmt='jpeg', jpg_qual=75):
    """Loop through and save tile images in outdir."""
    # TODO: add tile filters and tile normalization
    meta_file = os.path.join(outdir, "tile_metadata.csv")
    meta_fields = ('col', 'row', 'tile_width', 'tile_height',
                   'slide_level', 'slide_x', 'slide_y',
                   'slide_width', 'slide_height')
    with open(meta_file, 'w+', newline='') as f:
        # TODO: write json instead?
        meta = csv.writer(f)
        meta.writerow(meta_fields)
        for address, tile in dz.itertiles():
            col, row = address
            filename = "{x}_{y}.{ext}".format(x=col, y=row, ext=tile_fmt)
            filename = os.path.join(outdir, filename)
            px, py = tile.size  # size of tile
            (x, y), lvl, (w, h) = dz.get_tile_coordinates(address)
            tile_stats = {'col': col, 'row': row,
                          'tile_width': px, 'tile_height': py,
                          'slide_level': lvl, 'slide_x': x, 'slide_y': y,
                          'slide_width': w, 'slide_height': h}
            tile.save(filename, quality=jpg_qual)
            meta.writerow(tile_stats[field] for field in meta_fields)


if __name__ == '__main__':
    argv = docopt(__doc__, version='0.1')
    # TODO: add (-f --force) option to skip overwritting files
    # TODO: add (-v --verbose) option to give some feedback

    # Type-check options
    tile_fmt = argv['--format']
    assert tile_fmt in ('jpeg', 'jpg', 'png'), \
        'Image format must be jpeg or png. "{}" given'.format(tile_fmt)

    jpg_qual   = _check_option('quality',       argv, 'integer', int)
    tile_size  = _check_option('size',          argv, 'integer', int)
    tile_over  = _check_option('overlap',       argv, 'integer', int)
    tile_zoom  = _check_option('magnification', argv, 'float',   float)
    tile_limit = argv['--limit-bounds']

    slide = OpenSlide(argv['<slide>'])
    outdir = argv['<outdir>']

    dz = DeepZoomTiler(slide,
                       tile_size=tile_size,
                       overlap=tile_over,
                       limit_bounds=tile_limit,
                       magnification=tile_zoom)

    # write dz parameters to tiling_metadata.json
    os.makedirs(outdir, exist_ok=True)
    dz_meta = {'size': dz.tile_size,
               'overlap': dz.overlap,
               'bounds': dz.limit_bounds,
               'grid': dz.shape,
               'dimenions': dz.dimensions,
               'magnification': dz.magnification,
               'tile_format': tile_fmt,
               'jpeg_qual': jpg_qual}
    dz_meta_path = os.path.join(outdir, 'tiling_params.json')
    with open(dz_meta_path, 'w+') as fid:
        json.dump(dz_meta, fid)

    main(dz, outdir, tile_fmt, jpg_qual)
