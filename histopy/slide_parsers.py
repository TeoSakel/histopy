"""
This module defines some useful classes and function to facilitate the parsing
of histopathology (Aperio?) slides.
"""
import os
import warnings
import math

import openslide
from openslide.deepzoom import DeepZoomGenerator


class OpenSlide(openslide.OpenSlide):
    """
    Simple wrapper around openslide.OpenSlide.

    It just remembers its input and stores the slide format.
    (useful for debugging)
    """

    # TODO: add information about Region of Interest (ROI)
    # see scratch/aperio_xml.py

    def __init__(self, filename, xml=None):
        """Initialize OpenSlide object."""
        super().__init__(filename)
        self._filename = os.path.abspath(filename)
        self._format = super().detect_format(filename)

    @property
    def filename(self):
        """Return path to slide."""
        return self._filename

    @property
    def format(self):
        """Return slide format."""
        return self._format


class DeepZoomTiler(DeepZoomGenerator):
    """
    Wrapper aroudn openslide.deepzoom.DeepZoomGenerator.

    Allows the user to select a tile level/magnification and iterate over the tiles.

    DeepZoom format details
    -----------------------
    DeepZoom is a xml-based file format that allows viewing of large images.
    A deepzoom file organizes the slide into a pyramid of increasing resolution slides.
    At the top of the pyramid (level=0) is a single tile depicting the whole slide.
    The whole slide at the "objective resolution" is at the bottom of the pyramid
    (level=self.level_count) partitioned into slides of size `tile_size` in pixels.
    Each level doubles the magnification of the previous level.
    See detailed descriptio here:
    https://docs.microsoft.com/en-us/previous-versions/windows/silverlight/dotnet-windows-silverlight/cc645077(v=vs.95)
    """

    def __init__(self, slide, tile_size=254, overlap=1, limit_bounds=False,
                 zoom_level=None, magnification=None):
        """
        Initialize DeepZoomTiler based on user specification.

        Parameters
        ----------
        slide : OpenSlide
            the slide used to generate the DeepZoom object.
        tile_size : int
            the size (in pixels) of the tiles at the bottom of the DeepZoom pyramid.
        overlap : int
            the number of extra pixels to add to each interior edge of a tile.
        limit_bounds : bool
            whether to render only the non-empty slide region.
        zoom_level : int
            the selected magnification level
        magnification : float
            the selected magnification (overwrites zoom_level)

        Attributes
        ----------
        All the input parameters are stored as attributes. In addition to:
        level_magnification : tuple of floats
            all the available magnification levels to choose from
        shape : (int, int)
            the grid shape for the current zoom level
        dimensions : (int, int)
            the slide dimentions (in pixels) for the current zoom level

        Yields
        ------
        PIL.Image
            tiles per row and column

        """
        if not isinstance(slide, openslide.OpenSlide):
            slide = OpenSlide(slide)
        super().__init__(slide, tile_size, overlap, limit_bounds)
        viewer_size = math.log2(tile_size + 2*overlap)
        if not viewer_size.is_integer():
            warnings.warn("For better viewer performance tile_size + 2*overlap"
                          "should be power of 2")
        # Memorize inputs
        self.slide = slide  # TODO: is this necessary?
        self._tile_size = tile_size
        self._overlap = overlap
        self._limit_bounts = limit_bounds
        self._zoom_level = self.level_count - 1  # default_choice [0-index]
        if zoom_level is not None:
            self.zoom_level = zoom_level
        # Extra functionality
        # original (objective) magnification
        try:
            self.obj_power = float(slide.properties[openslide.PROPERTY_NAME_OBJECTIVE_POWER])
        except KeyError:
            self.obj_power = 1.0
        level_magnification = [self.obj_power / pow(2, n) for n in range(self.level_count)]
        level_magnification.reverse()  # high = high resolution
        self._level_magnification = tuple(level_magnification)
        if magnification is not None:
            self.magnification = magnification

    @property
    def tile_size(self):
        """Return the tile size specified by the user."""
        return self._tile_size

    @property
    def overlap(self):
        """Return overlap value specified by the user."""
        return self._overlap

    @property
    def limit_bounds(self):
        """Return boolean whether the non-empty slide region is rendered."""
        return self._limit_bounds

    @property
    def zoom_level(self):
        """Return selected zoom level (counting from 0)."""
        return self._zoom_level

    @zoom_level.setter
    def zoom_level(self, value):
        assert 0 <= value <= self.level_count - 1, \
                "Value outside range of levels [0, {}]".format(self.level_count)
        self._zoom_level = value

    @property
    def level_magnification(self):
        """List all the magnification levels available."""
        return self._level_magnification

    @property
    def magnification(self):
        """Return selected magnification level."""
        return self.level_magnification[self.zoom_level]

    @magnification.setter
    def magnification(self, value):
        """Select a magnification level closest to "value"."""
        # each level doubles the magnification until the objective magnification
        # so to reach the required "valued" we have to compure the number of
        # halvings (/2) from the objective (obj_power / value)
        offset = round(self.obj_power / value / 2)
        self.zoom_level = self.level_count - offset - 1  # 0-indexed

    @property
    def shape(self):
        """Return tile-grid dimensions."""
        return self.level_tiles[self.zoom_level]

    @property
    def dimensions(self):
        """Return slide dimensions at current zoom level (in pixels)."""
        return self.level_dimensions[self.zoom_level]

    def __len__(self):
        """Return number of tiles in the current level."""
        cols, rows = self.shape
        return rows * cols

    def get_tile(self, address, level=None):
        """
        Return a PIL.Image of the tile at the specified address and level of the pyramid.

        This function is a wrapper around:
            `openslide.deepzoom.DeepZoomGenerator.get_tile`

        Parameters
        ----------
        address : (int, int)
            tuple of tile location on the grid. Format is (col, row)
        level : int
            level of Deep Zoom pyramid. If not specified, self.zoom_level is used.

        Returns
        -------
        PIL.Image.Image
            A pillow image of the tile

        """
        # TODO: rewrite as __getitem__ ?
        if level is None:
            level = self.zoom_level
        return super().get_tile(level, address)

    def get_tile_coordinates(self, address, level=None):
        """
        Return coordinates of the speciefied tile in original slide.

        The result can be used with `OpenSlide.read_region`
        This function is a wrapper around:
            `openslide.deepzoom.DeepZoomGenerator.get_tile_coordinates`

        Parameters
        ----------
        address : (int, int)
            tuple of tile location on the grid. Format is (col, row)
        level : int
            level of Deep Zoom pyramid. If not specified, self.zoom_level is used.

        Returns
        -------
        location : (int, int)
            tuple giving the top left pixel in the level 0 reference frame (original slide).
        level : int
            level number of the reference frame
        size : (int, int)
            tuple giving the region size

        """
        if level is None:
            level = self.zoom_level
        return super().get_tile_coordinates(level, address)

    def get_tile_dimensions(self, address, level=None):
        """
        Return the size of the specified tuple in pixels.

        This function is a wrapper around:
            `openslide.deepzoom.DeepZoomGenerator.get_tile_dimensions`

        Parameters
        ----------
        address : (int, int)
            tuple of tile location on the grid. Format is (col, row)
        level : int
            level of Deep Zoom pyramid. If not specified, self.zoom_level is used.

        Returns
        -------
        (int, int)
            a (pixel_x, pixel_y) tuple

        """
        if level is None:
            level = self.zoom_level
        return super().get_tile_dimensions(level, address)

    def itertiles(self, level=None):
        """
        Iterate over the tiles of specified level (zoom_level if None).

        Parameters
        ----------
        level : int
            zoom level to use. If None self.zoom_level is used

        Yields
        ------
        ((int, int) PIL.Image)
            a tuple of the tile address (col, row) and the tile as PIL.Image

        """
        cols, rows = self.shape
        level = self.zoom_level if level is None else level
        for col in range(cols):
            for row in range(rows):
                address = (col, row)
                yield (address, super().get_tile(level, address))

    # TODO: if slide has ROIs, compute overlap for tile and store in
    # a dictionary of the form (lvl, col, row): (overlap_0, overlap_1, ...)

