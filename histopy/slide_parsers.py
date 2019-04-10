"""
This module defines some useful classes and function to facilitate the parsing
of histopathology (Aperio?) slides.
"""
import os
import warnings
import math
import xml.etree.ElementTree as ET


import openslide
from openslide.deepzoom import DeepZoomGenerator
from shapely.geometry import Polygon, box
from shapely.ops import cascaded_union


__all__ = ['OpenSlide', 'DeepZoomTiler']


class OpenSlide(openslide.OpenSlide):
    """
    Wrapper around openslide.OpenSlide.

    If slide annotation is provided in the form of xml file, it used to
    populate the `self.annotation` list of annotations. Each annotation
    has two main properties:
    1. attrib: a dictionary of attributes
    2. polygon: a shapely.geometry.Polygon used to define the region of
                interest in terms of pixels.
    """

    def __init__(self, filename, xml=None):
        """Initialize OpenSlide object."""
        super().__init__(filename)
        self._filename = os.path.abspath(filename)
        self._format = super().detect_format(filename)
        mpp_x = self.properties[openslide.PROPERTY_NAME_MPP_X]
        mpp_y = self.properties[openslide.PROPERTY_NAME_MPP_Y]
        self._resolution = float(mpp_x), float(mpp_y)

        self.annotations = []
        if xml is not None:
            self.parse_annotations(xml)

    @property
    def filename(self):
        """Return path to slide."""
        return self._filename

    @property
    def format(self):
        """Return slide format."""
        return self._format

    @property
    def resolution(self):
        """Microns per Pixel"""
        return self._resolution

    def parse_annotations(self, xml):
        """
        Parse slide annotation from xml file.

        Annotations are appended to `self.annotations` as `ROIAnnotation` objects.

        Parameters
        ----------
        str
            path to an Aperio xml file with annotation data
        """
        tree = ET.parse(xml)
        root = tree.getroot()
        for annot in root.iterfind('./Annotation'):
            self.annotations.append(ROIAnnotation(annot))


class DeepZoomTiler(DeepZoomGenerator):
    """
    Wrapper aroudn openslide.deepzoom.DeepZoomGenerator.

    Major differences compared to DeepZoomGenerator:
        1. Has a default level attribute (`self.zoom_level`) for all
           DeepZoomGenerator "get" functions.
        2. Defines an iterator over tiles of any level (`self.itertiles`)
        3. Can compute overlap of tiles with slide annotations

    DeepZoom format details
    -----------------------
    DeepZoom is a xml-based file format that allows viewing of large images.
    A deepzoom file organizes the slide into a pyramid of increasing resolution slides.
    At the top of the pyramid (level=0) is a single tile depicting the whole slide.
    The whole slide at the "objective resolution" is at the bottom of the pyramid
    (level=self.level_count) partitioned into slides of size `tile_size` in pixels.
    Each level doubles the magnification of the previous level.
    See detailed description:
        - https://docs.microsoft.com/en-us/previous-versions/windows/silverlight/dotnet-windows-silverlight/cc645077(v=vs.95)
        - https://ysbecca.github.io/programming/2018/05/22/py-wsi.html
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
        annotations : dict of dict
            outer keys are slide coordinates (x, y, level)
            inner keys are slide annotation Names
            inner values are tile overlap with corresponding annotation.
            It's populated with calls to self.get_tile_annotation
        """
        if not isinstance(slide, openslide.OpenSlide):
            slide = OpenSlide(slide)
        super().__init__(slide, tile_size, overlap, limit_bounds)
        viewer_size = math.log2(tile_size + 2*overlap)
        if not viewer_size.is_integer():
            warnings.warn("For better viewer performance tile_size + 2*overlap"
                          "should be power of 2")
        # Memorize inputs
        self.slide = slide
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

        self.annotations = {}

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
        PIL.Image
            A pillow image of the tile

        """
        if level is None:
            level = self.zoom_level
        return super().get_tile(level, address)

    def __getitem__(self, address):
        """
        Return tile from address (x, y) in current zoom level.

        Convienient wrapper around `self.get_tile(address, level)`.

        Parameters
        ----------
        address : (int, int)
            tuple of tile location on the grid. Format is (col, row)

        Returns
        -------
        PIL.Image
            A pillow image of the tile
        """
        return super().get_tile(self.zoom_level, address)

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

    def get_tile_annotation(self, address, level=None):
        """
        Return tile overlap with all the annotations of the slide.

        The results are only computed the first time it's called on slide and
        are stored in `self.annotations`. Successive calls just retrieved memorized outcome.

        Parameters
        ----------
        address : (int, int)
            tuple of tile location on the grid. Format is (col, row)
        level : int
            level of Deep Zoom pyramid. If not specified, self.zoom_level is used.

        Returns
        -------
        dict:
            keys: annotation Name or Id if Name is empty
            value: overlaping area / total tile area
        """
        x, y = address
        level = self.zoom_level if level is None else level
        try:
            return self.annotations[x, y, level]
        except KeyError:
            pass
        (left, top), level, (width, height) = self.get_tile_coordinates(address, level)
        right, bottom = left + width, top - height
        tile = box(left, bottom, right, top)
        overlaps = {}
        for annot in self.slide.annotations:
            key = annot.attrib['Name'] if annot.attrib['Name'] else annot.attrib['Id']
            overlaps[key] = tile.intersection(annot.regions).area / tile.area
        self.annotations[x, y, level] = overlaps
        return overlaps

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
        # TODO: re-implement using async
        cols, rows = self.shape
        level = self.zoom_level if level is None else level
        for col in range(cols):
            for row in range(rows):
                address = (col, row)
                yield (address, super().get_tile(level, address))


class XmlNode:
    _float_keys = ()
    _bool_keys = ()
    _int_keys = ()

    def __init__(self, xmlnode):
        self.attrib = self.parse_attributes(xmlnode)

    def parse_attributes(self, node):
        """
        Convert region attributes to their natural type, based on educated guesses.
        """
        attrib = node.attrib
        for key in self._float_keys:
            try:
                attrib[key] = float(attrib[key])
            except KeyError:
                pass

        for key in self._int_keys:
            try:
                attrib[key] = int(attrib[key])
            except KeyError:
                pass

        for key in self._bool_keys:
            try:
                attrib[key] = attrib[key] == '1'
            except KeyError:
                pass

        return attrib


class ROIAnnotation(XmlNode):
    """
    A class to store information from '/Annotations/Annotation node of Aperio xml files.

    The class has 2 main attributes:
        - `attrib` which is a dictionary of the actual xml attributes
        - `regions` which is a `shapely` Polygon or MultiPolygon that defines
           the regions of interest
    """
    _bool_keys = ('ReadOnly', 'NameReadOnly', 'Incremental',
                  'LineColorReadOnly', 'Visible', 'Selected')

    def __init__(self, xmlnode):
        super().__init__(xmlnode)

        # Redefine Name if Attribute is uniquely specified
        attribute = xmlnode.findall('./Attributes/Attribute')
        if len(attribute) == 1 and attribute[0].attrib['Name']:
            if self.attrib['Name']:
                warnings.warn('In the given xml, both "Annotation" and '
                              '"Annotation/Attributes/Attribute" specify '
                              'a "Name" attribute. The latter will be used '
                              'but for no good reason')
            self.attrib['Name'] = attribute[0].attrib['Name']
        elif len(attribute) > 1:
            warnings.warn('Multilple Annotation/Attributes were unexpectedly found.  '
                          'All of them are ignored. Procceed at your own risk...')
        # TODO: better warning messages

        regions = [Region(n) for n in xmlnode.iterfind('./Regions/Region')]
        self.regions = self._combine_regions(regions)

    def __len__(self):
        """Number of Regions"""
        return len(self.regions)

    def __repr__(self):
        annotId = self.attrib['Id']
        label = self.attrib['Name'] if self.attrib['Name'] else None
        return "<ROIAnnotation {}:'{}' at {}>".format(annotId, label, id(self))

    def _combine_regions(self, regions):
        negative = [roi for roi in regions if roi.attrib['NegativeROA']]
        positive = [roi for roi in regions if not roi.attrib['NegativeROA']]

        for pos in positive:
            for neg in negative:
                # TODO: check negatives 'InputRegionId'?
                if pos.polygon.intersects(neg.polygon):
                    pos.polygon = pos.polygon - neg.polygon

        return cascaded_union([pos.polygon for pos in positive])


class Region(XmlNode):
    _float_keys = ('Zoom', 'ImageFocus', 'Length', 'Area', 'LengthMicrons', 'AreaMicrons')
    _bool_keys = ('Selected', 'NegativeROA', 'Analyze')

    def __init__(self, xmlnode):
        super().__init__(xmlnode)
        self.polygon = self.parse_vertices(xmlnode)

    def __len__(self):
        """Number of vertices"""
        return len(self.polygon.exterior.coords)

    def __repr__(self):
        rtype = 'Negative' if self.attrib['NegativeROA'] else 'Positive'
        return '<{} Region at {}>'.format(rtype, id(self))

    def parse_vertices(self, node):
        # TODO: check whether these are absolute coord or depend on Zoom/ImageFocus
        axes = ('X', 'Y', 'Z')
        vertices = [tuple(int(v.attrib[k]) for k in axes)
                    for v in node.iterfind('./Vertices/Vertex')]

        if len(set(z for x, y, z in vertices)) > 1:
            warnings.warn("3D slides are not supported yet. Z dimention will be ingored")

        polygon = Polygon([(x, y) for x, y, z in vertices])

        return polygon
