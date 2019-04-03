"""
Rip-off from: https://github.com/CODAIT/deep-histopath/blob/master/deephistopath/preprocessing.py

Original code distributed under ASF v2.0
"""
from functools import partial
import numpy as np

from scipy.ndimage.morphology import binary_fill_holes
from skimage.color import rgb2gray
from skimage.feature import canny
from skimage.morphology import binary_closing, binary_dilation, disk

# TODO: implement Vahadane et al: https://ieeexplore.ieee.org/document/7460968
# see: https://github.com/Peter554/StainTools/blob/master/staintools/stain_extraction/vahadane_stain_extractor.py
# and: https://github.com/JisuiWang/Stain_Normalization


def tissue_cover_gray(img, thres=0.5):
    """
    Estimate tissue coverage based on gray coverage.

    Naive filter to estimate tissue coverage of image. Tissue is assumed to be
    darker than background. The image is converted to grayscale and
    binarized based on a user-defined threshold and the mean is return.

    Parameters
    ----------
    img : RGB image as numpy array
    thres : float
        threshold to binarize image

    Returns
    -------
    float : the percentage of image covered by tissue

    """
    img = np.asarray(img)  # skimage.color.rgb2gray != PIL.Image.convert('L')
    img = rgb2gray(img)
    return (img < thres).mean()


def tissue_cover_edge(img, detector=partial(canny, sigma=2.),
                      closing_kernel=disk(10), dilation_kernel=disk(10)):
    """
    Calculate tissue coverage based on an edge-heuristic.

    The edges of the tissue are identified using the detector argument.
    Then a hole-filling heuristic is applied to cover holes in the tissue.
    This in effect computes a "tissue convex-hull".
    The output is the percent of the image covered by this hull.

    Parameters
    ----------
    img : RGB image
    detector : callable
        a function that takes as input the image and
        retuns the binary image of edges
    closing_kernel : an image kernel
        used to apply closing operation (see fill_minor_holes)
    dilation_kernel : an image kernel
        used to apply dilation operation (see fill_minor_holes)

    Returns
    -------
    float : the percentage of image covered by tissue

    """
    img = np.asarray(img)
    img = 1 - rgb2gray(img)  # 1: dense tissue, 0: background
    img = detector(img)      # detect edges of tissue
    img = fill_minor_holes(img, closing_kernel, dilation_kernel)
    return img.mean()


def tissue_cover_dens(img, beta=0.15, closing_kernel=disk(10), dilation_kernel=disk(10)):
    """
    Calculate tissue coverage based on an optical density heuristic.

    The optical density of each channel is computed as:
        od = -np.log10(x/255 + 1e-8)
    Then a hole-filling heuristic is applied to cover holes in the tissue.
    This in effect computes a "tissue convex-hull".
    The output is the percent of the image covered by this hull.

    Parameters
    ----------
    img : RGB image
    beta : float
        used to threshold the optical density.
        Values greater than beta are considered tissue.
    closing_kernel : an image kernel
        used to apply closing operation (see fill_minor_holes)
    dilation_kernel : an image kernel
        used to apply dilation operation (see fill_minor_holes)

    Returns
    -------
    float : the percentage of image covered by tissue

    """
    img = np.asarray(img)
    img = -np.log10(img/255 + 1e-8)
    img = img.min(axis=2) >= beta
    img = fill_minor_holes(img, closing_kernel, dilation_kernel)
    return img.mean()


def fill_minor_holes(img, closing_kernel=disk(10), dilation_kernel=disk(10)):
    """
    Run simple heuristic to fill small holes in tissue.

    It's 3 morphological operation applied in order:
        1. binary_closing: used to connect tissue edges
        2. binary dilation: makes tissue edges bolder (holes smaller)
        3. binary_fill_holes: fill in small holes
    In effect this algorithm identifies a "tissue convex hull".

    Parameters
    ----------
    img : an binary image as a numpy array
    closing_kernel : image kernel used for closing
    dilation_kernel : image kernel used for dilation

    Returns
    -------
    numpy.ndarray :
        transformed binary image of the same size

    """
    img = binary_closing(img, closing_kernel)    # connects edges
    img = binary_dilation(img, dilation_kernel)  # make tissue bolder
    img = binary_fill_holes(img)                 # fills in holes
    return img


# https://ieeexplore.ieee.org/abstract/document/5193250
def normalize_staining(img, beta=0.15, alpha=1, stain_ref=None, max_sat_ref=None):
    """
    Normalize the staining of H&E histology slides.

    Parameters
    ----------
    img : an RGB image
    beta : float
        threshold to remove pixels if channels have density/absorbanse less than it
    alpha : float between (0, 100)
        percentile used for define robust extermes of angles
    stain_ref : numpy array of shape (3, 2)
        reference values for each channel (3) for the 2 stains
    max_sat_ref : numpy aarray of length 2
        reference maximum saturation values for each stain

    Returns
    -------
    A (slide_num, sample) tuple, where the sample is a 3D NumPy array
    of shape (H,W,C) that has been stain normalized.

    References
    ----------
    - Macenko, Marc, et al. "A method for normalizing histology slides
    for quantitative analysis." Biomedical Imaging: From Nano to Macro,
    2009.  ISBI'09. IEEE International Symposium on. IEEE, 2009.
      - http://wwwx.cs.unc.edu/~mn/sites/default/files/macenko2009.pdf
    - https://github.com/mitkovetta/staining-normalization

    """
    assert 0 < alpha < 100, "alpha must be in the (0, 100) interval"
    # Constants
    Io = 255.  # normalizing constant
    # Reference stain vectors and stain saturations.  We will normalize all slides
    # to these references.  To create these, grab the stain vectors and stain
    # saturations from a desirable slide.
    if stain_ref is None:
        stain_ref = np.array([[0.54598845, 0.32211600],
                              [0.72385198, 0.76419107],
                              [0.42182333, 0.55879629]])
    else:
        assert stain_ref.shape == (3, 2)

    if max_sat_ref is None:
        max_sat_ref = np.array([0.82791151, 0.61137274]).reshape(2, 1)
    else:
        assert max_sat_ref.shape == (2, 1)

    fimg = np.asarray(img)
    h, w, c = fimg.shape
    fimg = fimg.reshape(-1, c).astype(float)  # shape (K=H*W, C)

    # Convert RGB to Optical Density.
    OD = -np.log10(fimg/Io + 1e-8)

    # Remove data with OD intensity less than beta (transparent pixels)
    OD_thresh = OD[np.all(OD >= beta, 1), :]  # shape (K, C)
    if OD_thresh.shape[0] < 3:
        return np.zeros_like(img) + 255  # return all white image

    # Calculate eigenvectors.
    # Note: We can either use eigenvector decomposition, or SVD.
    _, _, V = np.linalg.svd(OD_thresh, full_matrices=False)

    # Extract two largest eigenvectors.
    # Note: We swap the sign of the eigvecs here to be consistent
    # with other implementations.  Both +/- eigvecs are valid, with
    # the same eigenvalue, so this is okay.
    top_eigvecs = V[:2, :].T * -1  # shape (C, 2)

    # Project thresholded optical density values onto plane spanned by
    # 2 largest eigenvectors.
    proj = np.dot(OD_thresh, top_eigvecs)  # shape (K, 2)

    # Calculate angle of each point wrt the first plane direction.
    # Note: the parameters are `np.arctan2(y, x)`
    angles = np.arctan2(proj[:, 1], proj[:, 0])  # shape (K,)

    # Find robust extremes (a and 100-a percentiles) of the angle.
    min_angle = np.percentile(angles, alpha)
    max_angle = np.percentile(angles, 100-alpha)

    # Convert min/max vectors (extremes) back to optimal stains in OD space.
    # This computes a set of axes for each angle onto which we can project
    # the top eigenvectors.  This assumes that the projected values have
    # been normalized to unit length.
    extreme_angles = np.array([[np.cos(min_angle), np.cos(max_angle)],
                               [np.sin(min_angle), np.sin(max_angle)]])
    stains = np.dot(top_eigvecs, extreme_angles)  # shape (C, 2)

    # Merge vectors with hematoxylin first, and eosin second, as a heuristic.
    if stains[0, 0] < stains[0, 1]:
        stains[:, [0, 1]] = stains[:, [1, 0]]  # swap columns

    # Calculate saturations of each stain.
    # Note: Here, we solve
    #    OD = VS
    #     S = V^{-1}OD
    # where `OD` is the matrix of optical density values of our image,
    # `V` is the matrix of stain vectors, and `S` is the matrix of stain
    # saturations.  Since this is an overdetermined system, we use the
    # least squares solver, rather than a direct solve.
    sats, _, _, _ = np.linalg.lstsq(stains, OD.T, rcond=None)

    # Normalize stain saturations to have same pseudo-maximum based on
    # a reference max saturation.
    max_sat = np.percentile(sats, 99, axis=1, keepdims=True)
    sats = sats / max_sat * max_sat_ref

    # Compute optimal OD values.
    OD_norm = np.dot(stain_ref, sats)

    # Recreate image.
    img_norm = 10**(-OD_norm) * Io - 1e-8
    img_norm = img_norm.T.reshape(h, w, c)
    img_norm = np.clip(np.round(img_norm), 0, 255).astype(np.uint8)

    return img_norm
