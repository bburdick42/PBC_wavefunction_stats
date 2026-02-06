import numpy as np
from numpy import linalg as LA
from scipy.ndimage import center_of_mass, variance

#########################################################################################################

def length_to_index(vector, array, lengths):
    """Turn a vector of lengths into a vector of indices in the array"""
    dims = array.shape
    n = np.round(dims * vector / lengths)
    return n.astype(int)

def index_to_length(index, array, lengths):
    """Turn a vector of indices into a vector of lengths corresponding to a position in the array"""
    dims = array.shape
    d = np.array([lengths[i]/dims[i] for i in range(3)])
    return d * index

def calc_periodic_mean(array, lengths):
    """Calculate the pseudo center of mass according to this paper: https://pubs.aip.org/aip/jcp/article/162/20/204103/3347482/On-the-estimation-of-center-of-mass-in-periodic"""
    dims = array.shape
    d = np.array([lengths[i]/dims[i] for i in range(3)])
    grids = np.meshgrid(*[np.arange(d) for d in dims], indexing="ij")

    abs_array = np.abs(array)**2

    xi_bar   = np.array([np.sum(np.cos(2*np.pi * grids[i] * d[i] / lengths[i]) * abs_array) for i in range(3)])
    zeta_bar = np.array([np.sum(np.sin(2*np.pi * grids[i] * d[i] / lengths[i]) * abs_array) for i in range(3)])

    return lengths * (np.arctan2(-zeta_bar, -xi_bar) + np.pi) / (2*np.pi) 

def shift_periodic_mean_to_center(array, lengths):
    """Shift the array periodically so that the pseudo center of mass is centered in the array"""
    dims = array.shape
    center = lengths/2

    r = calc_periodic_mean(array, lengths)

    shift_vec = center - r
    shift_ind = length_to_index(shift_vec, array, lengths)

    return np.roll(array, shift=shift_ind, axis=(0,1,2))

def get_periodic_com(array, lengths):
    """Compute the center of mass position vector under periodic boundary conditions.

    This function follows section IV of: https://pubs.aip.org/aip/jcp/article/162/20/204103/3347482/On-the-estimation-of-center-of-mass-in-periodic"""
    dims = array.shape
    d = lengths/dims
    center = lengths/2

    r = calc_periodic_mean(array, lengths)
    shift_vec = center - r
    
    shifted_array = np.abs(shift_periodic_mean_to_center(array, lengths))**2
    centered_com = center_of_mass(shifted_array) * d

    return centered_com - shift_vec

def periodic_distance(grid, length_comp):
    """Calculate the distance between two points under periodic boundary conditions"""
    return np.mod(grid + length_comp/2, length_comp) - length_comp/2

def get_wf_spread(wf, lengths):
    """Calculate the length of the vector of NORMALIZED standard deviations for the wavefunction under periodic boundary conditions."""
    dims = wf.shape
    d = lengths/dims
    grids = np.meshgrid(*[np.arange(d) for d in dims], indexing="ij")

    prob = np.abs(wf)**2

    center = get_periodic_com(wf, lengths)
    sep = np.array([periodic_distance(grids[i] * d[i] - center[i], lengths[i]) for i in range(3)])

    sigma_vec = np.sqrt(np.array([np.sum(prob * sep[i]**2) for i in range(3)])) / lengths

    return LA.norm(sigma_vec)

#########################################################################################################

