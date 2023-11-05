import numpy as np
import scipy.ndimage as nd
import matplotlib.pyplot as plt
import sys


def get_mask(shape):
    """
    Apply round mask with radius N on square array with dimensions (npix, npix).
    """
    xx, yy = np.meshgrid(
        np.linspace(-1, 1, shape[1], True), 
        np.linspace(-1, 1, shape[0], True), indexing='ij')
    mask = (xx**2 + yy**2) <= 1.
    return mask


def forwardproject(sample, angles):
    """
    Simulate data aquisition in tomography from line projections.
    Forwardproject a given input sample slice to obtain a simulated sinogram.

    Hints
    -----
    Use scipy.ndimage.rotate(..., reshape=False) to simulate the sample
    rotation.
    Use numpy.sum() along one axis to simulate the line projection integral.
    """
    nproj = len(angles)  # calculate number of projections
    sh = sample.shape  # calculate shape of sample

    # define empty sinogram container, angles along x-axis
    sinogram = np.zeros([nproj, sh[1]], dtype=np.float32)

    # loop over all projections
    for i, ang in enumerate(angles):
        sys.stdout.write("\rSimulating: %03i/%i" % (i + 1, nproj))
        sys.stdout.flush()
        im_rot = nd.rotate(sample, ang, reshape=False, order=1)
        sinogram[i] = np.sum(im_rot, axis=0)
    return sinogram


def filter_ramlak(sinogram):
    """
    Filter a given sinogram using a ramp filter

    Hints
    -----
    First define a ramp filter in Fourier domain (you can use np.fft.fftfreq).
    Filter the sinogram in Fourier space using the convolution theorem.
    """
    # sinogram dimensions
    npix = sinogram.shape[1]

    # Generate basic ramp filter (hint: there is the function np.fft.fftfreq.
    # Try it and see what it does. Watch out for a possible fftshift)
    ramp_filter = np.pi * np.abs(np.fft.fftfreq(npix))

    # filter the sinogram in Fourier space in detector pixel direction
    # Use the np.fft.fft along the axis=1
    sino_ft = np.fft.fft(sinogram, axis=1)

    # Multiply the ramp filter onto the 1D-FT of the sinogram and transform it
    # back into spatial domain
    sino_filtered = np.real(np.fft.ifft(sino_ft * ramp_filter, axis=1))

    return sino_filtered


def backproject(sinogram, angles):
    """
    Backproject a given sinogram.
    Hints
    -----
    Perform the backprojection inversely to the way we did the
    forwardprojection, by smearing each projection in the sinogram back along
    the axis that you summed before in forwardproject(),
    then rotating the resulting backprojection
    to get the right backprojection angle.
    Use scipy.ndimage.rotate(..., ..., reshape=False)
    """
    # calculate number of projections, and pixels
    nproj, npix = sinogram.shape

    # define empty container for reconstruction of sample
    reconstruction = np.zeros([npix, npix], dtype=np.float32)

    # loop over all projections
    for i, ang in enumerate(angles):
        sys.stdout.write("\rReconstructing: %03i / %i" % (i + 1, nproj))
        sys.stdout.flush()

        backprojection = np.tile(sinogram[i], (npix, 1))
        backprojection /= nproj  # Just normalization

        rotated_backprojection = nd.rotate(
            backprojection, -ang, reshape=False, order=1)

        # Add the rotated backprojection
        reconstruction += rotated_backprojection

    reconstruction *= get_mask(reconstruction.shape)

    return reconstruction
    
def art(iters, shape, detector_pixel_lenght, fbp, angles, sinogram, pixel_size):
    '''
    We perform an algebraic reconstruction with iters number of iterations.
    '''
    # Prepare mask and renormalization term (Don't worry about this too much)
    mask = get_mask(shape)

    # The number of image pixels contributing to the respective ray (rot. symmetric)
    mask_proj = np.sum(mask, axis=0)

    # Norm is one over the number of image pixels. This accounts for division by 0.
    renorm = np.zeros(detector_pixel_lenght, dtype=np.float32)
    renorm[mask_proj != 0] = 1. / mask_proj[mask_proj != 0]

    # Start with an initial guess of zeros for the tomogram or the FBP
    initial_tomo = fbp.copy()

    # Initialize
    tomo = initial_tomo.copy()
    error = []

    # Main loop over the iterations
    for i in range(iters):
        
        # Forwardproject your tomorgam for the current angle
        proj = forwardproject(tomo, angles)
                
        # Calculate the difference between the forward projection and 
        #the sinogram for the current angle    
        diff = sinogram/pixel_size - proj
            
        # Accumulate the error to the total error
        error.append(np.sum(diff**2))
            
        # Back-project the renormalized difference for the current angle
        # Hint - multiply the difference by renorm calculated above    
        tomo_update = backproject(renorm * diff, angles)
            
        # Update the tomogram
        tomo += tomo_update
        
        print('\rIteration %i completed!'.ljust(80) % (i + 1))
        print('Error: %.2f' % error[-1])
        
    return tomo, error