import numpy as np
from astropy.io import fits
from scipy.interpolate import RectBivariateSpline

def index_PSF(PSF, x, y):
    return PSF[x+y*3]

def get_native_PSF(filt, x, y, the_path):
    x = float(np.clip(x, 0, 1014))
    y = float(np.clip(y, 0, 1014))

    f = fits.open("%sPSFSTD_WFC3IR_%s.fits" % (the_path, filt))
    PSF = f[0].data
    f.close()

    if x < 507:
        sx = x/507.
        minx = 0
    else:
        sx = (x - 507.)/507.
        minx = 1
        
    if y < 507:
        sy = y/507.
        miny = 0
    else:
        sy = (y - 507.)/507.
        miny = 1

    out_PSF = 0.
    for dx in [0, 1]:
        for dy in [0, 1]:
            this_x = minx + dx
            this_y = miny + dy
            this_w = (sx*(dx == 1) + (1 - sx)*(dx == 0))*(sy*(dy == 1) + (1 - sy)*(dy == 0))
            print ("x", x, "y", y, "this_x", this_x, "this_y", this_y, "this_w", this_w)
            out_PSF += index_PSF(PSF, x = this_x, y = this_y)*this_w
    return out_PSF

def get_sampled_PSF(filt, x, y, subsample, the_path = "./"):
    native_PSF = get_native_PSF(filt, x, y, the_path)
    orig_sub = np.arange(len(native_PSF), dtype=np.float64)*0.25
    orig_sub -= np.median(orig_sub)

    ifn = RectBivariateSpline(orig_sub, orig_sub, native_PSF, kx = 3, ky = 3, s=0)
    new_sub = np.arange(len(native_PSF)*subsample/4., dtype=np.float64)/subsample
    new_sub -= np.median(new_sub)

    return ifn(new_sub, new_sub)

