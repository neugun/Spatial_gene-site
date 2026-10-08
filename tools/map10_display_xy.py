"""Single orientation convention for map10 spatial display layers.

IMPORTANT: Stage631 x_um / y_um are ALREADY x-flipped and y-flipped relative
to raw array pixel coordinates (confirmed against Stage393/Stage838 convention).
Therefore they must be used directly, without applying flip_xy() again.
flip_xy() below is ONLY for genuinely unflipped raw input coordinates.
Never use image-level mirroring for correcting scientific source orientations.
Raw cells, expression and region identities are frozen.
"""
import numpy as np

CONVENTION = "x flipped, y flipped"

def flip_xy(x, y, *, x_bounds=None, y_bounds=None):
    xa=np.asarray(x, dtype=float)
    ya=np.asarray(y, dtype=float)
    if x_bounds is None:
        x_bounds=(float(np.nanmin(xa)),float(np.nanmax(xa)))
    if y_bounds is None:
        y_bounds=(float(np.nanmin(ya)),float(np.nanmax(ya)))
    xmin,xmax=map(float,x_bounds)
    ymin,ymax=map(float,y_bounds)
    assert np.isfinite([xmin,xmax,ymin,ymax]).all()
    assert xmax>xmin and ymax>ymin
    return xmin+xmax-xa, ymin+ymax-ya

def verify_flip_xy(x_raw,y_raw,x_display,y_display,*,x_bounds=None,y_bounds=None,atol=1e-6):
    expect_x,expect_y=flip_xy(x_raw,y_raw,x_bounds=x_bounds,y_bounds=y_bounds)
    assert np.allclose(x_display,expect_x,rtol=0,atol=atol),"x flipped transform mismatch"
    assert np.allclose(y_display,expect_y,rtol=0,atol=atol),"y flipped transform mismatch"
    return True


def from_stage631_preflipped(x_um, y_um):
    """Stage631 pre-flipped canonical plotting frame, no additional flip."""
    x=np.asarray(x_um,dtype=float)
    y=np.asarray(y_um,dtype=float)
    assert x.shape==y.shape and np.isfinite(x).all() and np.isfinite(y).all()
    return x.copy(),y.copy()
