"""Single orientation convention for map10 spatial display layers.

Every 2D tissue map uses independent X and Y flips from the section's raw
physical coordinate bounds. Never use matplotlib axis inversion or pixel rotations.
This is a *display* transform; raw cell centroids, expression and IDs stay frozen.
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
