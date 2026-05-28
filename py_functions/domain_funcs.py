import numpy as np


def rotated_box(lon0, lat0, width, height, angle_deg):
    """Return the four corners of a rotated rectangular polygon centered at (lon0, lat0).

    Parameters
    ----------
    lon0, lat0   : center longitude and latitude
    width, height: full width and height in degrees
    angle_deg    : counter-clockwise rotation angle in degrees
    """
    angle = np.deg2rad(angle_deg)
    dx, dy = width / 2, height / 2
    corners = np.array([[-dx, -dy], [dx, -dy], [dx, dy], [-dx, dy]])
    R = np.array([[np.cos(angle), -np.sin(angle)],
                  [np.sin(angle),  np.cos(angle)]])
    rotated = corners @ R.T
    rotated[:, 0] += lon0
    rotated[:, 1] += lat0
    return rotated


def regional_weighted_mean(da, regions):
    """Area-weighted mean of da for all regions, returning a DataArray with a 'region' dim.

    Parameters
    ----------
    da      : xr.DataArray with 'lat' and 'lon' coordinates
    regions : regionmask.Regions object

    Returns
    -------
    xr.DataArray with a named 'region' coordinate; use .sel(region=name) to extract one.
    """
    mask3d = regions.mask_3D(da)
    weights = np.cos(np.deg2rad(da.lat))
    result = da.where(mask3d).weighted(weights).mean(("lat", "lon"))
    return result.assign_coords(region=mask3d.names)


def nam_regions():
    """Return the standard NAM south/north domain polygons and regionmask object.

    Returns
    -------
    coords_map : dict with keys 'south' and 'north', values are (4,2) corner arrays
                 in (-180:180) longitude convention — ready for matplotlib patches
    regions    : regionmask.Regions with names ['south', 'north']
                 coordinates are in (0:360) convention for regionmask masking
    """
    coords_map = {}
    coords_calc = {}

    # North: axis-aligned box
    coords_map['north'] = rotated_box(-110.5, 32.5, 8, 6, angle_deg=0)
    coords_calc['north'] = [((lon + 360) % 360, lat) for lon, lat in coords_map['north']]

    # South: tilted bottom from a rotated 35° box, horizontal top flush with north's southern edge
    _south_rot = rotated_box(-107, 25, 4.5, 10, angle_deg=35)
    coords_map['south'] = np.array([
        _south_rot[0],           # rotated SW corner
        _south_rot[1],           # rotated SE corner
        coords_map['north'][1],  # top-right = north's SE corner
        coords_map['north'][0],  # top-left  = north's SW corner
    ])
    coords_calc['south'] = [((lon + 360) % 360, lat) for lon, lat in coords_map['south']]

    import regionmask
    regions = regionmask.Regions(
        outlines=[coords_calc['south'], coords_calc['north']],
        names=["south", "north"],
        abbrevs=["south", "north"],
    )
    return coords_map, regions
