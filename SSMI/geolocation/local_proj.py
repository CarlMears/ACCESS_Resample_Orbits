import numpy as np
import pyproj

def lat_lon_to_xy(lat, lon, lat0, lon0):
    # Define ellipsoid parameters (WGS84)
    a = 6378137  # semi-major axis in meters
    f = 1 / 298.257223563  # flattening factor
    e_squared = 2 * f - f ** 2
    
    # Convert coordinates to radians
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    lat0_rad = np.radians(lat0)
    lon0_rad = np.radians(lon0)

    # Calculate difference in longitude
    delta_lon = lon_rad - lon0_rad
    delta_lat = lat_rad - lat0_rad
     # Calculate auxiliary values
    e_squared = 2 * f - f ** 2
    rho_EW = a / np.sqrt(1 - e_squared * np.sin(lat0_rad) ** 2)
    rho_NS = a * (1 - e_squared) / ((1 - e_squared * np.sin(lat0_rad) ** 2) ** (3/2))
    
    # Calculate x-coordinate
    x = rho_EW * np.cos(lat0_rad) * np.sin(delta_lon)

    # Calculate y-coordinate
    y = rho_NS * np.tan(lat_rad - lat0_rad)
    y1 = rho_NS * np.sin(delta_lat)

    return x, y, x1, y1

def xy_to_lat_lon_local(x, y, lat0, lon0):

    a = 6378137  # semi-major axis in meters
    f = 1 / 298.257223563  # flattening factor
    e_squared = 2 * f - f ** 2
    # Convert reference latitude to radians
    lat0_rad = np.radians(lat0)
    lon0_rad = np.radians(lon0)

    # Calculate auxiliary values
    e_squared = 2 * f - f ** 2
    rho_EW = a / np.sqrt(1 - e_squared * np.sin(lat0_rad) ** 2)
    rho_NS = a * (1 - e_squared) / ((1 - e_squared * np.sin(lat0_rad) ** 2) ** (3/2))

    # Calculate latitude
    lat_rad = lat0_rad + np.arcsin((y / rho_NS))
    lat = np.degrees(lat_rad)

    # Calculate longitude
    lon_rad = lon0_rad + np.arcsin(x / ((rho_EW)*np.cos(lat0_rad)))
    lon = np.degrees(lon_rad)

    return lat, lon

lat0 = 45.0
lon0 = 34.0

lat = 46.0
lon = 39.0



x, y,x1,y1 = lat_lon_to_xy(lat, lon, lat0, lon0)

print(x, y)


lat1, lon1 = xy_to_lat_lon_local(x, y, lat0, lon0)
print(lat, lon)
print(lat1, lon1)
print()


