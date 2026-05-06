#!/usr/bin/env python3
"""
motion_to_kml.py - Convert a gps-sim ECEF motion file to KML for Google Earth

Usage:
    python3 motion_to_kml.py --input <motion.csv> [--output <file.kml>]

Input:
    CSV file with one ECEF X,Y,Z position (meters) per line at 10 Hz,
    as produced by tle_to_motion.py:  time,x,y,z

Output:
    KML file viewable in Google Earth showing the orbital path as a
    3D arc at actual altitude.

Dependencies:
    None (stdlib only)
"""

import argparse
import math
import sys
import xml.etree.ElementTree as ET


# ---------------------------------------------------------------------------
# WGS84 ellipsoid parameters
# ---------------------------------------------------------------------------

WGS84_A = 6378137.0              # semi-major axis (meters)
WGS84_F = 1.0 / 298.257223563   # flattening
WGS84_B = WGS84_A * (1.0 - WGS84_F)  # semi-minor axis
WGS84_E2 = 1.0 - (WGS84_B ** 2) / (WGS84_A ** 2)  # first eccentricity squared


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert gps-sim ECEF motion CSV to KML for Google Earth."
    )
    parser.add_argument(
        "--input", required=True, metavar="FILE",
        help="Path to motion CSV file (time,x,y,z in ECEF meters)"
    )
    parser.add_argument(
        "--output", metavar="FILE", default=None,
        help="Output KML filename (default: same basename as input with .kml)"
    )
    parser.add_argument(
        "--downsample", type=int, default=10, metavar="N",
        help="Keep every Nth point (default: 10, i.e. 1 Hz from 10 Hz input)"
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# ECEF to geodetic conversion (iterative Bowring method, WGS84)
# ---------------------------------------------------------------------------

def ecef_to_geodetic(x, y, z):
    """
    Convert ECEF coordinates (meters) to geodetic latitude, longitude (degrees)
    and altitude (meters above WGS84 ellipsoid).
    """
    lon = math.atan2(y, x)

    # Distance from Z-axis
    p = math.sqrt(x ** 2 + y ** 2)

    # Initial estimate of latitude using Bowring's method
    theta = math.atan2(z * WGS84_A, p * WGS84_B)
    e_prime2 = (WGS84_A ** 2 - WGS84_B ** 2) / (WGS84_B ** 2)

    lat = math.atan2(
        z + e_prime2 * WGS84_B * math.sin(theta) ** 3,
        p - WGS84_E2 * WGS84_A * math.cos(theta) ** 3,
    )

    # Iterate for precision (converges in 2-3 iterations)
    for _ in range(5):
        sin_lat = math.sin(lat)
        N = WGS84_A / math.sqrt(1.0 - WGS84_E2 * sin_lat ** 2)
        lat = math.atan2(z + WGS84_E2 * N * sin_lat, p)

    # Altitude
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    N = WGS84_A / math.sqrt(1.0 - WGS84_E2 * sin_lat ** 2)

    if abs(cos_lat) > 1e-10:
        alt = p / cos_lat - N
    else:
        alt = abs(z) - WGS84_B

    return math.degrees(lat), math.degrees(lon), alt


# ---------------------------------------------------------------------------
# CSV reading
# ---------------------------------------------------------------------------

def read_motion_csv(filepath):
    """
    Read a gps-sim motion CSV (time,x,y,z) and return list of (t, x, y, z).
    """
    points = []
    with open(filepath, "r") as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            if len(parts) != 4:
                print(f"  WARNING: skipping line {lineno}: expected 4 fields, "
                      f"got {len(parts)}", file=sys.stderr)
                continue
            try:
                t, x, y, z = float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])
                points.append((t, x, y, z))
            except ValueError:
                print(f"  WARNING: skipping line {lineno}: non-numeric data",
                      file=sys.stderr)
    return points


# ---------------------------------------------------------------------------
# KML generation
# ---------------------------------------------------------------------------

def build_kml(geo_points, name):
    """
    Build a KML document with a styled LineString at actual altitude.
    geo_points: list of (lat, lon, alt_m)
    """
    kml = ET.Element("kml", xmlns="http://www.opengis.net/kml/2.2")
    doc = ET.SubElement(kml, "Document")

    ET.SubElement(doc, "name").text = name

    # Style for the orbit line
    style = ET.SubElement(doc, "Style", id="orbitStyle")
    line_style = ET.SubElement(style, "LineStyle")
    ET.SubElement(line_style, "color").text = "ff0000ff"  # red (aabbggrr)
    ET.SubElement(line_style, "width").text = "3"

    # Placemark with LineString
    pm = ET.SubElement(doc, "Placemark")
    ET.SubElement(pm, "name").text = "Orbital Path"
    ET.SubElement(pm, "styleUrl").text = "#orbitStyle"

    ls = ET.SubElement(pm, "LineString")
    ET.SubElement(ls, "extrude").text = "0"
    ET.SubElement(ls, "tessellate").text = "0"
    ET.SubElement(ls, "altitudeMode").text = "absolute"

    # Build coordinate string: lon,lat,alt (KML order)
    coord_lines = []
    for lat, lon, alt in geo_points:
        coord_lines.append(f"{lon:.6f},{lat:.6f},{alt:.1f}")
    ET.SubElement(ls, "coordinates").text = "\n".join(coord_lines)

    return kml


def write_kml(kml_element, filepath):
    tree = ET.ElementTree(kml_element)
    ET.indent(tree, space="  ")
    tree.write(filepath, xml_declaration=True, encoding="UTF-8")
    print(f"  Written: {filepath}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    print(f"\nInput     : {args.input}")

    # Read motion CSV
    points = read_motion_csv(args.input)
    if not points:
        sys.exit("ERROR: No valid data points found in input file.")
    print(f"  Points  : {len(points)} total")

    # Downsample
    sampled = points[::args.downsample]
    print(f"  Sampled : {len(sampled)} points (every {args.downsample}th)")

    # Convert ECEF to geodetic
    geo_points = []
    for t, x, y, z in sampled:
        lat, lon, alt = ecef_to_geodetic(x, y, z)
        geo_points.append((lat, lon, alt))

    # Sanity check on first point
    lat0, lon0, alt0 = geo_points[0]
    print(f"  First pt: lat={lat0:.4f} lon={lon0:.4f} alt={alt0/1000:.1f} km")

    # Output filename
    if args.output:
        outfile = args.output
    else:
        base = args.input.rsplit(".", 1)[0] if "." in args.input else args.input
        outfile = base + ".kml"

    # Build and write KML
    name = args.input.rsplit("/", 1)[-1].rsplit(".", 1)[0]
    kml = build_kml(geo_points, name)
    write_kml(kml, outfile)

    duration = points[-1][0] - points[0][0]
    print(f"  Duration: {duration:.1f}s ({duration/60:.1f} min)")
    print(f"\nDone. Open {outfile} in Google Earth to view the orbital path.\n")


if __name__ == "__main__":
    main()
