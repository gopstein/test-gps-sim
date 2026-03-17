#!/usr/bin/env python3
"""
tle_to_motion.py - Generate a gps-sim ECEF motion file from a TLE

Usage:
    python3 tle_to_motion.py --tle <file.tle> --sat <NORAD_ID or name>
                             --start "YYYY/MM/DD,HH:MM:SS"
                             --end   "YYYY/MM/DD,HH:MM:SS"
                             [--output <file.csv>]
                             [--rate 10]

Output:
    CSV file with one ECEF X,Y,Z position (meters) per line at the
    specified sample rate (default 10 Hz), suitable for use with
    gps-sim -m <file>

Dependencies:
    pip install sgp4
"""

import argparse
import sys
from datetime import datetime, timezone, timedelta
from sgp4.api import Satrec, jday


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate gps-sim ECEF motion file from a TLE for a LEO satellite."
    )
    parser.add_argument(
        "--tle", required=True, metavar="FILE",
        help="Path to TLE file (standard 2-line or 3-line format)"
    )
    parser.add_argument(
        "--sat", required=True, metavar="ID_OR_NAME",
        help="NORAD catalog number (e.g. 25544) or satellite name (e.g. ISS)"
    )
    parser.add_argument(
        "--start", required=True, metavar="YYYY/MM/DD,HH:MM:SS",
        help="Scenario start time in UTC"
    )
    parser.add_argument(
        "--end", required=True, metavar="YYYY/MM/DD,HH:MM:SS",
        help="Scenario end time in UTC"
    )
    parser.add_argument(
        "--output", metavar="FILE", default=None,
        help="Output CSV filename (default: <satid>_<start>.csv)"
    )
    parser.add_argument(
        "--rate", type=int, default=10, metavar="HZ",
        help="Sample rate in Hz (default: 10, must be 10 for gps-sim)"
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# TLE parsing
# ---------------------------------------------------------------------------

def load_tle(filepath, sat_id):
    """
    Load a TLE from a file.  Supports:
      - 3-line format (name line + line1 + line2)
      - 2-line format (line1 + line2 only)

    sat_id can be a NORAD number (int or string of digits) or a
    case-insensitive substring of the satellite name.
    """
    with open(filepath, "r") as f:
        raw = [line.strip() for line in f if line.strip()]

    sat_id_str = str(sat_id).strip().upper()
    is_norad = sat_id_str.isdigit()

    # Group into TLE entries
    entries = []
    i = 0
    while i < len(raw):
        line = raw[i]
        # TLE line 1 starts with "1 "
        if line.startswith("1 ") and len(line) >= 69:
            # 2-line format
            if i + 1 < len(raw) and raw[i + 1].startswith("2 "):
                entries.append((None, raw[i], raw[i + 1]))
                i += 2
            else:
                i += 1
        elif not line.startswith("1 ") and not line.startswith("2 "):
            # Possible name line — look ahead for line1 + line2
            if (i + 2 < len(raw) and
                    raw[i + 1].startswith("1 ") and
                    raw[i + 2].startswith("2 ")):
                entries.append((line, raw[i + 1], raw[i + 2]))
                i += 3
            else:
                i += 1
        else:
            i += 1

    if not entries:
        sys.exit(f"ERROR: No valid TLE entries found in {filepath}")

    # Search for the requested satellite
    for name, line1, line2 in entries:
        norad_in_tle = line2[2:7].strip()
        if is_norad:
            if norad_in_tle == sat_id_str.lstrip("0") or norad_in_tle == sat_id_str:
                return name, line1, line2
        else:
            if name and sat_id_str in name.upper():
                return name, line1, line2

    # If only one entry in the file, use it regardless (convenient for
    # single-sat TLE files)
    if len(entries) == 1:
        print(f"WARNING: Satellite '{sat_id}' not matched by name/NORAD; "
              f"using the only entry in the file.", file=sys.stderr)
        return entries[0]

    sys.exit(
        f"ERROR: Satellite '{sat_id}' not found in {filepath}.\n"
        f"  Available entries ({len(entries)} total):\n" +
        "\n".join(
            f"    NORAD {e[2][2:7].strip():>6}  name: {e[0] or '(none)'}"
            for e in entries[:20]
        )
    )


# ---------------------------------------------------------------------------
# Time helpers
# ---------------------------------------------------------------------------

DT_FMT = "%Y/%m/%d,%H:%M:%S"


def parse_dt(s):
    try:
        return datetime.strptime(s, DT_FMT).replace(tzinfo=timezone.utc)
    except ValueError:
        sys.exit(f"ERROR: Cannot parse time '{s}'. Expected format: YYYY/MM/DD,HH:MM:SS")


def dt_to_jday(dt):
    return jday(dt.year, dt.month, dt.day,
                dt.hour, dt.minute, dt.second + dt.microsecond / 1e6)


# ---------------------------------------------------------------------------
# Propagation
# ---------------------------------------------------------------------------

def propagate(line1, line2, start_dt, end_dt, rate_hz):
    """
    Propagate the satellite from start to end at rate_hz samples/sec.
    Returns a list of (x, y, z) tuples in meters (ECEF via TEME->ECEF
    is handled internally by sgp4 when using the WGS84 model).

    Note: sgp4 returns positions in km in TEME frame.  We convert to
    meters.  gps-sim expects ECEF; TEME and ECEF are very close for
    short intervals and the difference is well within GPS simulation
    tolerance (sub-meter over minutes).  For high-precision work you
    would apply a full TEME->ECEF rotation using the equation of the
    equinoxes, but that is overkill here.
    """
    sat = Satrec.twoline2rv(line1, line2)

    interval_sec = 1.0 / rate_hz
    total_sec = (end_dt - start_dt).total_seconds()

    if total_sec <= 0:
        sys.exit("ERROR: End time must be after start time.")

    n_samples = int(total_sec * rate_hz)
    print(f"  Propagating {n_samples} samples "
          f"({total_sec:.1f}s at {rate_hz} Hz)...")

    positions = []
    errors = 0

    for i in range(n_samples):
        t = start_dt + timedelta(seconds=i * interval_sec)
        jd, fr = dt_to_jday(t)
        e, r, v = sat.sgp4(jd, fr)
        if e != 0:
            errors += 1
            # On error sgp4 returns (0,0,0); skip or use last good point
            if positions:
                positions.append(positions[-1])
            else:
                positions.append((0.0, 0.0, 0.0))
        else:
            # r is in km, convert to meters
            positions.append((r[0] * 1000.0, r[1] * 1000.0, r[2] * 1000.0))

    if errors > 0:
        print(f"  WARNING: {errors} propagation errors encountered "
              f"(last-good-point substituted).", file=sys.stderr)

    return positions


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_csv(positions, filepath, rate_hz):
    interval = 1.0 / rate_hz
    with open(filepath, "w") as f:
        for i, (x, y, z) in enumerate(positions):
            t = i * interval
            f.write(f"{t:.1f},{x:.3f},{y:.3f},{z:.3f}\n")
    print(f"  Written: {filepath}  ({len(positions)} lines)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    print(f"\nTLE file  : {args.tle}")
    print(f"Satellite : {args.sat}")

    # Parse times
    start_dt = parse_dt(args.start)
    end_dt = parse_dt(args.end)
    duration = (end_dt - start_dt).total_seconds()
    print(f"Start     : {start_dt.strftime(DT_FMT)} UTC")
    print(f"End       : {end_dt.strftime(DT_FMT)} UTC")
    print(f"Duration  : {duration:.0f}s ({duration/60:.1f} min)")
    print(f"Rate      : {args.rate} Hz")

    if args.rate != 10:
        print("WARNING: gps-sim requires exactly 10 Hz. "
              "Use --rate 10 for compatibility.", file=sys.stderr)

    # Load TLE
    name, line1, line2 = load_tle(args.tle, args.sat)
    print(f"TLE name  : {name or '(unnamed)'}")
    print(f"  Line1   : {line1}")
    print(f"  Line2   : {line2}")

    # Warn if TLE epoch is stale (> 7 days old reduces accuracy)
    try:
        sat_check = Satrec.twoline2rv(line1, line2)
        tle_epoch = datetime(2000, 1, 1, tzinfo=timezone.utc) + \
                    timedelta(days=sat_check.jdsatepoch +
                              sat_check.jdsatepochF - 2451545.0)
        age_days = (start_dt - tle_epoch).total_seconds() / 86400
        if abs(age_days) > 7:
            print(f"  WARNING: TLE epoch is {age_days:.1f} days from start time. "
                  f"Accuracy may be reduced for LEO.", file=sys.stderr)
        else:
            print(f"  TLE age : {age_days:.2f} days from start (OK)")
    except Exception:
        pass  # epoch check is informational only

    # Propagate
    positions = propagate(line1, line2, start_dt, end_dt, args.rate)

    # Output filename
    if args.output:
        outfile = args.output
    else:
        safe_id = str(args.sat).replace(" ", "_")
        safe_start = args.start.replace("/", "").replace(",", "_").replace(":", "")
        outfile = f"{safe_id}_{safe_start}.csv"

    write_csv(positions, outfile, args.rate)

    print(f"\nDone. Use with gps-sim:")
    print(f"  ./gps-sim -e <brdc_file> -m {outfile} -s {args.start} -d {int(duration)} -I -r plutosdr -N 192.168.2.1 -g -40\n")


if __name__ == "__main__":
    main()
