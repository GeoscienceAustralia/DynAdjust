#!/usr/bin/env -S uv run --script
# /// script
# dependencies = []
# ///
"""Convert DynaML XML station and measurement files to JSONL format.

Usage:
    xml2jsonl.py input.xml [input2.xml ...] [-o output.jsonl]
    xml2jsonl.py stations.xml measurements.xml  # auto-names .jsonl

If -o is not given, output is written to input.jsonl (replacing .xml).
"""
import argparse
import json
import sys
import xml.etree.ElementTree as ET


def text(elem, tag, default=""):
    """Get text content of a child element."""
    child = elem.find(tag)
    if child is None:
        return default
    return child.text if child.text else default


def num(s):
    """Convert a string to a JSON number if possible, else return as string.

    Preserves exact decimal representation by using float for values that
    round-trip losslessly, otherwise keeps the original string.
    """
    if not s:
        return s
    try:
        f = float(s)
        # Verify round-trip: float -> str -> float matches
        if f == float(repr(f)):
            return f
        return s
    except (ValueError, OverflowError):
        return s


def numtext(elem, tag, default=""):
    """Get text of a child element and convert to number if possible."""
    v = text(elem, tag, default)
    return num(v)


# Measurement types that use Second station
TYPES_WITH_SECOND = set("ABKCEMSDGLVZX")
# Measurement types that use Third station
TYPES_WITH_THIRD = set("A")
# Measurement types that use Value/StdDev
TYPES_WITH_VALUE = set("ABKCEMDHLVZRIJPQS")
# Measurement types that use InstHeight/TargHeight
TYPES_WITH_HEIGHTS = set("SVZ")


def convert_station_file(xmlpath, jsonlpath):
    """Convert a DynaML station XML file to JSONL."""
    tree = ET.parse(xmlpath)
    root = tree.getroot()
    count = 0

    with open(jsonlpath, "w") as f:
        # Header
        hdr = {
            "DnaXmlFormat": {
                "type": root.get("type", ""),
                "referenceframe": root.get("referenceframe", ""),
                "epoch": root.get("epoch", ""),
            }
        }
        f.write(json.dumps(hdr, separators=(",", ":")) + "\n")

        for stn_elem in root.findall("DnaStation"):
            coord_elem = stn_elem.find("StationCoord")
            coord = {}
            if coord_elem is not None:
                coord = {
                    "Name": text(coord_elem, "Name"),
                    "XAxis": numtext(coord_elem, "XAxis"),
                    "YAxis": numtext(coord_elem, "YAxis"),
                    "Height": numtext(coord_elem, "Height"),
                }
                hz = text(coord_elem, "HemisphereZone")
                if hz:
                    coord["HemisphereZone"] = hz

            stn = {
                "Name": text(stn_elem, "Name"),
                "Constraints": text(stn_elem, "Constraints"),
                "Type": text(stn_elem, "Type"),
                "StationCoord": coord,
                "Description": text(stn_elem, "Description"),
            }
            f.write(json.dumps({"DnaStation": stn}, separators=(",", ":")) + "\n")
            count += 1

    return count


def parse_covariance(cov_elem):
    """Parse a GPSCovariance or PointCovariance element."""
    return {
        "m11": numtext(cov_elem, "m11"),
        "m12": numtext(cov_elem, "m12"),
        "m13": numtext(cov_elem, "m13"),
        "m21": numtext(cov_elem, "m21"),
        "m22": numtext(cov_elem, "m22"),
        "m23": numtext(cov_elem, "m23"),
        "m31": numtext(cov_elem, "m31"),
        "m32": numtext(cov_elem, "m32"),
        "m33": numtext(cov_elem, "m33"),
    }


def parse_baseline_elem(child):
    """Parse a GPSBaseline XML element into a dict."""
    bsl = {
        "X": numtext(child, "X"),
        "Y": numtext(child, "Y"),
        "Z": numtext(child, "Z"),
        "SigmaXX": numtext(child, "SigmaXX"),
        "SigmaXY": numtext(child, "SigmaXY"),
        "SigmaXZ": numtext(child, "SigmaXZ"),
        "SigmaYY": numtext(child, "SigmaYY"),
        "SigmaYZ": numtext(child, "SigmaYZ"),
        "SigmaZZ": numtext(child, "SigmaZZ"),
    }
    msr_id = text(child, "MeasurementID")
    if msr_id:
        bsl["MeasurementID"] = int(msr_id)
    covs = [parse_covariance(c) for c in child.findall("GPSCovariance")]
    if covs:
        bsl["GPSCovariance"] = covs
    return bsl


def parse_clusterpoint_elem(child):
    """Parse a Clusterpoint XML element into a dict."""
    pt = {
        "X": numtext(child, "X"),
        "Y": numtext(child, "Y"),
        "Z": numtext(child, "Z"),
        "SigmaXX": numtext(child, "SigmaXX"),
        "SigmaXY": numtext(child, "SigmaXY"),
        "SigmaXZ": numtext(child, "SigmaXZ"),
        "SigmaYY": numtext(child, "SigmaYY"),
        "SigmaYZ": numtext(child, "SigmaYZ"),
        "SigmaZZ": numtext(child, "SigmaZZ"),
    }
    msr_id = text(child, "MeasurementID")
    if msr_id:
        pt["MeasurementID"] = int(msr_id)
    covs = [parse_covariance(c) for c in child.findall("PointCovariance")]
    if covs:
        pt["PointCovariance"] = covs
    return pt


def convert_measurement_file(xmlpath, jsonlpath):
    """Convert a DynaML measurement XML file to JSONL."""
    tree = ET.parse(xmlpath)
    root = tree.getroot()
    count = 0

    with open(jsonlpath, "w") as f:
        # Header
        hdr = {
            "DnaXmlFormat": {
                "type": root.get("type", ""),
                "referenceframe": root.get("referenceframe", ""),
                "epoch": root.get("epoch", ""),
            }
        }
        f.write(json.dumps(hdr, separators=(",", ":")) + "\n")

        for msr_elem in root.findall("DnaMeasurement"):
            msr_type = text(msr_elem, "Type")
            m = {"Type": msr_type}

            # Source
            source = text(msr_elem, "Source")
            m["Source"] = source

            # Ignore
            m["Ignore"] = text(msr_elem, "Ignore")

            # Epoch (all types may have it)
            epoch = text(msr_elem, "Epoch")
            if epoch:
                m["Epoch"] = epoch

            obs_epoch = text(msr_elem, "EpochOfObservation")
            if obs_epoch:
                m["EpochOfObservation"] = obs_epoch

            # ReferenceFrame (GPS types)
            if msr_type in ("G", "X", "Y"):
                m["ReferenceFrame"] = text(msr_elem, "ReferenceFrame")

            # MeasurementID and ClusterID (all types may have them)
            msr_id = text(msr_elem, "MeasurementID")
            if msr_id:
                m["MeasurementID"] = int(msr_id)
            cluster_dbid = text(msr_elem, "ClusterID")
            if cluster_dbid:
                m["ClusterID"] = int(cluster_dbid)

            if msr_type in ("G", "X", "Y"):
                m["Vscale"] = numtext(msr_elem, "Vscale")
                m["Pscale"] = numtext(msr_elem, "Pscale")
                m["Lscale"] = numtext(msr_elem, "Lscale")
                m["Hscale"] = numtext(msr_elem, "Hscale")

            if msr_type == "G":
                # G type: single baseline, First/Second at measurement level
                m["First"] = text(msr_elem, "First")
                m["Second"] = text(msr_elem, "Second")
                bsl_elem = msr_elem.find("GPSBaseline")
                if bsl_elem is not None:
                    m["GPSBaseline"] = [parse_baseline_elem(bsl_elem)]

            elif msr_type == "X":
                # X type: multiple baselines, First/Second per baseline
                m["Total"] = int(text(msr_elem, "Total"))
                baselines = []
                cur_first = cur_second = ""
                for child in msr_elem:
                    if child.tag == "First":
                        cur_first = child.text or ""
                        if "First" not in m:
                            m["First"] = cur_first
                    elif child.tag == "Second":
                        cur_second = child.text or ""
                        if "Second" not in m:
                            m["Second"] = cur_second
                    elif child.tag == "GPSBaseline":
                        bsl = parse_baseline_elem(child)
                        bsl["First"] = cur_first
                        bsl["Second"] = cur_second
                        baselines.append(bsl)
                m["GPSBaseline"] = baselines

            elif msr_type == "Y":
                # Y type: clusterpoints with per-point First
                m["Coords"] = text(msr_elem, "Coords")
                m["Total"] = int(text(msr_elem, "Total"))
                points = []
                cur_first = ""
                for child in msr_elem:
                    if child.tag == "First":
                        cur_first = child.text or ""
                        if "First" not in m:
                            m["First"] = cur_first
                    elif child.tag == "Clusterpoint":
                        pt = parse_clusterpoint_elem(child)
                        pt["First"] = cur_first
                        points.append(pt)
                m["Clusterpoint"] = points

            elif msr_type == "D":
                # D type: direction set
                m["First"] = text(msr_elem, "First")
                second = text(msr_elem, "Second")
                if second:
                    m["Second"] = second
                m["Value"] = numtext(msr_elem, "Value")
                m["StdDev"] = numtext(msr_elem, "StdDev")
                m["Total"] = int(text(msr_elem, "Total"))

                directions = []
                for dir_elem in msr_elem.findall("Directions"):
                    d = {
                        "Ignore": text(dir_elem, "Ignore"),
                        "Target": text(dir_elem, "Target"),
                        "Value": numtext(dir_elem, "Value"),
                        "StdDev": numtext(dir_elem, "StdDev"),
                    }
                    directions.append(d)
                m["Directions"] = directions

            else:
                # Scalar types: A, B, C, E, H, I, J, K, L, M, P, Q, R, S, V, Z
                m["First"] = text(msr_elem, "First")
                if msr_type in TYPES_WITH_SECOND:
                    m["Second"] = text(msr_elem, "Second")
                if msr_type in TYPES_WITH_THIRD:
                    m["Third"] = text(msr_elem, "Third")
                if msr_type in TYPES_WITH_VALUE:
                    m["Value"] = numtext(msr_elem, "Value")
                    m["StdDev"] = numtext(msr_elem, "StdDev")
                if msr_type in TYPES_WITH_HEIGHTS:
                    ih = numtext(msr_elem, "InstHeight")
                    if ih:
                        m["InstHeight"] = ih
                    th = numtext(msr_elem, "TargHeight")
                    if th:
                        m["TargHeight"] = th

            f.write(
                json.dumps({"DnaMeasurement": m}, separators=(",", ":")) + "\n"
            )
            count += 1

    return count


def detect_file_type(xmlpath):
    """Detect whether an XML file is a station or measurement file."""
    tree = ET.parse(xmlpath)
    root = tree.getroot()
    file_type = root.get("type", "").lower()

    if "station" in file_type:
        return "station"
    if "measurement" in file_type:
        return "measurement"

    # Fall back to content detection
    if root.find("DnaStation") is not None:
        return "station"
    if root.find("DnaMeasurement") is not None:
        return "measurement"

    return "unknown"


def main():
    parser = argparse.ArgumentParser(
        description="Convert DynaML XML files to JSONL format."
    )
    parser.add_argument("input", nargs="+", help="Input XML file(s)")
    parser.add_argument(
        "-o",
        "--output",
        help="Output JSONL file (only valid with single input)",
    )
    args = parser.parse_args()

    if args.output and len(args.input) > 1:
        print("Error: -o can only be used with a single input file.", file=sys.stderr)
        sys.exit(1)

    for xmlpath in args.input:
        if args.output:
            jsonlpath = args.output
        elif xmlpath.endswith(".xml"):
            jsonlpath = xmlpath[:-4] + ".jsonl"
        else:
            jsonlpath = xmlpath + ".jsonl"

        file_type = detect_file_type(xmlpath)

        if file_type == "station":
            count = convert_station_file(xmlpath, jsonlpath)
            print(f"{jsonlpath}: {count} stations")
        elif file_type == "measurement":
            count = convert_measurement_file(xmlpath, jsonlpath)
            print(f"{jsonlpath}: {count} measurements")
        else:
            print(f"Error: Cannot determine type of {xmlpath}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
