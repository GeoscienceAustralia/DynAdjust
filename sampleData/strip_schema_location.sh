#!/bin/bash
# Strip xsi:noNamespaceSchemaLocation attribute from DynaML XML files.
# Usage: strip_schema_location.sh <input.xml> <output.xml>
[ $# -ne 2 ] && { echo "Usage: $0 <input.xml> <output.xml>"; exit 1; }
[ -f "$1" ] || { echo "FAIL: $1 not found"; exit 1; }
sed 's/ xsi:noNamespaceSchemaLocation="DynaML\.xsd"//' "$1" > "$2"
