#!/bin/bash
# Check <Source> tags in exported DynaML XML. Exits 1 on any mismatch.
[ $# -lt 2 ] && { echo "Usage: $0 <xml_file> <source> ..."; exit 1; }
f="$1"; shift
[ -f "$f" ] || { echo "FAIL: $f not found"; exit 1; }
r=0
for s in "$@"; do
    [[ $s == EMPTY ]] && pat='<Source></Source>|<Source */>' || pat="<Source>$s</Source>"
    grep -Eq "$pat" "$f" && echo "PASS: $s in $f" || { echo "FAIL: $s not in $f"; r=1; }
done
exit $r
