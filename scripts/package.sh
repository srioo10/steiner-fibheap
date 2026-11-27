#!/usr/bin/env bash
# Builds the two submission files in the repository root:
#   CS23I1006_CS23I1002.pdf  the report
#   CS23I1006_CS23I1002.zip  README, Makefile, code, instances, results, scripts, report sources
# then unpacks the zip into a temporary folder and checks that it builds and runs there.
set -euo pipefail
cd "$(dirname "$0")/.."
NAME=CS23I1006_CS23I1002
make -s report
cp report/report.pdf "$NAME.pdf"
rm -f "$NAME.zip"
STAGE=$(mktemp -d)
mkdir -p "$STAGE/$NAME"
cp -r README.md Makefile include src tests scripts instances results "$STAGE/$NAME/"
mkdir -p "$STAGE/$NAME/report"
cp -r report/report.tex report/sections report/generated "$STAGE/$NAME/report/"
cp "$NAME.pdf" "$STAGE/$NAME/report/report.pdf"
(cd "$STAGE" && zip -qr "$OLDPWD/$NAME.zip" "$NAME")
rm -rf "$STAGE"

# smoke test from the archive alone
CHECK=$(mktemp -d)
unzip -q "$NAME.zip" -d "$CHECK"
(cd "$CHECK/$NAME" && make -s all && make -s test > /dev/null &&
  for a in greedy1_kmb greedy2_sph greedy3_zelikovsky dp_dreyfus_wagner bnb_steiner; do
    ./build/$a instances/steinlib/B/b01.stp | cut -d, -f1,2,6,7,10,11
  done)
rm -rf "$CHECK"
ls -lh "$NAME.pdf" "$NAME.zip"
