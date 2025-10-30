#!/usr/bin/env bash
# Usage: bash summarize_kraken_species.sh input.kraken.species output.summary

set -euo pipefail

input="$1"
output="$2"

v=$((grep '  Afrotheria$' "$input" || echo 0 0) | awk '{{print $2}}')
echo "$v" "$input"

w=$((grep '  Primates$' "$input" || echo 0 0) | awk '{{print $2}}')
echo "$v: $w" "$input"

x=$((grep '  Laurasiatheria$' "$input" || echo 0 0) | awk '{{print $2}}')
echo "$v: $w: $x" "$input"

y=$((grep '  Glires$' "$input" || echo 0 0) | awk '{{print $2}}')
echo "$v: $w: $x: $y" "$input"

m=$((grep '  Mammalia$' "$input" ||echo 0 0) | awk '{{print $2-'$v'-'$w'-'$x'-'$y'}}')
echo " $v: $w: $x: $y: $m : " "$input"

b=$((grep '  Bacteria$' "$input" || echo 0 0) | awk '{{print $2}}')
echo "$v: $w: $x: $y : $m : $b" "$input"

s=$((grep '  Sauropsida$' "$input" || echo 0 0) | awk '{{print $2}}')
echo "$v: $w: $x: $y : $m : $b : $s :" "$input"

r=$(grep '	root$' "$input" | awk '{{print $2-'$m'-'$b'-'$s'-'$v'-'$w'-'$x'-'$y'}}')
echo "$v: $w: $x: $y : $m : $b : $s : $r :" "$input"

u=$(grep '	unclassified$' "$input" | awk '{{print $2}}')
echo "$v: $w: $x: $y : $m : $b : $s : $r : $u : " "$input"

echo "Mammalia $m" >> "$output"
echo "Bacteria $b" >> "$output"
echo "Sauropsida $s" >> "$output"
echo "root $r" >> "$output"
echo "unclassified $u" >> "$output"
echo "Afrotheria $v" >> "$output"
echo "Primates $w" >> "$output"
echo "Laurasiatheria $x" >> "$output"
echo "Glires $y" >> "$output"
