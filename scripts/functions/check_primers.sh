#!usr/bin/env bash
#Check if primers include Is instead of Ns
#usage source check_primers.sh "${primers_file}"


if echo "$1" | grep  -q I; then echo
echo "found I , replacing them with Ns"
NEWFASTA=$(echo "$1" | sed 's/I/N/g')

primers_file="${NEWFASTA}"

cat "${primers_file}"

else
  echo "primers look good"
fi
