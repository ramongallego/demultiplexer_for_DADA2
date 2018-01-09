#!usr/bin/env/bash
# a function to retrieve a filename without the path
# useful for messages
# Usage source strip_path.sh file varname
# returns short_file as the object with the new name

#echo $1
#echo $2
suffix=$2

path_to_delete="$(dirname "${1}" | sed 's_/_\\/_g')"
declare short_$suffix=$(echo "$1" | sed "s/$path_to_delete//g")
varname=short_$suffix
