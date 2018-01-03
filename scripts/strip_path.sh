path_to_delete="$(dirname "${0}" | sed 's_/_\\/_g')"
echo "$0"
echo "${path_to_delete}"

#newpath=$(echo "$path_to_delete" | sed 's_/_\\/_g')

echo "$0" | sed "s/$path_to_delete//g"
