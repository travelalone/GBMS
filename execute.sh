
awk 'NR>3{print $1,$2,$3}' atom.vert > coordinate
awk 'NR>3{print $1,$2,$3}' atom.face > triangle
