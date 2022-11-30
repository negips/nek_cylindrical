#!/bin/sh

genbox << EOF
section.box
EOF
#
reatore2 << EOF
box
cyl
EOF
#
rm -f box.rea
rm -f cyl.rea
genmap << EOF
cyl
0.00001
EOF

