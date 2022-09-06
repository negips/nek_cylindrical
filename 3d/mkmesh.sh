#!/bin/sh

genbox << EOF
section2d.box
EOF
#
reatore2 << EOF
box
three
EOF
#
rm -f box.rea
rm -f three.rea
genmap << EOF
three
0.00001
EOF

