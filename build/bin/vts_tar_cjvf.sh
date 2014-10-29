find . \( -name '*.vts' -o -name '*.pvts' \) | 
xargs rm -f ; snac2vtk . 1 10000000 ; find . \( -name '*.vts' -o -name '*.pvts' \) -print | 
tar cjvf lowc2.tar.bz2 --files-from - ; find . \( -name '*.vts' -o -name '*.pvts' \) | 
xargs rm -f