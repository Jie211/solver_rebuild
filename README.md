# rebuild soler now

## MEMO

1. 16/8/12 read optiong from command line
2. 16/8/29 show input data
3. 16/9/17 add vpgcr, vpgmres

## TODO
1. Select matrix dir
2. Select which bvec in matrix to use
4. kskipcr still have bug.
5. CUDA support

## BUG
1. in io.c get_mat_data, S_ISDIR get worng return ???
  - fixed
2. kskipcr still bugging
  - delete kskipcr now
