# # # #
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/path/to/capstone/install
# # # #

flags="-g -O0"
cmake .. \
  -DCMAKE_INSTALL_PREFIX:PATH="/path/to/install/directory" \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DCMAKE_EXE_LINKER_FLAGS="${flags}" \
  -DSCOREC_CXX_OPTIMIZE=ON \
  -DENABLE_CAPSTONE=ON \
  -DENABLE_FPP=ON
