

cd ann_mwrapper
ann_compile_mex
cd ..

mex mexSource/mexComputeFeature.cpp -output mex/mexComputeFeature
mex mexSource/mexComputeFeature3D.cpp -output mex/mexComputeFeature3D
mex mexSource/mexTensorMatching.cpp -output mex/mexTensorMatching


