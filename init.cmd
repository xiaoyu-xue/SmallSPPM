cd external
cmake -S . -B ./build
cd build
msbuild ext_build.sln /property:Configuration=Debug
msbuild ext_build.sln /property:Configuration=Release