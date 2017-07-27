#!/bin/bash
echo generating...

mkdir -p build/debug
cd build/debug
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-g2" ../..
cd -

mkdir -p build/release
cd build/release
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3" ../..
cd -

if [[ "$OSTYPE" == "darwin"* ]]; then
  mkdir -p build/xcode
  cd build/xcode
  cmake -G Xcode ../..
fi

echo "done..."
