# ADCS Unix Testing

## Build 
```
git clone https://github.com/BrownSpaceEngineering/adcs-unix-testing.git
cd adcs-unix-testing
git checkout cmsis-dsp-integrate-finish
git submodule update --init --recursive
cmake -S . -B build
cmake --build build -j4
```
