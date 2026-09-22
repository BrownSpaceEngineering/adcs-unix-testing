# ADCS Unix Testing

## Build 
```bash
git clone https://github.com/BrownSpaceEngineering/adcs-unix-testing.git
cd adcs-unix-testing
git submodule update --init --recursive
cmake -S . -B build
cmake --build build -j4
```
