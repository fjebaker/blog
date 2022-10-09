

```
export FC
export CC
export CXX

./configure --enable-xs-models-only --enable-mac_arm_build > config.txt 2>&1
```

```
HEADAS=/Users/lx21966/Developer/xspec/heasoft-6.30.1/aarch64-apple-darwin21.5.0 julia --project=. ./loading.jl
```