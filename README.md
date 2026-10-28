# Communications libraries and blocks

This communications toolkit provides DSP and higher level communications
blocks such as filters, basic arithmetic, waveform generators, symbol
manipulation, encoding/decoding, frame synchronization, error correction,
and MAC layer support.

This branch (`pothos2`) builds against the modernized core
(https://github.com/markvin/PothosCore, branch `pothos2`). The modules use
only the public API of the core.

## Dependencies

- an installed Pothos core (`find_package(Pothos)`);
- nlohmann_json 3.11 or higher and xsimd 14.2 or higher (both are headers
  only).

CMake finds the dependencies with `find_package`. They can come from the
system packages or from Conan 2 (`conanfile.txt`).

The filter and window blocks use the filter design code in
`filter/design/`. This code comes from spuce (https://github.com/hayguen/spuce,
commit `ad08df7`). Comms does not need spuce.

## Build

With system packages (for example Arch Linux: `nlohmann-json xsimd`):

```sh
cmake -S . -B build -G Ninja -DCMAKE_PREFIX_PATH=<core prefix>
cmake --build build
ctest --test-dir build --output-on-failure
```

With Conan 2:

```sh
conan install . --output-folder=build/conan -s build_type=Release -s compiler.cppstd=23 --build=missing
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_TOOLCHAIN_FILE=build/conan/conan_toolchain.cmake \
    -DCMAKE_PREFIX_PATH=<core prefix>
cmake --build build
ctest --test-dir build --output-on-failure
```

`<core prefix>` is the installation prefix of the Pothos core. Do not set
it when the core is installed in a default location of the system.

Install the modules:

```sh
cmake --install build
```

The core loads modules only from its own prefix (and from
`POTHOS_PLUGIN_PATH`). Thus the default install prefix of this project is
the prefix of the core. When you give `CMAKE_INSTALL_PREFIX`, set
`POTHOS_PLUGIN_PATH` to `<prefix>/lib/Pothos/modules<ABI version>`.

## Tests

The self-tests of comms use blocks of PothosBlocks (for example
`/blocks/feeder_source` and `/blocks/collector_sink`). Before you run the
tests, install PothosBlocks (https://github.com/markvin/PothosBlocks,
branch `pothos2`) into the prefix of the core. The build of comms does not
need PothosBlocks.

Each module lists its self-tests (`TESTS` of `pothos_comms_module()`).
The test `comms.stage` installs the modules into `build/stage`. Then
`comms.selftests` runs each self-test with the PothosUtil of the core and
`POTHOS_PLUGIN_PATH` set to the staged modules.

`comms.selftests` stops with an error in these conditions:

- PothosBlocks is not installed in the prefix of the core;
- the comms modules are also installed in the prefix of the core (the
  core would load the two copies). Remove the installed copy, or test
  against another installation of the core.

## Documentation

* https://github.com/pothosware/PothosComms/wiki

## Licensing information

Use, modification and distribution is subject to the Boost Software
License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt). Some files have other licenses, as
marked in their headers: the SIMD kernels in `math/SIMD` are
BSD-3-Clause, and `filter/design/SincBandStop.cpp` is BSD-3-Clause
(Cardinal Peak, LLC).
