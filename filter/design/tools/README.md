# Reference data generator

`GenerateReferenceData.cpp` makes `filter/design/ReferenceData.hpp`. The
data comes from spuce `ad08df7` (https://github.com/hayguen/spuce). Before
its own filter design code, comms used this version of spuce. The tests of
the filter design code (`filter/design/TestFilterDesign.cpp`) compare their
results with this data.

The build of comms does not compile this program. Use it only when the
case set changes. Commit the new `ReferenceData.hpp` with the change.

## Procedure

1. Build and install spuce `ad08df7` outside of the comms repository:

   ```sh
   git clone https://github.com/hayguen/spuce.git && cd spuce && git checkout ad08df7
   cmake -S . -B build -G Ninja -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DENABLE_GUI=OFF -DCMAKE_BUILD_TYPE=Release
   cmake --build build && cmake --install build --prefix $PWD/install
   ```

2. Build the generator:

   ```sh
   g++ -std=c++23 -O2 -I<spuce>/install/include filter/design/tools/GenerateReferenceData.cpp \
       <spuce>/install/lib/libspuce.a -o GenerateReferenceData
   ```

3. Run the generator:

   ```sh
   ./GenerateReferenceData > filter/design/ReferenceData.hpp
   ```

## Notes

- spuce writes outside of its buffers for some Remez band edges. The
  generator does not record these cases. It writes each skipped case to
  stderr. `remezGridIsSafe()` tells the condition.
- To find more cases with undefined behavior, build spuce and the generator
  with `-fsanitize=address,undefined`. Then run the generator with
  `ASAN_OPTIONS=detect_leaks=0` (`gaussian_fir` of spuce has a memory leak).
  The output must be the same as the output of the Release build.
- For an unknown IIR type or band, `design_iir` writes a message to
  `std::cout`. The generator keeps these messages out of the header and
  sets `unknown` in the case.
