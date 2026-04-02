Build the IAMR-style fluid executable with CMake by enabling the `fluid/`
subproject.

From the repo root:

```bash
cmake -S . -B build/release -DELFF_BUILD_FLUID=ON
cmake --build build/release --target iamr
```

For a standalone configure directly in `fluid/`, the local install prefix
`../install/release` is added to `CMAKE_PREFIX_PATH` automatically when it
exists.
