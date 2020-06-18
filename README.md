# Hilbert Curves in C++

This is a single-header (and `.cpp`) C++17 library for creating Hilbert Curves from Z-Order (Morton) Curves.

Right now, there is only an implementation for 3D Hilbert curves.
You will find many other implementations for 2D variants though.

Note that in order to convert to- and from 3D-points, you will also need to compute a morton index for 3D-points.
You can do so using the `ileave3` function from my [bitmanip](https://github.com/Eisenwave/bitmanip) library.
You can convert from morton indices using `dileave3`.

## License

See the `LICENSE` file.
