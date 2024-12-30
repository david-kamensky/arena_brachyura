# Arena Brachyura

This is a minimal software-rendered FPS game I wrote as an exercise to try out the internet's "[most admired](https://github.blog/2023-08-30-why-rust-is-the-most-admired-language-among-developers/)" programming language, Rust. It started as a direct port of a smaller "2.5D" game I'd written in college while learning C/C++, but I ended up generalizing it in a number of ways, including 3D mouselook and dynamic lighting. I also overhauled the art with various crab-like monsters inspired by Rust's unofficial mascot, hence "*Brachyura*" (the scientific name for crabs) in the title.

## How to compile and run

Install [Rust](https://www.rust-lang.org/tools/install), [SDL2](https://wiki.libsdl.org/SDL2/Installation), and [SDL2\_image](https://github.com/libsdl-org/SDL_image/releases), clone the repository, `cd` into its root directory, and run `cargo build --release`. (The default debug build is slow to the point of being effectively unplayable, at least on my computer.) The directory `<repository>/target/release/` will contain an executable, `arena_brachyura`. Download an archive of the game data [here](https://www.dropbox.com/scl/fi/sx9hzzt497yiz6mcp3w6y/arena_brachyura_data_7_20_2024.zip?rlkey=wegn0nce782q0ipadcz7namuv&st=thv20x4y&dl=1), unzip it, and move the `data/` directory into the same directory as the executable. The executable can be run from any working directory, as long as it and `data/` are in the same parent directory (i.e., the executable will look for `data/` in its own parent directory).

## Troubleshooting

* It may be necessary to manually set the version of SDL2 in `Cargo.toml` if an older version is installed on your system (e.g., through the package manager of an older version of Ubuntu Linux).
* Some fullscreen rendering issues have been observed on newer versions of Ubuntu Linux when using the default Wayland display server protocol.  (The same issues appear in other fullscreen applications running at less than the desktop resolution, and therefore appear to be system bugs outside the scope of this game.)  The game can be run in a window with the command line parameter `--windowed`, or Wayland can be disabled in favor of X by editing system configuration files (see, e.g., [this tutorial](https://linuxconfig.org/how-to-enable-disable-wayland-on-ubuntu-22-04-desktop)).
