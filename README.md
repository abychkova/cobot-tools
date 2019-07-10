# cobot-tools
Tools for computational biology

# Software required to use cobot-tools

## ViennaRNA

ViennaRNA is a library that already contains lots of useful algorithms that work
with RNA sequences. We use it in order to make tools from cobot-tools better.

Installation guide for Linux can be found [here](https://www.tbi.univie.ac.at/RNA/).

Some issues occur when you try to install ViennaRNA for MacOS. If you want to avoid
them, use the following guide.

### Installation of ViennaRNA for MacOS

Installation of ViennaRNA for MacOS can be done in following steps:

  1. Download library. Here is the [link](https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.13.tar.gz) 
  for version 2.4.13.
  2. Unpack it.
  3. Go to the directory, where archive has been unpacked.
  4. Run command `./configure --without-perl --without-python`.
  5. Remove field `libRNA_la_LDFLAGS` and flag `-static` from the file `src/ViennaRNA/Makefile.am`.
  6. Due to us having changed .am-file, it is needed to install following packages: autoconf and automake. In order to do it, run command: `brew install autoconf automake`.
  7. Run in the root of library command `autoreconf`.
  8. Compile library using command `make -j8` (`j8` - flag that allows us to compile code using 8 jobs).
  9. Run command `make install`.
  
# Tools of cobot-tools

## Sequence.Primer.Optimisation.designPrimer

Given 'DNA' sequence and position in that sequence designs forward primer for that sequence. 
Primer will start at the given position.


