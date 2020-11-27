# ulc-codec
ulc-codec (Ultra-Low-Complexity Codec) is intended to be a low-bitrate audio codec, providing ultra-low-complexity decoding.

## Getting started

### Prerequisites
None (so far). Perhaps GCC if building from source.

### Installing
Run ```make all``` to build the file-based encoding and decoding tools (```ulcencode``` and ```ulcdecode```).

You could also ```make encodetool``` or ```make decodetool```.

## Usage
For the time being, both encoding and decoding tools operate on raw 16-bit audio (with interleaved channels).
Work is being considered for implementing reading/writing of more common formats.

Additionally, the core encoding/decoding routines can theoretically work with any data they are fed, allowing for easier integration with non-file-based blocks of audio in the future.

### Encoding
```ulcencodetool Input.raw Output.ulc RateHz RateKbps [-nc:1] [-blocksize:2048]```

This will take ```Input.raw``` (with a playback rate of ```RateHz```) and encode it into the output file ```Output.ulc```, at a coding rate of ```RateKbps```. ```-nc:X``` sets the number of channels, ```-blocksize:X``` sets the size of each block (ie. the number of coefficients per block).

### Decoding
```ulcdecodetool Input.ulc Output.raw```

This will take ```Input.ulc``` and output ```Output.raw```.

## Possible issues
* Syntax is flexible enough to cause buffer overflows.
* No block synchronization (if an encoded file is damaged, there is no way to detect where the next block lies)
    * It should be possible to prepend each block with the two-byte nybble sequence ```0h,0h,0h,0h```. Such a sequence could only happen at the start of a block (set quantizer to 2<sup>0</sup>, followed by three zero coefficients) and never in any other place (as four zero coefficients would be coded as ```8h,1h```), avoiding false-positives.
* The psychoacoustic model used is somewhat bare-bones, so as to avoid extra complexity and memory usage. As an example, blocks are processed with no memory of prior blocks, which could cause some inefficiency in coding (such as not taking advantage of temporal masking effects). However, it does appear to work very well for what it *does* do.
* Encode/decode tools compile with SSE+SSE2/AVX+AVX2/FMA enabled by default. If the encoder crashes/doesn't work, change these flags in the ```Makefile```.

## Technical details
* Target bitrate: 32..256kbps+ (44.1kHz, M/S stereo)
    * No hard limits on playback rate or coding bitrate
* MDCT-based encoding (using sine window)
    * Window switching is combined with so-called 'overlap switching', the latter of which varies the size of the overlap segment of transient \[sub]blocks. The idea is to center the transient within a subblock, at which point overlap switching takes over to emphasize it, without having to switch to use small windows for the entire block, overall resulting in improved quality compared to the more-common '1 long block or N short blocks' strategy.
* Non-linear coefficient quantization for greater control over dynamic range
* Extremely simple nybble-based syntax (no entropy-code lookups needed)

## Authors
* **Ruben Nunez** - *Initial work* - [Aikku93](https://github.com/Aikku93)

## License
ulc-codec is released under the GPLv3 license. See the LICENSE file for full terms.

## Acknowledgements
* Huge thanks to Dennis K (`DekuTree64`) and `musicalman` for listening to (and for playing programming-rubber-ducky to) all my rants and thought processes as I worked through understanding audio codec design, as well as testing out experimental builds and offering immensely helpful critiques. Without them, this codec would never have gotten as far as it has.

## Gameboy Advance player

As a proof of concept of the decoding complexity, a Gameboy Advance demonstration may be found in the [ulcplayer-gba](https://github.com/Aikku93/ulcplayer-gba) repository.
