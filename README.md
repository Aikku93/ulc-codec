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
```ulcencodetool Input.raw Output.ulc RateHz RateKbps [nChan=1]```

This will take ```Input.raw``` (with a playback rate of ```RateHz``` and ```nChan``` channels) and encode it into the output file ```Output.ulc```, at a coding rate of ```RateKbps```.

### Decoding
```ulcdecodetool Input.ulc Output.raw```

This will take ```Input.ulc``` and output ```Output.raw```.

## Possible issues
* Syntax is flexible enough to cause buffer overflows
* No block synchronization (if an encoded file is damaged, there is no way to detect where the next block lies)
    * It should be possible to prepend each block with the two-byte nybble sequence ```0h,0h,0h,0h```. Such a sequence could only happen at the start of a block (set quantizer to 2<sup>0</sup>, followed by three zero coefficients) and never in any other place (as four zero coefficients would be coded as ```8h,1h```), avoiding false-positives.
* Due to the extremely simplified quantization model, a large transform size is almost essential to avoid excessive quality degradation at low bitrates (eg. 32kbps @ 44.1kHz). However, an MP3-esque transform (N=1024; 512 coefficients) will sound only slightly inferior compared to it at 'average' bitrates (eg. 128kbps @ 44.1kHz).
* The psychoacoustic model used is a very crude approximation to the usual models, as I didn't want to overcomplicate it by involving critical bands and the like. However, it does seem to work well enough to remove inaudible details as needed.
* Build tools compile with SSE+SSE2/AVX+AVX2/FMA enabled by default. If the encoder crashes/doesn't work, change these flags in the ```Makefile```.

## Technical details
* Target bitrate: 16..256kbps+ (44.1kHz, M/S stereo)
    * No hard limits on playback rate or coding bitrate
* MDCT-based encoding (using sine window)
    * Encoding/decoding tools use N=4096 (2048 coefficients) with 37.5% overlap, but can use any sensible N=2<sup>n</sup> with any overlap (provided that the number of overlap samples is a multiple of 16)
* Extremely simple nybble-based syntax (no entropy-code lookups needed)
* Transient pre-echo reduction formula (more important at ultra-low bitrates)

## Authors
* **Ruben Nunez** - *Initial work* - [Aikku93](https://github.com/Aikku93)

## License
ulc-codec is released under the GPLv3 license. See the LICENSE file for full terms.

## Acknowledgements
* Huge thanks to Dennis K (DekuTree64) for listening to (and for playing programming-rubber-ducky to) all my rants and thought processes as I worked through understanding audio codec design
* Special thanks to [No!ze Freakz](https://soundcloud.com/user-462957379) for permission to use their music for demonstration purposes

## Gameboy Advance player

As a proof of concept of the decoding complexity, a Gameboy Advance demonstration may be found in the [ulcplayer-gba](https://github.com/Aikku93/ulcplayer-gba) repository.
