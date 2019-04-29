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
* Aside from a somewhat arbitrary high-frequency boost in the frequency selection criteria (and larger quantizer bands in the high-frequency range), there is no psychoacoustic optimizations, which might be resulting in lower quality for a given rate
    * The reason this model works at all is that low frequencies mask more neighbouring bands than high frequencies do (most audio contains a lot of low-frequency content relative to the high frequencies). However, this breaks down when there isn't much low-frequency masking going on, combined with with decent amount of high-frequency activity; this can lead to completely losing important low-frequency bands. For that reason, a suitable psychoacoustic model is currently being researched.
* Syntax is flexible enough to cause buffer overflows
* No block synchronization (if an encoded file is damaged, there is no way to detect where the next block lies)
    * It should be possible to prepend each block with the two-byte nybble sequence ```0h,0h,0h,0h```. Such a sequence could only happen at the start of a block (set quantizer to 2<sup>0</sup>, followed by three zero coefficients) and never in any other place (as four zero coefficients would be coded as ```8h,1h```), avoiding false-positives.

## Technical details
* Target bitrate: 16..256kbps+ (44.1kHz, M/S stereo)
    * No hard limits on playback rate or coding bitrate
* MDCT-based encoding (using sine window)
    * Encoding/decoding tools use N=4096 (2048 coefficients) with 37.5% overlap, but can use any sensible N=2<sup>n</sup> with any overlap (provided that the number of overlap samples is a multiple of 16)
        * Due to the extremely simplified quantization model combined with an improper psychoacoustic model, a large transform size is almost essential to avoid excessive quality degradation at low bitrates
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

As a proof of concept of the decoding complexity, a Gameboy Advance demonstration may be found in the ulcplayer-gba folder. CPU usage is around 65% for 32768Hz @ 128kbps (M/S stereo). Note that this is entirely a proof of concept; decode time for N=4096 (default for encoding tools) is 2-3 frames, so usage in real applications would need some form of threading to avoid excessive lag.

To use this player, you must:
* Provide your own ```SoundData.ulc``` in the ```source/res``` folder
* Modify the ```PATH``` variable in the ```Makefile``` to point to your build tools
* Compile with a suitable ARM assembler+linker (wholly written in assembly; no compiler needed)

The player supports both mono and stereo files (requires a rebuild; mono/stereo(+M/S) toggle found in ```source/ulc/ulc_Specs.inc```). M/S stereo isn't the best quality (as the transform is performed on clipped 8-bit samples), but this doesn't appear to cause noticeable quality degradation.

### Pre-built tracks

[No!ze Freakz - Freedom (64kbps @ 32.768kHz, M/S stereo)](ulcplayer-gba/No!ze%20Freakz%20-%20Freedom%20(64k).gba)
