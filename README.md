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
* Syntax is flexible enough to cause buffer overflows
* No block synchronization (if an encoded file is damaged, there is no way to detect where the next block lies)

## Technical details
* Target bitrate: 16..256kbps+ (44.1kHz, M/S stereo)
    * No hard limits on playback rate or coding bitrate
* MDCT-based encoding (50% overlap, using sine window)
    * Encoding/decoding tools use N=4096 (2048 coefficients), but can use any sensible 2^n
        * Due to extremely simplified quantization model, a large transform size is almost essential to avoid excessive quality degradation at low bitrates
* Extremely simple nybble-based syntax (no entropy-code lookups needed)
* Transient pre-echo reduction formula (more important at ultra-low bitrates)

## Authors
* **Ruben Nunez** - *Initial work* - [Aikku93](https://github.com/Aikku93)

## License
ulc-codec is released under the GPLv3 license. See the LICENSE file for full terms.

## Acknowledgements
* Huge thanks to Dennis K (DekuTree64) for listening to (and for playing programming-rubber-ducky to) all my rants and thought processes as I worked through understanding audio codec design
