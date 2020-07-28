### Overview
***

Each block is coded sequentially, with each channel also coded sequentially. For example:

```Block0[Chan0,Chan1], Block1[Chan0,Chan1]...```

The encoding tools align each block (that is, including all channels) to byte boundaries, but this is not strictly necessary; any alignment will do.

These blocks code MDCT coefficients that have been pre-normalized. So to decode, one must simply apply an unnormalized IMDCT to return the decoded PCM stream.

Each block always starts with a window control code before any coefficients are coded. This controls the window lengths and shapes, and is explained in the `Block header` section.

NB: The encoder is expected to handle all scaling, such that the inverse transform needs no scaling whatsoever (not even MDCT normalization).

### Block header
***

Each block starts with a nybble that encodes overlap scaling, with the high bit being a decimation (window switch) toggle.

| Nybble sequence    | Window 0 | Window 1 | Window 2 | Window 3 |
| ------------------ | -------- | -------- | -------- | -------- |
| ```0h..7h```       | N*       | -        | -        | -        |
| ```8h..Fh,0010b``` | N/2*     | N/2      | -        | -        |
| ```8h..Fh,0011b``` | N/2      | N/2*     | -        | -        |
| ```8h..Fh,0100b``` | N/4*     | N/4      | N/2      | -        |
| ```8h..Fh,0101b``` | N/4      | N/4*     | N/2      | -        |
| ```8h..Fh,0110b``` | N/2      | N/4*     | N/4      | -        |
| ```8h..Fh,0111b``` | N/2      | N/4      | N/4*     | -        |
| ```8h..Fh,1000b``` | N/8*     | N/8      | N/4      | N/2      |
| ```8h..Fh,1001b``` | N/8      | N/8*     | N/4      | N/2      |
| ```8h..Fh,1010b``` | N/4      | N/8*     | N/8      | N/2      |
| ```8h..Fh,1011b``` | N/4      | N/8      | N/8*     | N/2      |
| ```8h..Fh,1100b``` | N/2      | N/8*     | N/8      | N/4      |
| ```8h..Fh,1101b``` | N/2      | N/8      | N/8*     | N/4      |
| ```8h..Fh,1110b``` | N/2      | N/4      | N/8*     | N/8      |
| ```8h..Fh,1111b``` | N/2      | N/4      | N/8      | N/8*     |

Overlap scaling (bit0..2 of the first nybble) is applied to the \[sub]block denoted with an asterisk, and results in an overlap of ```50% * 2^-Scale``` for that \[sub]block (that is, ```SubBlockSize * 2^-Scale``` samples of overlap).

Note that the transient subblock index can be easily obtained using a population count (POPCNT) on bit1..3 of the second nybble (minus 1 for 0-based indexing).

### Block syntax
***

| Nybble sequence       | Explanation        | Effect                                       |
| --------------------- | ------------------ | -------------------------------------------- |
| ```-7h..+7h```        | Normal coefficient | Insert ```Coef[n++] = Sign[Nybble]*Nybble^2 * Quantizer``` |
| ```8h,0h,0h..Dh```    | Quantizer change   | Set ```Quantizer = 2^-(5+X)```               |
| ```8h,0h,Eh,0h..Ch``` | Quantizer change   | Set ```Quantizer = 2^-(5+14+X)```            |
| ```8h,0h,Fh```        | Stop               | Stop reading coefficients; fill rest with 0  |
| ```8h,1h..Eh```       | Zeros run (short)  | Insert zeros                                 |
| ```8h,Fh,Yh,Xh```     | Zeros run (long)   | Insert zeros                                 |

#### ```-7h..+7h```: Normal coefficient

This is a quantized coefficient. To unpack, take the signed squared value (ie. ```x*ABS(x)```), scale by the current quantizer, and move to the next coefficient.

If we've received as many coefficients as there are in a transform block (eg. 2048 coefficients), this channel's coefficients are finished, and we can move on to inverse transforming.

#### ```8h,0h,0h..Dh```, ```8h,0h,Eh,0h..Ch```: Quantizer change

This modifies the current quantizer.

When the third nybble is between ```0h..Dh```, this signals a normal quantizer change; at average audio levels, these are the most common quantization scalers.

```Quantizer = 2^-(5 + [third nybble])```

(Note the bias of 5)

When the nybble is ```Eh```, however, this signals that the quantizer needs to be smaller still, and so is followed by an extra nybble to afford a so-called 'extended-precision quantizer', allowing a maximum quantization scale of 2^-31.

```Quantizer = 2^-(5 + 14 + [fourth nybble])```

This can be considered to result in a floating-point value containing a sign bit, 3 bits of mantissa, and 5 bits of exponent, with the minimum possible coefficient value being ±1.0×2^-31.

Each channel begins with a nybble specifying the initial quantizer, akin to a silent ```8h,0h``` at the start (extended-precision quantizers are also allowed for this first quantizer band).

#### ```8h,0h,Fh```: Stop

This signals that the remaining coefficients for this channel should be zeroed out.

A channel's block may begin with ```[8h,0h,]Fh```; this means that the channel is silent.

#### ```8h,1h..Eh```: Zeros run (short)

This inserts coefficients with a value of 0. As these are fairly common after quantization, we use a compact format to code runs of them.

To unpack: ```n = [second nybble] + 2```. This allows 3..16 zeros to be coded with two nybbles.

#### ```8h,Fh,Yh,Xh```: Zeros run (long)

Further refinement on the above. This in particular is more common at lower bitrates, where quantization becomes too aggressive and/or we must cull large swathes of frequency bands in order to meet the target bitrate.

Because a short run can already code up to 16 zeros at once, a long run codes 17..272 zeros at once.

To unpack: ```n = ([second nybble]<<4 | [third nybble]) + 17```

### Inverse transform process
***

Inverse transform takes the above-decoded coefficients and applies the IMDCT to give the output stream.

#### Deinterleaving

Decimated subblocks (ie. when using window switching) have their coefficient interleaved/shuffled to obtain a more compact bitstream. The interleaving pattern groups coefficients of larger subblocks together so that they share the same frequency bands as that of the smallest subblock (for example, if the smallest subblock is N/8 and the largest is N/2, then 4 coefficients of the N/2 subblock are grouped as one), and then ordered by window position (eg. for A=N/2,B=N/4,C=N/8,D=N/8, the following coefficients are coded one after the other: 4 coefficients of A, 2 coefficients of B, 1 coefficient of C, 1 coefficient of D). Deinterleaving must simply undo this ordering to restore all of the subblocks, which may then proceed through IMDCT one after the other.

#### IMDCT

The IMDCT used in this codec follows the normal IMDCT formula found on any maths book, but without any scaling factor in front (as it is moved to the MDCT in the encoder) and with a shifted basis that results in phase inversion (as this results in programming optimizations), ie.

    y[n] =  Sum[X[k]*Cos[(n+1/2 + N/2 + N*2)(k+1/2)*Pi/N], {k,0,N-1}]
         = -Sum[X[k]*Cos[(n+1/2 + N/2      )(k+1/2)*Pi/N], {k,0,N-1}]
    
The windowing function used is a sine window, with the overlap amount based on the block size scaled by the nybble at the start of the block (see the `Overview` section). To use a different window, the MDCT and IMDCT functions of the source code must be modified to accomodate such (this is not too difficult, and only involves loading the sine/cosine values with appropriate data). Note that using a sine window allows reuse of the DCT coefficients table, whereas a different window cannot reuse these coefficients and so needs double the storage space.

When using window switching, it is important to note that a subblock's overlap may be larger than allowed by the last subblock. When this happens, the number of overlap samples must be clipped to the size of the previous, smaller subblock. This unfortunately results in an additional block delay for decoding (on top of the MDCT delay), as the encoder must have knowledge about the next \[sub]block to account for this.
