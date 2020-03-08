### Overview
***

Each block is coded sequentially, with each channel also coded sequentially. For example:

```Block0[Chan0,Chan1], Block1[Chan0,Chan1]...```

The encoding tools align each block (that is, including all channels) to byte boundaries, but this is not strictly necessary; any alignment will do.

These blocks code MDCT coefficients that have been pre-normalized. So to decode, one must simply apply an unnormalized IMDCT to return the decoded PCM stream.

NB: The encoder is expected to handle all scaling, such that the inverse transform needs no scaling whatsoever (not even MDCT normalization).

### Syntax
***

| Nybble sequence    | Explanation        | Effect                                       |
| ------------------ | ------------------ | -------------------------------------------- |
| ```-7h..+7h```     | Normal coefficient | Insert ```Coef[n++] = Sign[Nybble]*Nybble^2 * Quantizer``` |
| ```8h,0h,0h..Dh``` | Quantizer change   | Set ```Quantizer = 2^-(4+X)```               |
| ```8h,0h,Eh,Xh```  | Quantizer change   | Set ```Quantizer = 2^-(4+14+X)```            |
| ```8h,0h,Fh```     | Stop               | Stop reading coefficients; fill rest with 0  |
| ```8h,1h..Bh```    | Zeros run (short)  | Insert zeros                                 |
| ```8h,Ch..Fh,Xh``` | Zeros run (long)   | Insert zeros                                 |

#### ```-7h..+7h```: Normal coefficient

This is a quantized coefficient. To unpack, take the (signed) squared value, scale by the current quantizer, and move to the next coefficient.

If we've received as many coefficients as there are in a transform block (eg. 2048 coefficients), this channel's coefficients are finished, and we can move on to inverse transforming.

#### ```8h,0h,0h..Dh```, ```8h,0h,Eh,Xh```: Quantizer change

This modifies the current quantizer.

When the third nybble is between ```0h..Dh```, this signals a normal quantizer change; at average audio levels, these are the most common quantization scalers.

```Quantizer = 2^-(4 + [third nybble])```

(Note the bias of 4: This is because no coefficient can be larger than 1.0, and the possibly (non-zero) quantized coefficients are ```{1,4,9,16,25,36,49}```, the log2 of the mean of which is approximately 4. If the bias wasn't used, there would be quantizer values that would never be used by the encoding tools, representing loss of efficiency)

When the nybble is ```Eh```, however, this signals that the quantizer needs to be smaller still, and so is followed by an extra nybble to afford a so-called 'extended-precision quantizer', allowing a maximum quantization scale of 2^-(4 + 14 + 15).

```Quantizer = 2^-(4 + 14 + [fourth nybble])```

This can be considered to result in a floating-point value containing a sign bit, 3 bits of mantissa, and circa 5 bits of exponent, with the minimum possible coefficient value being ±1.0×2^-33 (approximately 0.00000000012).

Each channel begins with a nybble specifying the initial quantizer, akin to a silent ```8h,0h``` at the start (extended-precision quantizers are also allowed for this first quantizer band).

#### ```8h,0h,Fh```: Stop

This signals that the remaining coefficients for this channel should be zeroed out.

A channel's block may begin with ```[8h,0h,]Fh```; this means that the channel is silent.

#### ```8h,1h..Bh```: Zeros run (short)

This inserts coefficients with a value of 0. As these are fairly common after quantization, we use a compact format to code runs of them.

To increase density, these runs are always in multiples of 2, starting at 4 (eg. 4, 6, 8, 10...24). This results in an average code length of 2.5 nybbles per zero run, if we assume that there's a 50:50 chance of even:odd runs (8h,Xh = 2 nybbles, 8h,Xh,0h: 3 nybbles).

To unpack: ```n = [second nybble] * 2 + 2```

#### ```8h,Ch..Fh,Xh```: Zeros run (long)

Further refinement on the above. This in particular is more common at lower bitrates, where quantization becomes too aggressive and/or we must cull large swathes of frequency bands in order to meet the target bitrate.

Like the above short run, this also runs in multiples of 2 (average length of 3.5 nybbles). However, we supply an extra nybble to extend the range much further.

Because a short run can already code up to 24 zeros at once, a long run codes 26,28,30...152 zeros at once.

To unpack: ```n = (([second nybble] - Ch)<<4 | [third nybble]) * 2 + 26```

### Inverse transform process
***

Inverse transform takes the above-decoded coefficients and applies the IMDCT to give the output stream.

#### IMDCT

The IMDCT used in this codec follows the normal IMDCT formula found on any maths book, but without any scaling factor in front (as it is moved to the MDCT in the encoder), ie.

    y[n] = Sum[X[k]*Cos[(n+1/2 + N/2)(k+1/2)*Pi/N], {k,0,N-1}]
    
The windowing function used is a sine window. To use a different window, the MDCT and IMDCT functions of the source code must be changed.

Note that when using less than 50% overlap, the MDCT/IMDCT transforms can be slightly improved in efficiency by offsetting the basis functions by -1 in ```n``` and ```k``` (ie. ```n-1/2``` and ```k-1/2```), as this avoids having to negate coefficients in the non-overlap region (for the overlapping region, this simply involves changing the quadrature oscillator signs). Note that while this results in negated MDCT coefficients as output, the phase-inversion is corrected in the IMDCT stage so long as both MDCT and IMDCT use the same shifted basis.
