### Overview
***

Each block is coded sequentially, with each channel also coded sequentially. For example:

```Block0[Chan0,Chan1], Block1[Chan0,Chan1]...```

The encoding tools align each block (that is, including all channels) to byte boundaries, but this is not strictly necessary; any alignment will do.

These blocks code MDCT coefficients that have been pre-filtered with a pre-echo reduction formula. Therefore, to decode, we simply apply the inverse of said formula and then perform an IMDCT.

NB: The encoder is expected to handle all scaling, such that the inverse transform needs no scaling whatsoever.

### Syntax
***

| Nybble sequence    | Explanation        | Effect                                      |
| ------------------ | ------------------ | ------------------------------------------- |
| ```-7h..+7h```     | Normal coefficient | Insert ```Coef[n++] = Nybble * Quantizer``` |
| ```8h,0h,0h..Eh``` | Quantizer change   | Set ```Quantizer = 2^X```                   |
| ```8h,0h,Fh```     | Stop               | Stop reading coefficients; fill rest with 0 |
| ```8h,1h..Bh```    | Zeros run (short)  | Insert zeros                                |
| ```8h,Ch..Fh,Xh``` | Zeros run (long)   | Insert zeros                                |

#### ```-7h..+7h```: Normal coefficient

This is a quantized coefficient. To unpack, just scale by the current quantizer, and move to the next coefficient.

If we've received as many coefficients as there are in a transform block (eg. 2048 coefficients), this channel's coefficients are finished, and we can move on to inverse transforming.

#### ```8h,0h,0h..Eh```: Quantizer change

This modifies the current quantizer. Simple enough.

```Quantizer = 2^[third nybble]```

#### ```8h,0h,Fh```: Stop

This signals that the remaining coefficients for this channel should be zeroed out.

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

Inverse transform takes the above-decoded coefficients and applies the inverse of a pre-echo reduction formula, followed by IMDCT.

#### Pre-echo reduction formula

The pre-echo reduction formula used in this codec is a basic sum/difference transform.

The transform (and its inverse, due to its involutive nature) is as follows:

    X[2*n+0] = x[2*n+0] + x[2*n+1]
    X[2*N+1] = x[2*n+0] - x[2*n+1]
    
#### IMDCT

The IMDCT used in this codec follows the normal IMDCT formula found on any maths book, but without any scaling factor in front (as it is moved to the MDCT in the encoder), ie.

    y[n] = Sum[X[k]*Cos[(n+1/2 + N/2)(k+1/2)*Pi/N], {k,0,N-1}]
    
The windowing function used is a sine window with 50% overlap, resulting in what is known as the modulated lapped transform (MLT).
