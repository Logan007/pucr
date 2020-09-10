# Implicit Huffman Trees

As for entropy encoding, Huffmann would come in handy but making the receiver need to reconstruct the tree. So how could we avoid transmitting the tree itself?

Luckily, Move-to-Front transform usually outputs small values representing the locally most common values. In understanding that those should be represented by the shortest symbols, we could rearrange the regular `[$00..$FF]` byte alphabet as follows:

                        /\
                       /  \
                   0  /    \  1
                     /      \
                    /        \
               [$00..$3F]    /\
                            /  \
                        0  /    \  1
                          /      \
                         /        \
                   [$40..$7F]  [$80..$FF]

The first quarter of 6-bit symbols get a `0` prefix making them 7 bit long. The second quarter get a `10` prefix applied making those 6-bit symbols 8 bit long. Finally, the last half of – hopefully less often – 7-bit symbols get the `11` prefix applied, making those 9 bit long.

If Move-to-Front is extremely successful, it might be promising to recursively apply the same scheme to the first quarter, i.e. the leftmost branch (only!), leaving untouched the right side of the tree:

                        /           \
                       /             \
                   0  /               \  1
                     /                 \
                    /                   \
                   /\                   /\
                  /  \                 /  \
              0  /    \ 1          0  /    \  1
                /      \             /      \
               /        \           /        \
          [$00..$0F]    /\    [$40..$7F]  [$80..$FF]
                       /  \
                    0 /    \ 1
                     /      \
                    /        \
              [$10..$1F]  [$20..$3F]

In this case, the additional level would result in the following symbol lengths:

    00  + [$00..$0F](4 bits) = 6 bits
    010 + [$10..$1F](4 bits) = 7 bits
    011 + [$20..$3F](5 bits) = 8 bits

This could be continued up to level 4, in which case only the symbol `$00` would be encoded by 4 consecutive `0` bits and the other symbols longer respectively.

This system might not optimally relflect the actual entropy. On the other hand, to describe the tree, we only need to transmit the _level_ used. With a view to small packets to compress, maybe this already is quite sufficent despite non-optimality.

# Decoding

A directly deduced recursive way of an implementation:

         uint32_t get_huffed_value_rec (uint32_t no_of_bits, uint32_t rec_lvl) {

             if ( get_bits (1) == 0 )
                     if ( rec_lvl == 1 ) {
                             return (get_bits (no_of_bits - 2));
                     } else {
                             return (get_huffed_value_rec (no_of_bits - 2, rec_lvl - 1));
                     }
             else if ( get_bits(1) == 0 )
                     return (get_bits (no_of_bits - 2) | (1 << (no_of_bits - 2)));
             else
                     return (get_bits (no_of_bits - 1) | (1 << (no_of_bits - 1)));
         }

Note that the recursive call happens for a `0` bit prefix only. In other cases, no recursive calls happen. Also, recursion comes to an end as soon as `rec_lvl` gets too small. `no_of_bits` is the number of unencoded symbol bits – 8 for characters.

To avoid recursive function calls, we could linearize that function:

         uint32_t get_huffed_value_lin (uint32_t no_of_bits, uint32_t rec_lvl) {

             while ( rec_lvl > 0) {
                 if ( get_bits (1) == 0 )
                         if ( rec_lvl == 1 ) {
                                 return (get_bits (no_of_bits - 2));
                         } else {
                                 no_of_bits -= 2;
                                 rec_lvl--;
                                 continue;
                         }
                 else if ( get_bits(1) == 0 )
                         return (get_bits (no_of_bits - 2) | (1 << (no_of_bits - 2)));
                 else
                         return (get_bits (no_of_bits - 1) | (1 << (no_of_bits - 1)));
             }
             // rec_lvl == 0
             return (get_bits (no_of_bits));
         }

Recursion is replaced by a loop in which the parameters are adjusted the same way as in the recursive call. Also, we take care of the counted down `rec_lvl == 0` outside the loop.

While shifting towards a more branchless bit-reader design, I came up with the following well working version:

         uint32_t inline up_get_huffed_value (uint32_t no_of_bits, int32_t rec_lvl) {

             uint32_t input;

             uint32_t *buffer;
             buf = &data[bit_count >>3];
             inp = *buffer;

             input = be32toh (input);
             input <<= bit_count & 7;

             int32_t n = __builtin_clz(input);
             n = min (n, rec_lvl); // max. number of consecutive '0'-bits to read

             uint32_t m = gt (rec_lvl, n); // '1' if rec_lvl != n   '0' if rec_lvl == n
                                           // i.e. is there a 'right' subtree to take care of?

             uint32_t r = (no_of_bits - (n << 1)); // remaining bits to read... each leading '0' read
                                                   // here is good for '00' in the original prefix
             r -= (m << 1);                        // -0 if rec_lvl    or -2 for subtree indication otherwise

             input = rotl32 (input, n + (m << 1)); // n plus 0 or 2 bits more

             r += (m & input); // '0' for left sub-tree (10) or '1' for right sub-tree (11)

             bit_count += n + (m << 1) + r;

             input = rotl32 (input, r);
             input &= (1 << r) - 1;
             input |= m << r;

             return (input);
         }

`gt` and `min` denote implementations of _greater than_ as well as _minimum_ which rely on [sign-extension tricks](https://hbfs.wordpress.com/2008/08/05/branchless-equivalents-of-simple-functions/) using arithmetic right shift. Also, leftwise bit rotation [`rotl`](https://fgiesen.wordpress.com/2018/02/19/reading-bits-in-far-too-many-ways-part-1/) is used to the the _value of interest_ in the LSBs. 

It pays off having re-applied this scheme to the leftmost branch only as we just need to count the leading zeros using the corresponding intrinsic and thus finally allowing a branchless design. Will SSE allow for four parallel streams? To be tested somewhen…
