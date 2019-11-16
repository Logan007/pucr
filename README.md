# Pucr 

## Experimental Compression Playground

Looking out for some suitable compression for the lightweight VPN software [`n2n`](https://github.com/ntop/n2n), I happened to come across Pasi Ojala’s [`pucrunch`](http://a1bert.kapsi.fi/Dev/pucrunch/). Not only that it is extremly assymetric and thus lightweight in terms of decompression, documentation also is highly explainatory giving full particulars in a wonderfully comprehensable way.

Having enjoyed reading that extremly inspring [article](https://github.com/Logan007/Pucr/blob/master/pucrunch/README.md) (converted to markdown, see pucrunch [folder](https://github.com/Logan007/Pucr/tree/master/pucrunch)) more than just several times, I started off to create my own implementation `pucr` which uses some of the original code and also includes some optimizations.

## Use Case

`n2n` software provides virtual ethernet adapters that encryptedly tunnel traffic between participants (at the `edges` of the network) even through NATs using hole punching technqiues but also the help of a forwarding `supernode` if necessary. To reduce network traffic and also to performancewise save some encryption cost, optional compression was planned for (using `minilzo`).

Today, most `edge` nodes presumably are desktop-like computers with a heap of CPU horse power that easily can afford some compression of small sized ethernet packets which usually do not grow above 1492 bytes in size.

However, some tiny `edges` in the _Internet of Things_ might not allow CPU cycles for compression, those would send out uncompressed data. To make them at least accept compressed packets, it is crucial that decompression does not eat up too many of their precious CPU cycles. That is the reason for taking a closer look on `pucrunch` as one of the compression algorithms offering that certain asymmetry.

As soon as `pucr` reaches a usable state, I will have it added to my fork of `n2n`.

## Modifications vis-à-vis the Original

The following list of changes might be completed and explained textually more detailed at a later point in time:

- An entry to the RLE character table only saves bits if it gets used _more_ than once (on positions 1 to 3) or twice (on positions 4 to 15), respectively. With a view to relatively small packet sizes, this criterion becomes relevant and additionally is implemented.

- For each packet, the optimal bit-size of LZ2-offset is determined and used – in lieu of a fixed 8-bit size.

- General LZ offsets (for longer-than-two byte matches) are output twofold: The LSBs are stored in regular binary 2<sup>n</sup> coding whereas the exceeding MSBs are encoded using the variable length Elias Gamma coding (with inverted prefix). An optimization run to determine the optimal number _n_ of plainly encoded LSBs is performed per packet.

- While determining the RLE and LZ cost, the longest match and _all_ shorter matches are checked.

- All LITerals of the created graph – and only the LITerals – get a Move-to-Front encoding applied to hopefully minimize the use of _ESCaped_ LITerals by assembling most of the LITerals’ values in the lower range. Also, it could help to keep the number of escape bits low and thus shortening all other regular ESCape sequences. However, as it is context-dependant, there are some cases which achieve better compression skipping Move-to-Front. Thus, `pucr` figures out whether it is better to take advantage of it – or not.

- Inspired by the Wikipedi [article](https://en.wikipedia.org/wiki/Move-to-front_transform#Example), Move-to-Front received a _Second Line of Defence_. i.e. the index after which the LITerals are inserted (standard Move-to-Front always uses `0`) is parametrized. This is beneficial for repeating LITerals – those would get sorted below that certain index. That could be the case, if the same LITeral gets re-used at the next LITeral position but especially for non-ranked RLEs.

- Having Move-to-Front in place, a following Huffman encoding step on the LITerals is performed where beneficial. Not to interfere with the well working ESCape-system, only the trailing bits of the literals, i.e. the last `8 minus number_of_escape_bits` bits, are Huffman-encoded. Each distinct possible ESCape sequence, e.g. `00`, `01`, `10`, and `11` in case of two ESCape bits, gets their own Huffman-tree – but only, if it saves bits. First trials show that so far, only the `00...`-symbols seem to save enough to afford the costs of a tree; maybe that is another consequence of the preceding Move-to-Front. The tree bitwisely gets encoded to the header.

- However, a `fast_lane` parameter switches off some of the mentioned optimizations making the program use less loops or just use some sane default values. 

- The header contains parameters for the mentioned optimizations such as the _number of LSBs_ or _LZ offset_ etc. Fields are used bitwise, e.g. 4 bits only for the _number of ESCape bits_. It does not require neither CBM-specific parts nor the integrated decruncher.

## Status, To-Do and Thoughts

This all is _work in progress_! The code still is extremly polluted with `fprintf`s to `stderr` and other things. It has nearly no error checking and therefore is sensitive to malformed data. Also, the `fast_lane` is still hard-coded in  `main`. It can be found as the last parameter of `pucrunch_256_encode` where `0` is slowest, and `3` should be the fastest, `1` and `2` something in between – check it out.

The `ivanova.bin`-file now regularily gets compressed to __9245__ bytes including the header and a few more for the fast lanes.

One possible string matching speed-up that takes advantage of RLE is still lacking; it might follow as soon as a way is found how not to loose compression ratio. Maybe, this is a `fast_lane` candidate.

So far, only some advantage of `max_gamma` is taken yet. Is this is a high-priority todo? An implementation with manually set `max_gamma` (set to 7, see line 1155) for the LZ lengths showed only a tiny compression gain... Maybe better for RLE lengths?

Another finding while testing the use of LZ length history (encoding the historic lengths as the first LZ lengths, the non-historic ones becoming longer then), it always increased output file size for no matter what history lenght being used (1,2,3, or 4). Shall we try history for LZ offset? Hope is, that compression of repetitions of similar sequences, such as `12345678` and `1234X678` will benefit.

An LZ offset table similar to the already implemented RLE character table might be worth deliberating. It would require some changes in the code and especially the way LZ offsets are being encoded... not sure if it yields better comression. Maybe not to be tested anytime soon.

~~Talking about tables and regarding Elias-Gamma coded numbers as binary-tree-position (`0` bit for _branch to the left_, `1` bit for _brach to the right_), a ranking table for LITerals (to be spcific, their _reamainders_ following the possible ESCape bits) might lead to similar compression as the currently implemented Huffman coded LITerals. Drawback: In comparison to the Huffman tree, it probably is not the optimal encoding but on the other hand, we do not need to transmit a list of unranked LITerals (as we would need to in case of a Huffmann tree). A big advantage could be speed, as the ranking and the corresponding Elias-Gamma code implicitly already include the tree. The LITeral's treatment shall be aim of some more thoughts.~~ _This does not seems to work out too well. Maybe we need to try a more flat implicit tree?_

To overlap or not to overlap... Currently, LZ matches cannot overlap the pattern, `find_matches` and the graph optimizer just assume full matches to end before the pattern starts. This allows to encode LZ-offsets beginning with `0` for the position `pattern - match_length`. Overlapping matches would be beneficial only to represent unranked RLE longer than 2 or recurring patterns. A full implementation might eat up some of the speed-wins gained in `find_matches` which would require more flexibility , e.g. cannot presume equal values of `rle_count`. In addition, the `offset` pointing to the match needs more bits for encoding. First experiments do not look too promising. This point needs thoughts and will be reconsidered in conjunction with RLE.

The Move-to-Front encoding is quite fast and works extremly well for most part of the available, limited test set.

As of now, compression of 1492 byte sized packets is quickly done on an old i7-2860QM CPU. The roughly drafted code especially gains speed from `gcc`’s `-O3` option. However, plans are to look deeper in how `pthread` could be of additional help here.

With a view to the presumed usecase, the chosen data types limit data size to `64K - 1` bytes (or so). This may be broadended in future versions.

The code should compile straight away by `gcc -O3 pucr.c -o pucr` to an executable called `pucr` – flawlessly in current Arch Linux. It takes input from `stdin`, compresses, decompresses and outputs to `stdout` while stats are output to `stderr`. One possible healthy call could look as follows:

``./pucr < ivanova.bin > ivanova.bin.pu.upu``

Hopefully, the two files are identical then. Up to now, a following `diff ivanova.bin ivanova.bin.pu.upu` has never complained yet.

Any hint or support is welcome – just leave an _Issue_!
