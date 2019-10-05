# Pucr 

## Experimental Compression Playground

Looking out for some suitable compression for the lightweight VPN software [`n2n`](https://github.com/ntop/n2n), I happened to come across Pasi Ojala's [`pucrunch`](http://a1bert.kapsi.fi/Dev/pucrunch/). Not only that it is extremly assymetric and thus lightweight in terms of decompression, documentation also is highly explainatory giving full particulars in a wonderfully comprehensable way.

Having enjoyed reading that extremly inspring [article](https://github.com/Logan007/Pucr/blob/master/pucrunch/README.md) (converted to markdown, see pucrunch [folder](https://github.com/Logan007/Pucr/tree/master/pucrunch)) more than just several times, I started off to create my own implementation `pucr` which uses some of the original code and also includes some optimizations.

## Use Case

`n2n` software provides virtual ethernet adapters that encryptedly tunnel traffic between participants (at the `edges` of the network) even through NATs using hole punching technqiues but also the help of a forwarding `supernode` if necessary. To reduce network traffic and also to performancewise save some encryption cost, optional compression was planned for (using `minilzo`).

Today, most `edge` nodes presumably are desktop-like computers with a heap of CPU horse power that easily can afford some compression of small sized ethernet packets which usually do not grow above 1492 bytes in size.

However, some tiny `edges` in the _Internet of Things_ might not allow CPU cycles for compression, those would send out uncompressed data. To make them at least accept compressed packets, it is crucial that decompression does not eat up too many of their precious CPU cycles. That is the reason for taking a closer look on `pucrunch` as one of the compression algorithms offering that certain asymmetry.

As soon as `pucr` reaches a usable state, I will have it added to my fork of `n2n`.

## Modifications vis-à-vis the Original

The following list of changes might be completed and explained textually more detailed at a later point in time:

- An entry to the RLE table only saves bits if it gets used _more_ than once (on positions 1 to 3) or twice (on positions 4 to 15), respectively. With a view to relatively small packet sizes, this criterion becomes relevant and additionally is implemented.

- For each packet, the optimal bit-size of LZ2-offset is determined and used – in lieu of a fixed 8-bit size.

- General LZ offsets (for longer-than-two byte matches) are output twofold: The LSBs are stored in regular binary 2<sup>n</sup> coding whereas the exceeding MSBs are encoded using the variable length Elias Gamma coding (with inverted prefix). An optimization run to determine the optimal number _n_ of plainly encoded LSBs is performed per packet.

- While determining the RLE and LZ cost, the longest match and _all_ shorter matches are checked.

- All LITerals of the created graph – and only the LITerals – get a Move-to-Front encoding applied to hopefully minimize the use of _ESCaped_ LITerals by assembling most of the LITerals' values in the lower range. Also, it could help to keep the number of escape bits low and this shortening all other regular ESCape sequences.

- However, a `fast_lane` parameter switches off some of the mentioned optimizations making the program use less loops or just use some sane default values. 

- The header contains parameters for the mentioned optimizations such as the _number of LSBs_ or _LZ offset_ etc. Fields are used bitwise, e.g. 4 bits only for the _number of ESCape bits_. It does not require neither CBM-specific parts nor the integrated decruncher.

## Status, To-Do and Thoughts

This all is _work in progress_! The code still is extremly polluted with `fprintf`s to `stderr` and other things. It has nearly no error checking and therefore is sensitive to malformed data. Also, the `fast_lane` is still hard-coded in  `main`. It can be found as the last parameter of `pucrunch_256_encode` where `0` is slowest, and `3` should be the fastest, `1` and `2` something in between – check it out.

The `ivanova.bin`-file now regularily gets compressed to headerless ~~9567~~ __9320__ bytes (using Move-to-Front) and a few more for the fast lanes. By the way, the header would add 20 bytes in this case.

One possible string matching speedup that takes advantage of RLE is still lacking, will follow.

So far, no advantage of `max_gamma` is taken yet, this definitely is a todo.

Another optimization step from the original is still missing: After determining the best path through the graph (and thus the best LZ-mach lenghts which might be shorter than max), `pucr` should search for a maybe closer matches of that shorter length which could save a few more bits for offset. A quick proof of concept revealed _minus 24 bytes_ for the `ivanova.bin` file. Coming soon.

The Move-to-Front encoding is quite fast and works extremly well for most part of the available, limited test set. However, as it is context-dependant, there might be some cases which achieve better compression skipping Move-to-Front. Thus, a natural action item would be to let `pucr` figure out whether it is better to take advantage of it – or not. It is more an educated guess that even a CBM could easiliy perform Move-to-Front decoding for LITerals: It _temporarily_ requires _only_ 256 bytes of RAM for the alphabet and slightly more decompression code.

Having Move-to-Front encoding in place, a following huffman encoding step on the LITerals might be beneficial. Some early and dirty-coded tries look promising at least for `ivanova.bin` (maybe minus another 1908 _bits_ already including the tree) – not so much for the shorter ethernet packet sized files. To be refined and implemented until the end of the year.

As of now, compression of 1492 byte sized packets is quickly done on an old i7-2860QM CPU. However, I will try to look deeper in how `pthread` could be of additional help here.

With a view to the presumed usecase, the chosen data types limit data size to `64K - 1` bytes (or so). This may be broadended in future versions.

The code should compile straight away by `gcc pucr.c -o pucr` to an executable called `pucr` – flawlessly in current Arch Linux. It takes input from `stdin`, compresses, decompresses and outputs to `stdout` while stats are output to `stderr`. One possible healthy call could look as follows:

``./pucr < ivanova.bin > ivanova.bin.pu.upu``

Hopefully, the two files are identical then. So far, a following `diff ivanova.bin ivanova.bin.pu.upu` has never complained yet.

Any hint or support is welcome – just leave an _Issue_!
