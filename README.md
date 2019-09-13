# Pucr 

## Experimental Compression Playground

Looking out for some suitable compression for the lightweight VPN software [`n2n`](https://github.com/ntop/n2n), I happened to come across Pasi Ojala's [`pucrunch`](http://a1bert.kapsi.fi/Dev/pucrunch/). Not only that it is extremly assymetric and thus lightweight in terms of decompression, documentation also is highly explainatory giving full particulars in a wonderfully comprehensable way.

Having enjoyed reading that extremly inspring [article](https://github.com/Logan007/Pucr/blob/master/pucrunch/README.md) (converted to markdown, see pucrunch [folder](https://github.com/Logan007/Pucr/tree/master/pucrunch)) more than just several times, I started off to create my own implementation `pucr` which uses some of the original code and also includes some optimizations.

## Use Case

`n2n` software provides virtual ethernet adapters that encryptedly tunnel traffic between participants (at the `edges` of the network) even through NATs using hole punching technqiues but also the help of a forwarding `supernode` if necessary. To reduce traffic load and also to performancewise save some encryption cost, optional compression was planned for (using `minilzo`).

Today, most `edge` nodes presumably are desktop-like computers with a heap of CPU horse power that easily can afford some compression of small sized ethernet packets which usually do not grow above 1492 bytes in size.

However, some tiny `edges` in the _Internet of Things_ might not allow CPU cycles for compression, those would send out uncompressed data. To make them at least accept compressed packets, it is crucial that decompression does not eat up too many of thier precious CPU cycles. That is the reason for taking a closer look on compression algorithms offering that certain asymmetry.

## Modifications vis-à-vis the Original

The following list of changes might be completed and explained textually more detailed at a later point in time:

- An entry to the RLE table only saves bits if it gets used _more_ than once (on positions 1 to 3) or twice (on positions 4 to 15), respectively. With a view to relatively small packet sizes, this criterion becomes relevant and additionally is implemented.

- For each packet, the optimal bit-size of LZ2-offset is determined and used – in lieu of a fixated 8-bit size.

- General LZ offsets (for longer-than-two byte matches) are output twofold: The LSBs are stored in regular binary 2<sup>n</sup> coding whereas the MSBs are encoded using the variable length Elias Gamma coding (with inverted prefix). An optimization run to determine the optimal number of plainly encoded LSBs is performed.

- While determining the RLE and LZ costs, the longest match and _all_ shorter matches are checked.

- However, a `fast_lane` parameter switches off some of the mentioned optimizations making the program use less loops or just use some sane default values. ~~_Fun fact:_ The fastest `fast_lane` of `3` that just uses some default values results in even better compression for `ivanova.bin` (9620 bytes without header) than all slower modes! Maybe some potential left for optimizing the optimizing loops~~ _In fact, there was a bug in the lz offset lsb optimizer... The `ivanova.bin`-file now regularily gets compressed to 9605 (without header) bytes and to 9620 for all the fast lanes._

- The header contains parameters for the mentioned optimizations such as the _number of LSBs or LZ offset_ etc. Fields are used bitwise, e.g. 4 bits only for the _number of ESCape bits_. It does not require neither CBM-specific parts nor the integrated decruncher.

## Thoughts and Status

This all is _work in progress_! The code still is extremly polluted with `fprintf`s to `stderr` and other things. It has nearly no error checking and therefore is sensitive to malformed data. Also, the `fast_lane` is still hard-coded in  `main`. It can be found as the last parameter of `pucrunch_256_encode` where `0` is slowest, and `3` should be the fastest, `1` and `2` something in between – check it out.

One possible string matching speedup that takes advantage of RLE is still lacking, will follow.

So far, I haven't taken advantage of `max_gamma` yet, this definitely is a TODO.

So far, compression of 1492 byte sized packets is quickly done on my 8 year old i7-2860QM CPU. However, I will try to look deeper in how `pthread` could be of additional help here.

With a view to the presumed usecase, the chosen data types limit data size to `64K - 1` bytes (or so). This may be broadended in future versions.

The code should compile straight away by `gcc pucr.c -o pucr` to an executable called `pucr` – flawlessly in current Arch Linux. It takes input from `stdin`, compresses, decompresses and outputs to `stdout` while stats are output to `stderr`. One possible healthy call could look as follows:

``./pucr < ivanova.bin > ivanova.bin.pu.upu``

Hopefully, the two files are identical then. So far, a following `diff ivanova.bin ivanova.bin.pu.upu` has never complained yet.

Any hint or support is welcome – just leave an _Issue_!
