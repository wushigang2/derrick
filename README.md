# derrick
derrick is a tool for DNA storage encode/decode.

# Installation
```sh
git clone https://github.com/wushigang2/derrick.git
cd derrick
make
```

# run derrick
```sh
derrick
```
```sh
commands:
 encode      encode a file.
 decode      decode a file.
```

# encode a file
```sh
derrick encode
```
```sh
usage: derrick encode [options] <input file>
 -i <string> pi file, [NULL]
 -n <int>    n of rs(n,k), n = 2^N - 1, N = 8/10/12/14/16, [255]
 -k <int>    k of rs(n,k), [235]
 -s <int>    number of symbol(or call rs) per segment(or call block), [62]
 -v          verbose.
```

# decode a file
```sh
derrick decode
```
```sh
usage: derrick decode [options] <input file>
 -i <string> pi file, [NULL]
 -n <int>    n of rs(n,k), n = 2^N - 1, N = 8/10/12/14/16, [255]
 -k <int>    k of rs(n,k), [235]
 -s <int>    number of symbol(or call rs) per setment(or call block), [62]
 -M <int>    score for match, [2]
 -X <int>    penalty for mismatch, [-6]
 -O <int>    penalty for gap open, [-3]
 -E <int>    penalty for gap extension, [-2]
 -c <int>    max change number of rs soft decision, [0]
 -d <int>    max delete number of rs soft decision, [0]
 -m <string> sequencing mode: illumina/pacbio/nanopore, [pacbio/nanopore]
 -j <int>    jump mode of collision, [0]
             0: jump to last exceed.
             2: jump to first continuous exceed.
 -u <int>    search upper range of candidate, [255]
 -l <int>    search lower range of candidate, [32]
 -r <int>    search raise range of candidate, [0]
 -t <int>    how many seconds a block timeout, [3600]
 -R <string> ref file, [NULL]
 -b <int>    beg block, [0]
 -e <int>    end block, [0]
 -f <int>    shift or not, [1]
 -v          verbose.
```

# Contact
Shigang Wu wushigang@caas.cn
