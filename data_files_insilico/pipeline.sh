# encode
perl moban.pl -o vol3.txt -i children.mp4 -i pgzf -i ruan_group.jpg -i wtdbg2 -i wtdbg-cns -i wtpoa-cns > wcl.txt
derrick encode -i pi.txt -n 255 -k 211 -s 62 vol3.txt > ref.211.fa
gzip ref.211.fa
# decode
gunzip ref.211.fa.gz
derrick decode -i pi.txt -n 255 -k 211 -s 62 ref.211.fa > vol3.txt
perl mom.pl -i1 wcl.txt -i2 vol3.txt -o1 children.mp4 -o2 pgzf -o3 ruan_group.jpg -o4 wtdbg2 -o5 wtdbg-cns -o6 wtpoa-cns
