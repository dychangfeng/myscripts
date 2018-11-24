for i in $(ls ~/chroms)
    do  nohup python ~/myscripts/cruciform_DNA3.py -f ~/chroms/$i -p 2 >$i.loop6_7.txt &
done