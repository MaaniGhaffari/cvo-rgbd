for d in /data/rgbd_dataset/*/ ; do
    echo "$d"
    rm $d"assoc.txt"
    python3 associate.py $d"rgb.txt" $d"depth.txt" >> $d"assoc.txt"
done
