
cd src
echo "Making lpgenerator"
make
mv lpgenerator ../
cd ..
echo "Unzipping larger files"
unzip other.words.zip
cd trees
unzip alpha14.tree.zip
unzip coil14.tree.zip
echo "Building matcher"
cd ..
for i in 1 2 3 4 5 6 7 8; do
    cat segment$i >> matcher.txt
    rm segment$i
done

