for file in **/*.txt; do
    nfile=${file/ground_truth/withNativeDistances}
    mv "$file" "$nfile"
done
