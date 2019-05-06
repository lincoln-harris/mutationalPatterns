##fileformat=VCFv4.2
#!/bin/bash

for f in *; do
	echo $f
	sed -i '1i##fileformat=VCFv4.2' $f
done