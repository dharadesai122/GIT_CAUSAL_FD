path="/home/mschaudh/CausalInference_Experiments"
alpha='0.05'
region=""

for dir in CV_Run_*; do
	echo $dir
	cd $dir
	for d in CV_[1-8]; do
		echo $d
		cd $d
		mkdir Causal_Graphs
		for f in ClimData*.txt; do
			echo $f
		
			IFS='.' read -ra file <<< "$f"
		        i="${file}" 
        	        echo $i
        			tetrad_cmd="java -jar $path/tetradcmd-5.1.0-10.jar -data $f -datatype continuous -algorithm pc.stable -depth -1 -significance $alpha -graphxml $path/Causal_Graphs/${i}_${alpha}.xml"
				$tetrad_cmd > file.out
		done
		cd ..
	done
	cd ..
done