#!/bin/bash

if [ -z $1 ]; then
	echo "Pass version tag to check"
	exit 127
fi

trap ctrl_c INT

function ctrl_c() {
        echo "** Trapped CTRL-C"
        jobs -p | xargs kill
        exit 1
}

image=docker.io/avirshup/mst:workflows-$1

function runsteps(){
	cmd=$1
	testdir=$3
	input=$2
	finalizing=$4

	if [ -d $testdir ]; then
      	echo "Cleaning old test directory $dir"
      	rm -r $testdir
    fi
    mkdir $testdir

	mkdir $testdir/preprocess
	torun="$cmd $input --preprocess --outputdir $testdir/preprocess" 
	echo
	echo "> $torun"
	$torun 1>$testdir/preprocess/out 2>$testdir/preprocess/err

	if [ $? -eq 0 ]; then
		echo -e "\nPreprocessing succesful: $testdir/preprocess"
	else
		echo -e "\nERROR for $testdir/preprocess:"
		head -n 40 $testdir/preprocess/err
		echo "FAILED: $testdir"
		return 1
	fi

	mkdir $testdir/finish
	torun="$cmd --restart $testdir/preprocess/workflow_state.dill --outputdir $testdir/finish $finalizing"
	echo "> $torun"
	$torun 1>$testdir/finish/out 2>$testdir/finish/err

	if [ $? -eq 0 ]; then
		echo -e "\nSUCCESS: $testdir"
	else
		echo -e "\nERROR for $testdir/finish:"
		head -n 40 $testdir/finish/err
		echo "FAILED: $testdir"
		return 1
	fi

	return 0
}


function runsteps_docker(){
	cmd=$1
	testdir=$3
	input=$2
	finalizing=$4

	if [ -d $testdir ]; then
      	echo "Cleaning old test directory $dir"
      	rm -r $testdir
    fi
    mkdir $testdir

    dockercmd="docker run -v /var/run/docker.sock:/var/run/docker.sock"
    absdir=`pwd`/$testdir

	mkdir $testdir/preprocess
	torun="$dockercmd -v $absdir/preprocess:/outputs $cmd $input --preprocess --outputdir /outputs" 
	echo
	echo "> $torun"
	$torun 1>$testdir/preprocess/out 2>$testdir/preprocess/err

	if [ $? -eq 0 ]; then
		echo -e "\nPreprocessing succesful: $testdir/preprocess"
	else
		echo -e "\nERROR for $testdir/preprocess:"
		head -n 40 $testdir/preprocess/err
		echo "FAILED: $testdir"
		return 1
	fi

	mkdir $testdir/finish
	torun="$dockercmd -v $absdir/finish:/outputs -v $absdir/preprocess:/inputs $cmd --restart /inputs/workflow_state.dill --outputdir /outputs $finalizing"
	echo "> $torun"
	$torun 1>$testdir/finish/out 2>$testdir/finish/err

	if [ $? -eq 0 ]; then
		echo -e "\nSUCCESS: $testdir"
	else
		echo -e "\nERROR for $testdir/finish:"
		head -n 40 $testdir/finish/err
		echo "FAILED: $testdir"
		return 1
	fi

	return 0
}


echo 'Launching local VDE test'
runsteps "chemworkflow vde --localdocker" \
         "'{\"input\":\"CC\"}'" \
         testdir/test_local_vde &

echo 'Launching local minimize test'
runsteps "chemworkflow minimize --localdocker" \
         "'{\"input\":\"3AID\"}'" \
         testdir/test_local_minimize &

echo 'Launching docker VDE test'
runsteps_docker "$image vde --localdocker" \
			"'{\"input\":\"CC\"}'" \
             testdir/test_docker_vde &

echo 'Launching docker minimize test'
runsteps_docker "$image minimize --localdocker" \
              "'{\"input\":\"3AID\"}'" \
              testdir/test_docker_minimize &

wait

echo $done